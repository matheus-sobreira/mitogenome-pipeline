#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# Pilot QC — Analisa uma amostra piloto de reads e recomenda max_reads ideal
#
# Métricas extraídas:
#   1. Porcentagem de bases ≥ Q30 (qualidade Phred)
#   2. Porcentagem de reads com adaptador Illumina
#   3. Comprimento médio dos reads
#   4. Fração mitocondrial MEDIDA por mapeamento (bwa mem + samtools depth)
#   5. Breadth de cobertura e divergência dos reads mapeados
#
# Cálculo:
#   mito_reads = (target_cov × genome_size) / (read_len × 2)
#   adjusted   = mito_reads × quality_factor × adapter_factor
#   total      = adjusted / mito_fraction × safety_margin (1.5×)
#
# ── Sobre a medição da fração mitocondrial (Fase 3) ──────────────────────────
# A versão anterior estimava a fração com `grep` de k-mers da seed. Aquele
# estimador tinha erro de unidade — computava COBERTURA e a tratava como
# CONTAGEM DE READS — e só buscava a fita direta. Na prática saturava: nas duas
# execuções reais (A. leari e A. hyacinthinus) a fração caiu no piso do clamp e
# a recomendação no teto, ou seja, nenhum dos dois números veio de uma medição.
# Ver DEC-04 no vault.
#
# A medição agora é por profundidade (DEC-05):
#
#     fração = (profundidade_média × genome_avg) / (read_len × total_reads)
#
# Profundidade, e não contagem de reads mapeados, porque numa referência curta
# (seed COX1 de ~1.5 kb com reads de 150 bp) as bordas dominam a contagem.
#
# A referência pode ser (DEC-06):
#   - um mitogenoma completo (--pilot_reference) → medição direta
#   - a própria seed COX1 → a fórmula acima ESCALA implicitamente, ao usar
#     genome_avg no lugar do comprimento da referência. Isso assume cobertura
#     uniforme ao longo do mitogenoma; a premissa é declarada no relatório.
#
# Uso: pilot_qc.sh <R1> <R2> <ref> <genome_range> <target_cov> [mito_frac]
#                  [ref_mode] [threads] [bwa_opts]
# Saídas: pilot_report.txt, recommended_reads.txt, read_length.txt
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

R1="$1"
R2="$2"
REF="$3"
GENOME_RANGE="$4"
TARGET_COV="${5:-500}"
MITO_FRAC_OVERRIDE="${6:-0}"
REF_MODE="${7:-seed}"          # 'seed' (COX1, escalado) | 'genome' (completo)
THREADS="${8:-2}"
BWA_OPTS="${9:-}"

# Teto de segurança — reportado quando aplicado (DEC-10). O PISO foi removido:
# era da mesma família do piso de fração, existia para proteger de um estimador
# quebrado e, com medição real, só forçava download desnecessário.
MAX_READS="${PILOT_MAX_READS_CAP:-25000000}"

# ── Constantes do montador, medidas em 13/08/2026 (ver DEC-13) ──────────────
# MEM_GB      : memória que o NOVOPlasty receberá (resolvida no main.nf, para
#               que o Pilot QC e o montador usem o MESMO valor)
# READS_PER_GB: 1.691.018 reads/GB, desvio de 0,005% em 5 pontos (2–8 GB).
#               ⚠️ medido com reads de 150 bp; a memória guarda BASES, então a
#               constante quase certamente escala com o comprimento do read.
# RETENTION   : o NOVOPlasty não usa 100% do que recebe mesmo cabendo na
#               memória — retenção observada de 82% a 92%. 0,85 é conservador.
# TRIM_LOSS   : perda no Trim Galore (medida: 7,6% em A. leari)
MEM_GB="${NOVOPLASTY_MEM_GB:-7}"
READS_PER_GB="${NOVOPLASTY_READS_PER_GB:-1691018}"
RETENTION="${NOVOPLASTY_RETENTION:-0.85}"
TRIM_LOSS="${TRIM_LOSS:-0.10}"

# ── Parse genome range ──────────────────────────────────────────────────────
GENOME_MIN=$(echo "$GENOME_RANGE" | cut -d'-' -f1)
GENOME_MAX=$(echo "$GENOME_RANGE" | cut -d'-' -f2)
GENOME_AVG=$(( (GENOME_MIN + GENOME_MAX) / 2 ))

echo "[PILOT] Analisando amostra piloto..." >&2

# ── Métricas básicas ────────────────────────────────────────────────────────
TOTAL_READS=$(awk 'END{printf "%d", NR/4}' "$R1")
AVG_READ_LEN=$(awk 'NR%4==2 {total+=length($0); n++} END{printf "%d", total/n}' "$R1")

# Reads somando os dois mates — é esta a base do cálculo de fração, porque a
# profundidade vem de um BAM que contém R1 e R2.
TOTAL_READS_R2=$(awk 'END{printf "%d", NR/4}' "$R2")
TOTAL_READS_BOTH=$(( TOTAL_READS + TOTAL_READS_R2 ))

# ── Q30 (Phred+33: '?' = ASCII 63 = Q30) ───────────────────────────────────
Q30_PCT=$(awk 'NR%4==0 {
    for(i=1; i<=length($0); i++){
        c = substr($0,i,1); total++
        if (c >= "?") q30++
    }
} END {
    printf "%.1f", 100*q30/total
}' "$R1")

# ── Adaptador Illumina TruSeq (13-mer universal) ───────────────────────────
# Busca APENAS em linhas de sequência (NR%4==2) para evitar falsos positivos
# em linhas de qualidade (códigos ASCII coincidem com bases).
ADAPTER_KMER="AGATCGGAAGAGC"
ADAPTER_HITS=$(awk -v kmer="$ADAPTER_KMER" 'NR%4==2 && index($0,kmer){n++} END{print n+0}' "$R1")
ADAPTER_PCT=$(awk -v h="$ADAPTER_HITS" -v t="$TOTAL_READS" 'BEGIN{printf "%.1f", 100*h/t}')

# ── Fração mitocondrial ─────────────────────────────────────────────────────
MEAN_DEPTH="N/A"
BREADTH_PCT="N/A"
MAPPED_READS="N/A"
DIVERGENCE_PCT="N/A"
REF_LEN="N/A"

if [[ "$MITO_FRAC_OVERRIDE" != "0" ]] && [[ "$MITO_FRAC_OVERRIDE" != "0.0" ]]; then
    MITO_FRACTION="$MITO_FRAC_OVERRIDE"
    MITO_SOURCE="fornecido pelo usuário"
else
    echo "[PILOT] Mapeando amostra contra a referência com bwa mem..." >&2

    # Cópia local: os inputs do Nextflow são symlinks e os índices do bwa são
    # criados ao lado do caminho informado. O `-ef` evita o erro de copiar o
    # arquivo sobre si mesmo, caso a referência já tenha este nome.
    REF_LOCAL="pilot_ref.fa"
    [[ "$REF" -ef "$REF_LOCAL" ]] || cp -L "$REF" "$REF_LOCAL"

    bwa index "$REF_LOCAL" 2> bwa_index.log
    samtools faidx "$REF_LOCAL"

    REF_LEN=$(awk '{s+=$2} END{printf "%d", s}' "${REF_LOCAL}.fai")

    # shellcheck disable=SC2086
    bwa mem -t "$THREADS" $BWA_OPTS "$REF_LOCAL" "$R1" "$R2" 2> bwa_mem.log \
        | samtools sort -@ "$THREADS" -o pilot.bam - 2> samtools_sort.log
    samtools index pilot.bam

    # Reads mapeados primários (exclui não-mapeado, secundário, suplementar)
    MAPPED_READS=$(samtools view -c -F 0x904 pilot.bam)

    # Profundidade por posição, incluindo posições com zero (-a) e sem teto (-d 0)
    samtools depth -a -d 0 pilot.bam > depth.txt

    # Média de profundidade excluindo as bordas: numa referência curta, as
    # extremidades têm cobertura artificialmente baixa porque menos reads
    # conseguem sobrepô-las sem ultrapassar o limite. Só apara se sobrar
    # referência suficiente.
    read -r MEAN_DEPTH BREADTH_PCT <<< "$(awk -v rl="$AVG_READ_LEN" -v rlen="$REF_LEN" '
        BEGIN { trim = (rlen > 3*rl) ? rl : 0 }
        {
            pos = $2
            if (pos <= trim || pos > rlen - trim) next
            n++; sum += $3
            if ($3 > 0) cov++
        }
        END {
            if (n > 0) printf "%.4f %.1f", sum/n, 100*cov/n
            else       printf "0 0"
        }' depth.txt)"

    # Divergência média dos reads mapeados: NM / comprimento alinhado.
    # Distingue "não há mtDNA" de "a seed está distante demais" (DEC-06).
    DIVERGENCE_PCT=$(samtools view -F 0x904 pilot.bam | awk '
    {
        nm = -1
        for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) { nm = substr($i, 6) + 0; break }
        if (nm < 0) next
        cig = $6; alen = 0
        while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
            l  = substr(cig, RSTART, RLENGTH - 1) + 0
            op = substr(cig, RSTART + RLENGTH - 1, 1)
            if (op == "M" || op == "I" || op == "D" || op == "=" || op == "X") alen += l
            cig = substr(cig, RSTART + RLENGTH)
        }
        if (alen > 0) { sn += nm; sa += alen }
    }
    END { if (sa > 0) printf "%.2f", 100*sn/sa; else printf "0.00" }')

    # fração = (profundidade × genome_avg) / (read_len × reads_totais)
    MITO_FRACTION=$(awk -v d="$MEAN_DEPTH" -v g="$GENOME_AVG" \
                        -v rl="$AVG_READ_LEN" -v tr="$TOTAL_READS_BOTH" '
        BEGIN {
            if (tr <= 0 || rl <= 0) { printf "0"; exit }
            printf "%.8f", (d * g) / (rl * tr)
        }')

    if [[ "$REF_MODE" == "genome" ]]; then
        MITO_SOURCE="mapeamento bwa contra mitogenoma completo (${REF_LEN} bp)"
    else
        MITO_SOURCE="mapeamento bwa contra seed de ${REF_LEN} bp, escalado para ${GENOME_AVG} bp"
    fi
fi

# Sem sinal mitocondrial não há como recomendar volume: avisa alto e usa o teto.
NO_MITO_SIGNAL=0
if awk -v f="$MITO_FRACTION" 'BEGIN{exit !(f+0 <= 0)}'; then
    NO_MITO_SIGNAL=1
fi

# ── Cálculo da recomendação ─────────────────────────────────────────────────
# Reads mitocondriais necessários (paired-end)
MITO_READS=$(awk -v tc="$TARGET_COV" -v gs="$GENOME_AVG" -v rl="$AVG_READ_LEN" \
    'BEGIN{printf "%d", (tc * gs) / (rl * 2)}')

# Fatores de ajuste
Q30_FACTOR=$(awk -v q="$Q30_PCT" 'BEGIN{f=q/100; if(f<0.5) f=0.5; printf "%.3f", 1/f}')
# Adapter factor: reads com adaptador perdem ~30% das bases após trimming
# (maioria do adaptador aparece no fim, a maior parte das bases permanece)
ADAPTER_FACTOR=$(awk -v a="$ADAPTER_PCT" 'BEGIN{
    loss = a/100 * 0.3
    if (loss > 0.5) loss = 0.5
    printf "%.3f", 1/(1-loss)
}')

# Reads mito ajustados por qualidade
ADJUSTED=$(awk -v mr="$MITO_READS" -v qf="$Q30_FACTOR" -v af="$ADAPTER_FACTOR" \
    'BEGIN{printf "%d", mr * qf * af}')

if [[ "$NO_MITO_SIGNAL" -eq 1 ]]; then
    RAW="$MAX_READS"
else
    RAW=$(awk -v adj="$ADJUSTED" -v mf="$MITO_FRACTION" 'BEGIN{printf "%d", adj / mf}')
fi

RAW_WITH_MARGIN=$(awk -v r="$RAW" 'BEGIN{printf "%d", r * 1.5}')

# ── Janela de intake do montador ────────────────────────────────────────────
# O NOVOPlasty lê os reads sequencialmente até encher a memória e ignora o
# resto. Medido em 13/08/2026 (5 pontos, 2–8 GB, desvio de 0,005%):
#
#     reads_lidos = min( READS_PER_GB × memória ,  fornecidos × retenção )
#     cobertura   = reads_lidos × fração × read_len / genoma
#
# Consequência: baixar além da janela é desperdício puro — os reads são
# baixados, aparados e nunca lidos. Em A. leari isso eram 4× de download a mais.
# A retenção existe porque o NOVOPlasty descarta parte por subamostragem interna
# e filtro de pareamento, mesmo quando tudo caberia na memória. Ver DEC-13.
INTAKE_READS=$(awk -v pg="$READS_PER_GB" -v gb="$MEM_GB" 'BEGIN{printf "%d", pg * gb}')

# Spots a baixar para que a janela de intake fique cheia, compensando a
# retenção interna e a perda de trimming.
SPOTS_TO_FILL=$(awk -v ir="$INTAKE_READS" -v ret="$RETENTION" -v tl="$TRIM_LOSS" \
    'BEGIN{printf "%d", ir / ret / (1 - tl) / 2}')

# Cobertura máxima que esta arquitetura consegue entregar, independentemente
# do volume baixado.
COV_MAX=$(awk -v ir="$INTAKE_READS" -v mf="$MITO_FRACTION" -v rl="$AVG_READ_LEN" -v g="$GENOME_AVG" \
    'BEGIN{printf "%.1f", (ir * mf * rl) / g}')

# ── Recomendação: o menor entre o que o alvo pede e o que o montador lê ──────
# Piso removido (era da mesma família do piso de fração: existia para proteger
# de um estimador quebrado). O teto permanece como limite de segurança.
BOUND_BY="alvo"
RECOMMENDED=$RAW_WITH_MARGIN
if [[ "$SPOTS_TO_FILL" -lt "$RECOMMENDED" ]]; then
    RECOMMENDED="$SPOTS_TO_FILL"; BOUND_BY="intake"
fi
CAP_APPLIED=0
if [[ "$RECOMMENDED" -gt "$MAX_READS" ]]; then
    RECOMMENDED="$MAX_READS"; BOUND_BY="teto"; CAP_APPLIED=1
fi

RECOMMENDED=$(awk -v v="$RECOMMENDED" 'BEGIN{printf "%d", int((v + 500000) / 1000000) * 1000000}')

# Memória necessária para atingir o alvo, quando ele não é alcançável
MEM_NEEDED=$(awk -v tc="$TARGET_COV" -v cm="$COV_MAX" -v gb="$MEM_GB" \
    'BEGIN{ if (cm > 0) printf "%.1f", gb * tc / cm; else printf "?" }')

RECOMMENDED_M=$(awk -v r="$RECOMMENDED" 'BEGIN{printf "%.0f", r / 1000000}')
MITO_FRAC_PCT=$(awk -v mf="$MITO_FRACTION" 'BEGIN{printf "%.4f", mf * 100}')

# ── Bloco de alertas (DEC-10) ───────────────────────────────────────────────
ALERTS=""
if [[ "$NO_MITO_SIGNAL" -eq 1 ]]; then
    ALERTS="${ALERTS}  ⚠ NENHUM SINAL MITOCONDRIAL detectado na amostra.
    A recomendação abaixo é o TETO, não uma medição. Verifique se a seed
    corresponde ao táxon e se o dataset é WGS.
"
fi
if [[ "$CAP_APPLIED" -eq 1 ]]; then
    ALERTS="${ALERTS}  ⚠ TETO APLICADO: a recomendação bateu no limite de ${MAX_READS} reads.
    É o LIMITE, não o valor calculado.
"
fi
if [[ "$BOUND_BY" == "intake" ]]; then
    ALERTS="${ALERTS}  ℹ LIMITADO PELA JANELA DE INTAKE, não pelo alvo de cobertura.
    O cálculo pelo alvo pediria ${RAW_WITH_MARGIN} spots, mas o NOVOPlasty com
    ${MEM_GB} GB lê no máximo ${INTAKE_READS} reads — o excedente seria baixado,
    aparado e nunca lido. Baixar mais que ${RECOMMENDED} é desperdício.
"
fi

# Predição de viabilidade (DEC-13): o alvo é alcançável nesta máquina?
if awk -v cm="$COV_MAX" -v tc="$TARGET_COV" 'BEGIN{exit !(cm < tc)}'; then
    ALERTS="${ALERTS}  ⚠ ALVO INALCANÇÁVEL nesta configuração de memória.
    Cobertura máxima com ${MEM_GB} GB: ${COV_MAX}× · alvo pedido: ${TARGET_COV}×
    Nenhum volume de download levanta esse teto — quem limita é a memória.
    Para atingir ${TARGET_COV}× seriam necessários ~${MEM_NEEDED} GB.
"
fi

if [[ -z "$ALERTS" ]]; then
    ALERTS="  ✓ Recomendação derivada da medição, sem limites aplicados.
"
fi

# ── Relatório ───────────────────────────────────────────────────────────────
cat > pilot_report.txt << REPORT
╔══════════════════════════════════════════════════════════════════╗
║                       PILOT QC REPORT                           ║
╚══════════════════════════════════════════════════════════════════╝

── Amostra Piloto ─────────────────────────────────────────────────
  Reads analisados (R1):     ${TOTAL_READS}
  Reads analisados (R1+R2):  ${TOTAL_READS_BOTH}
  Comprimento médio:         ${AVG_READ_LEN} bp

── Métricas de Qualidade ──────────────────────────────────────────
  Bases ≥ Q30:               ${Q30_PCT}%
  Reads com adaptador (R1):  ${ADAPTER_PCT}% (${ADAPTER_HITS} hits)

── Fração Mitocondrial (medida) ───────────────────────────────────
  Método:                    ${MITO_SOURCE}
  Comprimento da referência: ${REF_LEN} bp
  Reads mapeados (primários):${MAPPED_READS}
  Profundidade média:        ${MEAN_DEPTH}×
  Breadth (posições ≥ 1×):   ${BREADTH_PCT}%
  Divergência média:         ${DIVERGENCE_PCT}%
  Fração medida:             ${MITO_FRACTION} (${MITO_FRAC_PCT}%)

  Métricas brutas, sem veredito automático. Breadth baixo com divergência
  alta sugere seed distante do táxon; breadth e fração baixos com
  divergência normal sugerem pouco mtDNA no dataset. Nenhum filtro de
  MAPQ é aplicado — reads de NUMTs não são discriminados nesta etapa.

── Parâmetros da Montagem ─────────────────────────────────────────
  Genoma estimado:           ${GENOME_AVG} bp
  Cobertura alvo:            ${TARGET_COV}×
  Reads mito necessários:    ${MITO_READS}
  Ajustado (qualidade):      ${ADJUSTED}
  Fator Q30:                 ${Q30_FACTOR}
  Fator adaptador:           ${ADAPTER_FACTOR}

── Janela de Intake do Montador ───────────────────────────────────
  Memória do NOVOPlasty:     ${MEM_GB} GB
  Reads que ele conseguirá ler: ${INTAKE_READS}
  Spots para saturá-la:      ${SPOTS_TO_FILL}
  ★ Cobertura MÁXIMA possível: ${COV_MAX}×   (alvo pedido: ${TARGET_COV}×)

  Baixar além de ${SPOTS_TO_FILL} spots não aumenta a cobertura: o montador
  lê sequencialmente até encher a memória e ignora o excedente.

── Recomendação ───────────────────────────────────────────────────
  Pelo alvo de cobertura:    ${RAW_WITH_MARGIN} spots (bruto ${RAW}, margem 1,5×)
  Pela janela de intake:     ${SPOTS_TO_FILL} spots
  Teto de segurança:         ${MAX_READS} spots
  Restrição que mandou:      ${BOUND_BY}
  ★ RECOMENDADO:             ${RECOMMENDED} reads (${RECOMMENDED_M}M)

${ALERTS}
  → Use: --sra_max_reads ${RECOMMENDED}
REPORT

cat pilot_report.txt >&2

echo "${RECOMMENDED}" > recommended_reads.txt
echo "${AVG_READ_LEN}" > read_length.txt
