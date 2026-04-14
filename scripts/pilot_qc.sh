#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# Pilot QC — Analisa uma amostra piloto de reads e recomenda max_reads ideal
#
# Métricas extraídas:
#   1. Porcentagem de bases ≥ Q30 (qualidade Phred)
#   2. Porcentagem de reads com adaptador Illumina
#   3. Comprimento médio dos reads
#   4. Fração mitocondrial estimada (via seed kmer grep)
#
# Cálculo:
#   mito_reads = (target_cov × genome_size) / (read_len × 2)
#   adjusted   = mito_reads × quality_factor × adapter_factor
#   total      = adjusted / mito_fraction × safety_margin (1.5×)
#
# Uso: pilot_qc.sh <R1> <R2> <seed> <genome_range> <target_cov> [mito_frac]
# Saídas: pilot_report.txt, recommended_reads.txt
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

R1="$1"
R2="$2"
SEED="$3"
GENOME_RANGE="$4"
TARGET_COV="${5:-500}"
MITO_FRAC_OVERRIDE="${6:-0}"

# ── Parse genome range ──────────────────────────────────────────────────────
GENOME_MIN=$(echo "$GENOME_RANGE" | cut -d'-' -f1)
GENOME_MAX=$(echo "$GENOME_RANGE" | cut -d'-' -f2)
GENOME_AVG=$(( (GENOME_MIN + GENOME_MAX) / 2 ))

echo "[PILOT] Analisando amostra piloto..." >&2

# ── Métricas básicas ────────────────────────────────────────────────────────
TOTAL_READS=$(awk 'END{printf "%d", NR/4}' "$R1")
AVG_READ_LEN=$(awk 'NR%4==2 {total+=length($0); n++} END{printf "%d", total/n}' "$R1")

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
if [[ "$MITO_FRAC_OVERRIDE" != "0" ]] && [[ "$MITO_FRAC_OVERRIDE" != "0.0" ]]; then
    MITO_FRACTION="$MITO_FRAC_OVERRIDE"
    MITO_SOURCE="fornecido pelo usuário"
    SEED_HITS="N/A"
else
    # Extrai k-mers de 20 bp a cada 50 bp ao longo da semente inteira
    # Mais k-mers = melhor sensibilidade mesmo com divergência inter-específica
    SEED_SEQ=$(grep -v '^>' "$SEED" | tr -d '\n\r\t ')
    SEED_LEN=${#SEED_SEQ}
    KMER_LEN=20
    TOTAL_HITS=0
    N_KMERS=0

    if [[ "$SEED_LEN" -ge "$KMER_LEN" ]]; then
        STEP=50
        POS=0
        while [[ "$POS" -le $(( SEED_LEN - KMER_LEN )) ]]; do
            KMER="${SEED_SEQ:$POS:$KMER_LEN}"
            HITS=$(grep -c "$KMER" "$R1" 2>/dev/null) || HITS=0
            TOTAL_HITS=$(( TOTAL_HITS + HITS ))
            N_KMERS=$(( N_KMERS + 1 ))
            POS=$(( POS + STEP ))
        done
    fi

    SEED_HITS="$TOTAL_HITS"

    if [[ "$TOTAL_HITS" -gt 0 ]]; then
        MITO_FRACTION=$(awk -v hits="$TOTAL_HITS" -v nk="$N_KMERS" -v tr="$TOTAL_READS" \
            'BEGIN{
                est = (hits / nk) * 2
                f = est / (tr * 2)
                if (f < 0.0005) f = 0.0005
                if (f > 0.05)   f = 0.05
                printf "%.6f", f
            }')
        MITO_SOURCE="estimativa por seed kmers"
    else
        MITO_FRACTION="0.001"
        MITO_SOURCE="padrão conservador (0.1%)"
    fi
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

# Total de reads necessários
RAW=$(awk -v adj="$ADJUSTED" -v mf="$MITO_FRACTION" 'BEGIN{printf "%d", adj / mf}')

# Margem de segurança (1.5×) + mínimo 5M + máximo 25M + arredonda p/ milhão
# Cap de 25M: 20M reads → 306× (A. leari), NOVOPlasty usa ~23% → 25M garante >150×
RECOMMENDED=$(awk -v r="$RAW" 'BEGIN{
    v = r * 1.5
    if (v < 5000000)  v = 5000000
    if (v > 25000000) v = 25000000
    v = int((v + 500000) / 1000000) * 1000000
    printf "%d", v
}')

RECOMMENDED_M=$(awk -v r="$RECOMMENDED" 'BEGIN{printf "%.0f", r / 1000000}')

# ── Relatório ───────────────────────────────────────────────────────────────
MITO_FRAC_PCT=$(awk -v mf="$MITO_FRACTION" 'BEGIN{printf "%.3f", mf * 100}')
RAW_WITH_MARGIN=$(awk -v r="$RAW" 'BEGIN{printf "%d", r * 1.5}')

cat > pilot_report.txt << REPORT
╔══════════════════════════════════════════════════════════════════╗
║                       PILOT QC REPORT                           ║
╚══════════════════════════════════════════════════════════════════╝

── Amostra Piloto ─────────────────────────────────────────────────
  Reads analisados (R1):     ${TOTAL_READS}
  Comprimento médio:         ${AVG_READ_LEN} bp

── Métricas de Qualidade ──────────────────────────────────────────
  Bases ≥ Q30:               ${Q30_PCT}%
  Reads com adaptador (R1):  ${ADAPTER_PCT}% (${ADAPTER_HITS} hits)

── Fração Mitocondrial ────────────────────────────────────────────
  Método:                    ${MITO_SOURCE}
  Seed kmer hits:            ${SEED_HITS}
  Fração estimada:           ${MITO_FRACTION} (${MITO_FRAC_PCT}%)

── Parâmetros da Montagem ─────────────────────────────────────────
  Genoma estimado:           ${GENOME_AVG} bp
  Cobertura alvo:            ${TARGET_COV}×
  Reads mito necessários:    ${MITO_READS}
  Ajustado (qualidade):      ${ADJUSTED}
  Fator Q30:                 ${Q30_FACTOR}
  Fator adaptador:           ${ADAPTER_FACTOR}

── Recomendação ───────────────────────────────────────────────────
  Cálculo bruto:             ${RAW} reads
  Com margem segurança:      ${RAW_WITH_MARGIN} reads
  ★ RECOMENDADO:             ${RECOMMENDED} reads (${RECOMMENDED_M}M)

  → Use: --sra_max_reads ${RECOMMENDED}
REPORT

cat pilot_report.txt >&2

echo "${RECOMMENDED}" > recommended_reads.txt
echo "${AVG_READ_LEN}" > read_length.txt
