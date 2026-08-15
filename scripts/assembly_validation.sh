#!/bin/bash
# assembly_validation.sh — validação pós-montagem por massa de profundidade (DEC-20)
#
# Mapeia a AMOSTRA PILOTO de volta contra a montagem circularizada e procura
# regiões cuja profundidade destoa da mediana do genoma:
#
#   razão << 1  → repeat SUPEREXPANDIDO na montagem (cópias a mais: os reads
#                 reais se espalham por cópias que não existem no genoma)
#   razão >> 1  → repeat COLAPSADO (cópias a menos: reads de várias cópias
#                 verdadeiras se empilham numa só)
#
# Motivação: em 14/08/2026 o NOVOPlasty circularizou A. leari com 18.058 bp —
# um VNTR de ~298 bp da região controle superexpandido (~4,8 cópias contra
# ~1–2 sustentadas pelos reads; razão de profundidade 0,22). O comprimento
# passava no genome_range; só a massa de profundidade denunciou. Ver DEC-20.
#
# A checagem de comprimento contra a referência (se fornecida) é informativa;
# o detector primário é a razão de profundidade, que LOCALIZA o problema.
#
# Uso: assembly_validation.sh <asm.fasta> <R1> <R2> [win] [low] [high] [ref|NONE] [threads]

set -euo pipefail
export LC_ALL=C

ASM=$1
R1=$2
R2=$3
WIN=${4:-200}
LOW=${5:-0.5}
HIGH=${6:-2.0}
REF=${7:-NONE}
THREADS=${8:-4}

# Bordas do círculo linearizado têm profundidade deprimida por construção
# (reads que atravessam a origem não mapeiam inteiros) — fora da varredura.
MARGIN=300

ASM_LEN=$(awk '!/^>/{l+=length($0)} END{print l}' "$ASM")

bwa index "$ASM" 2>/dev/null
bwa mem -t "$THREADS" "$ASM" "$R1" "$R2" 2>/dev/null \
    | samtools view -b -F 0x904 - 2>/dev/null \
    | samtools sort -@ 2 -o asm_val.bam - 2>/dev/null
samtools index asm_val.bam
samtools depth -a asm_val.bam > depth.tsv

MAPPED=$(samtools view -c asm_val.bam)
NPOS=$(wc -l < depth.tsv)

if [[ "$NPOS" -eq 0 || "$MAPPED" -lt 100 ]]; then
    {
        echo "VEREDITO: INDETERMINADO — só $MAPPED reads do piloto mapearam na montagem."
        echo "Sem massa para validar. Possíveis causas: piloto muito raso ou montagem espúria."
    } > validation_report.txt
    : > depth_windows.tsv
    echo "INDETERMINADO" > verdict.txt
    exit 0
fi

# Mediana (robusta ao próprio repeat) — sem 'exit' no awk consumidor: SIGPIPE
# sob pipefail mata o sort (armadilha documentada em 12/08).
MED=$(cut -f3 depth.tsv | sort -n | awk -v n="$NPOS" 'NR==int((n+1)/2){m=$0} END{print m}')
MEAN=$(awk -F'\t' '{s+=$3} END{printf "%.2f", s/NR}' depth.tsv)
BREADTH=$(awk -F'\t' '$3>=1{c++} END{printf "%.1f", 100*c/NR}' depth.tsv)

if [[ "$MED" -lt 2 ]]; then
    {
        echo "VEREDITO: INDETERMINADO — profundidade mediana do piloto na montagem é ${MED}×."
        echo "Abaixo de 2× a razão por janela é ruído. Aumente pilot_reads se quiser validar."
    } > validation_report.txt
    : > depth_windows.tsv
    echo "INDETERMINADO" > verdict.txt
    exit 0
fi

# Perfil por janela (arquivo completo, útil para figura)
awk -F'\t' -v win="$WIN" '
    { w = int(($2-1)/win); s[w]+=$3; c[w]++ }
    END {
        maxw = int((NR-1)/win) + 1
        for (w = 0; w <= maxw; w++)
            if (c[w] > 0) printf "%d\t%d\t%.2f\n", w*win+1, w*win+c[w], s[w]/c[w]
    }' depth.tsv > depth_windows.tsv

# Varredura: janelas fora de [LOW, HIGH]×mediana, com margens excluídas,
# fundidas em regiões contíguas de mesmo tipo
awk -F'\t' -v win="$WIN" -v med="$MED" -v lo="$LOW" -v hi="$HIGH" \
    -v margin="$MARGIN" -v len="$ASM_LEN" '
    $2 > margin && $2 <= len - margin {
        w = int(($2 - margin - 1) / win)
        s[w] += $3; c[w]++
    }
    END {
        maxw = int((len - 2*margin - 1) / win)
        inreg = 0
        for (w = 0; w <= maxw; w++) {
            if (c[w] == 0) continue
            m = s[w] / c[w]
            r = m / med
            flag = (r < lo) ? "SUPEREXPANSAO" : (r > hi ? "COLAPSO" : "")
            start = margin + w*win + 1
            end   = margin + w*win + c[w]
            if (flag != "") {
                if (inreg && flag == curflag) {
                    curend = end; sum += m; k++
                } else {
                    if (inreg) printf "%d\t%d\t%.2f\t%.2f\t%s\n", curstart, curend, sum/k, (sum/k)/med, curflag
                    curstart = start; curend = end; sum = m; k = 1; curflag = flag
                    inreg = 1
                }
            } else if (inreg) {
                printf "%d\t%d\t%.2f\t%.2f\t%s\n", curstart, curend, sum/k, (sum/k)/med, curflag
                inreg = 0
            }
        }
        if (inreg) printf "%d\t%d\t%.2f\t%.2f\t%s\n", curstart, curend, sum/k, (sum/k)/med, curflag
    }' depth.tsv > all_regions.tsv

# Uma janela ISOLADA fora da faixa é compatível com ruído amostral (a ~11× de
# piloto, ~15 reads/janela → CV ~26%) ou com viés de cobertura do sequenciamento
# — os mesmos dips aparecem em montagens sabidamente corretas. Anomalia de
# repeat real se sustenta por ≥ 2 janelas consecutivas (o VNTR do caso-motivo
# ocupou 7). Isoladas viram observação informativa, não alerta.
MINSPAN=$((WIN * 2))
awk -F'\t' -v ms="$MINSPAN" '($2-$1+1) >= ms' all_regions.tsv > flagged_regions.tsv
awk -F'\t' -v ms="$MINSPAN" '($2-$1+1) <  ms' all_regions.tsv > isolated_windows.tsv

NFLAG=$(wc -l < flagged_regions.tsv)
NISO=$(wc -l < isolated_windows.tsv)

# Comprimento vs referência (informativo)
LENLINE="(sem referência para comparar)"
if [[ "$REF" != "NONE" && -f "$REF" ]]; then
    REF_LEN=$(awk '!/^>/{l+=length($0)} END{print l}' "$REF")
    DELTA=$(awk -v a="$ASM_LEN" -v r="$REF_LEN" 'BEGIN{printf "%+.2f", 100.0*(a-r)/r}')
    LENLINE="montagem ${ASM_LEN} bp vs referência ${REF_LEN} bp (Δ ${DELTA}%)"
fi

VEREDITO="OK — nenhuma região com razão de profundidade fora de [${LOW}, ${HIGH}]"
VERDICT_WORD="OK"
if [[ "$NFLAG" -gt 0 ]]; then
    VEREDITO="ALERTA — ${NFLAG} região(ões) com razão de profundidade anômala"
    VERDICT_WORD="ALERTA"
fi
# Veredito legível por máquina — é o gatilho da rearbitragem (DEC-22)
echo "$VERDICT_WORD" > verdict.txt

{
    echo "╔══════════════════════════════════════════════════════════════════╗"
    echo "║              VALIDAÇÃO PÓS-MONTAGEM (massa de profundidade)      ║"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Montagem:              $(basename "$ASM") (${ASM_LEN} bp)"
    echo "Piloto mapeado:        ${MAPPED} reads"
    echo "Profundidade:          mediana ${MED}× · média ${MEAN}× · breadth ≥1×: ${BREADTH}%"
    echo "Comprimento:           ${LENLINE}"
    echo "Varredura:             janelas de ${WIN} bp, margens de ${MARGIN} bp excluídas"
    echo ""
    echo "VEREDITO: ${VEREDITO}"
    if [[ "$NISO" -gt 0 ]]; then
        echo ""
        echo "Janelas isoladas fora da faixa (informativo — compatível com ruído"
        echo "amostral ou viés de cobertura; NÃO geram alerta):"
        awk -F'\t' '{printf "  %6d–%-6d  %6.2f×  razão %.2f  (%s)\n", $1, $2, $3, $4, $5}' isolated_windows.tsv
    fi
    if [[ "$NFLAG" -gt 0 ]]; then
        echo ""
        echo "Regiões sinalizadas (início, fim, profundidade, razão, tipo):"
        awk -F'\t' '{printf "  %6d–%-6d  %6.2f×  razão %.2f  %s\n", $1, $2, $3, $4, $5}' flagged_regions.tsv
        echo ""
        echo "Interpretação:"
        echo "  SUPEREXPANSAO (razão < ${LOW}): a montagem tem cópias A MAIS de um repeat"
        echo "    — os reads reais se espalham por cópias inexistentes no genoma."
        echo "  COLAPSO (razão > ${HIGH}): a montagem tem cópias A MENOS"
        echo "    — reads de várias cópias verdadeiras empilhados numa só."
        echo "  A montagem circularizou, mas circular ≠ correto. Compare com execução"
        echo "  usando outro pool de reads (outra semente/volume) antes de usar."
    fi
} > validation_report.txt

rm -f depth.tsv asm_val.bam asm_val.bam.bai
