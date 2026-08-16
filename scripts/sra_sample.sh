#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# sra_sample.sh — Extrai reads de uma entrada SRA por amostragem estratificada
#
# Substitui o padrão `fastq-dump -X N`, que pega os N PRIMEIROS spots do run.
# Esse padrão é amostragem por conveniência: FASTQ saído do bcl2fastq vem
# ordenado por tile, então os primeiros reads são de um tile de borda (viés de
# qualidade documentado) e, em runs que concatenaram bibliotecas, vêm todos da
# primeira. Ver DEC-01/DEC-02 no vault.
#
# Modos:
#   head        — N spots CONTÍGUOS começando em <start_offset_frac>×total_spots
#                 (padrão 0 = primeiros N spots, comportamento antigo). O offset é
#                 a alavanca da rearbitragem sob download contíguo (DEC-22/DEC-23):
#                 em head a SEMENTE não muda quais spots entram — só o offset muda.
#                 Por isso a retentativa desloca o início do bloco (0→50→25→75%) em
#                 vez de re-semear, que não teria efeito nenhum aqui.
#   stratified  — N spots distribuídos em W janelas equidistantes ao longo do
#                 run, com deslocamento aleatório do início de cada janela
#   dense       — igual a stratified, com W alto; aproxima amostragem aleatória
#                 do run inteiro (usado como referência no teste C2ST)
#
# A aleatoriedade vem de um LCG determinístico semeado por <seed>, não de
# rand() do awk — implementações de rand() variam entre mawk e gawk, e a
# reprodutibilidade bit-a-bit é requisito para publicação.
#
# Uso: sra_sample.sh <source> <prefix> <n_reads> <mode> <windows> <seed> [start_offset_frac]
#   source            — acesso SRA (ex: SRR28399504) ou caminho de um .sra local
#   prefix            — prefixo dos arquivos de saída (<prefix>_1.fastq, <prefix>_2.fastq)
#   start_offset_frac — só em head: fração [0,1) do run onde o bloco contíguo
#                       começa (padrão 0). Ignorado em stratified/dense (lá a
#                       semente é que gera o deslocamento).
#
# Saídas: <prefix>_1.fastq, <prefix>_2.fastq, sampling_plan.txt
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SOURCE="$1"
PREFIX="$2"
N_READS="$3"
MODE="${4:-stratified}"
WINDOWS="${5:-10}"
SEED="${6:-42}"
START_OFFSET_FRAC="${7:-0}"   # só em head; ver docstring e DEC-23

# ── Total de spots do run ───────────────────────────────────────────────────
# Só metadados, sem baixar dados.
#
# vdb-dump --info vem primeiro porque rotula o campo ("spot : 118,543,210"),
# então o parsing não depende de posição de coluna. O sra-stat --quick é o
# reserva: sua saída é posicional (campo 3 = spot count) e um formato diferente
# do esperado produziria um total errado em silêncio — que é pior que não obter
# total nenhum. Sem nenhum dos dois, caímos no modo head com aviso explícito.
# Formatos verificados contra o SRA-Toolkit 3.0.10 em 12/08/2026 (SRR36182901):
#
#   vdb-dump --info   → linha "SEQ    : 12,016,481"   (NÃO existe linha "spot :")
#   sra-stat --quick  → "SRR36182901|BARCODE|12016481:3628977262:3628977262|:|:|:"
#                        o campo 3 é spots:bases:bases_bio, separado por DOIS-PONTOS,
#                        e há uma linha por read-group (somar).
#
# Sem 'exit' nos awk de propósito: com `set -o pipefail`, sair cedo mata o
# processo anterior com SIGPIPE e o `|| total=0` descartaria um valor já lido.
get_total_spots() {
    local src="$1" from_vdb=0 from_stat=0

    from_vdb=$(vdb-dump --info "$src" 2>/dev/null \
            | awk '/^SEQ[[:space:]]*:/ && !found {gsub(/[^0-9]/,"",$0); v=$0+0; found=1}
                   END {printf "%d", v+0}') || from_vdb=0

    from_stat=$(sra-stat --meta --quick "$src" 2>/dev/null \
            | awk -F'|' 'NF>=3 { split($3, a, ":"); if (a[1] ~ /^[0-9]+$/) s += a[1] }
                         END {printf "%d", s+0}') || from_stat=0

    # Discordância entre as duas fontes indica formato inesperado: um total
    # errado mas plausível é pior que total nenhum, então avisa alto.
    if [[ "$from_vdb" -gt 0 ]] && [[ "$from_stat" -gt 0 ]]; then
        if awk -v a="$from_vdb" -v b="$from_stat" \
               'BEGIN{d=(a>b?a-b:b-a); exit !(d > 0.01*a)}'; then
            echo "[SAMPLE] AVISO: vdb-dump (${from_vdb}) e sra-stat (${from_stat}) discordam" >&2
            echo "[SAMPLE]         em mais de 1%. Usando vdb-dump; confira o formato." >&2
        fi
    fi

    if   [[ "$from_vdb"  -gt 0 ]]; then echo "$from_vdb"
    elif [[ "$from_stat" -gt 0 ]]; then echo "$from_stat"
    else echo 0
    fi
}

TOTAL_SPOTS=$(get_total_spots "$SOURCE")
echo "[SAMPLE] Total de spots no run: ${TOTAL_SPOTS:-desconhecido}" >&2

if [[ "$TOTAL_SPOTS" -eq 0 ]]; then
    echo "[SAMPLE] AVISO: não foi possível obter a contagem de spots." >&2
    echo "[SAMPLE]         Caindo para o modo 'head' (primeiros ${N_READS} spots)." >&2
    echo "[SAMPLE]         A amostra NÃO é representativa do run inteiro." >&2
    MODE="head"
elif [[ "$TOTAL_SPOTS" -le "$N_READS" ]]; then
    echo "[SAMPLE] Run menor que a amostra pedida — extraindo o run inteiro." >&2
    MODE="head"
    N_READS="$TOTAL_SPOTS"
fi

# ── Plano de amostragem ─────────────────────────────────────────────────────
if [[ "$MODE" == "head" ]]; then
    # Bloco contíguo de N spots começando em start_offset_frac×total (DEC-23).
    # offset 0 (padrão) ⇒ spots 1..N, byte-idêntico ao comportamento antigo —
    # é isto que preserva a montagem bit-a-bit da tentativa 0. A rearbitragem
    # (DEC-22) passa offset > 0 para deslocar o início: com total desconhecido
    # (=0) não há como deslocar, então fica no início.
    START=1
    if [[ "$TOTAL_SPOTS" -gt 0 ]]; then
        # O clamp mantém o bloco de N spots dentro do run: se o offset jogaria o
        # fim além do total, recua o início para que caibam N spots (último bloco
        # possível). Aritmética em awk (offset é fração), impressa como inteiro.
        START=$(awk -v t="$TOTAL_SPOTS" -v n="$N_READS" -v o="$START_OFFSET_FRAC" 'BEGIN{
            s    = int(o * t) + 1
            last = t - n + 1; if (last < 1) last = 1
            if (s > last) s = last
            if (s < 1)    s = 1
            printf "%d", s
        }')
    fi
    END=$(( START + N_READS - 1 ))
    if [[ "$TOTAL_SPOTS" -gt 0 && "$END" -gt "$TOTAL_SPOTS" ]]; then
        END="$TOTAL_SPOTS"
    fi
    printf '%d\t%d\n' "$START" "$END" > sampling_plan.txt
else
    awk -v total="$TOTAL_SPOTS" -v w="$WINDOWS" -v n="$N_READS" -v seed="$SEED" '
    BEGIN {
        per    = int(n / w); if (per < 1) per = 1
        stride = int(total / w)

        # PRNG: Park-Miller "minimal standard" — x = 16807*x mod (2^31-1).
        # Escolhido em vez do LCG do glibc (mult. 1103515245) por uma razão
        # aritmética, não estética: awk faz contas em ponto flutuante de dupla
        # precisão, exato só até 2^53. O multiplicador do glibc produziria
        # produtos de até ~2^61, que são ARREDONDADOS — e o arredondamento
        # pode variar entre mawk e gawk, quebrando justamente a
        # reprodutibilidade que a semente existe para garantir.
        # Park-Miller mantém o produto em ~2^45.6, exato em qualquer awk.
        # rand() do próprio awk foi descartado pelo mesmo motivo.
        x = int(seed) % 2147483647
        if (x <= 0) x += 2147483646

        for (i = 0; i < w; i++) {
            base = i * stride + 1
            span = stride - per
            if (span < 0) span = 0
            x = (16807 * x) % 2147483647
            jitter = (span > 0) ? x % (span + 1) : 0
            start  = base + jitter
            end    = start + per - 1
            if (start > total) continue
            if (end   > total) end = total
            printf "%d\t%d\n", start, end
        }
    }' > sampling_plan.txt
fi

N_WINDOWS=$(wc -l < sampling_plan.txt)
if [[ "$MODE" == "head" ]]; then
    echo "[SAMPLE] Modo: head — bloco contíguo, offset ${START_OFFSET_FRAC} do run" >&2
else
    echo "[SAMPLE] Modo: ${MODE} — ${N_WINDOWS} janela(s), semente ${SEED}" >&2
fi

# ── Extração ────────────────────────────────────────────────────────────────
# Cada janela vai para um diretório próprio (o fastq-dump usa sempre o mesmo
# nome de saída) e é concatenada em seguida. Os spot IDs diferem entre janelas,
# então não há colisão de nome de read no arquivo final.
: > "${PREFIX}_1.fastq"
: > "${PREFIX}_2.fastq"

WIN_IDX=0
while IFS=$'\t' read -r START END; do
    WIN_IDX=$(( WIN_IDX + 1 ))
    echo "[SAMPLE]   janela ${WIN_IDX}/${N_WINDOWS}: spots ${START}–${END}" >&2

    rm -rf "win_${WIN_IDX}"
    mkdir -p "win_${WIN_IDX}"

    # Retentativa POR JANELA: os serviços do SRA falham transitoriamente
    # ("Failed to call external services" — 2 ocorrências em 14/08/2026) e,
    # sem isto, uma janela perdida derruba a extração inteira. Backoff
    # crescente; a última tentativa deixa o erro propagar (set -e).
    FETCH_OK=0
    for TRY in 1 2 3; do
        # --split-3 (não --split-files): um read só entra em _1/_2 se AMBOS os
        # mates existirem; singletons vão para <acc>.fastq (sem sufixo _N),
        # ignorado pelo glob abaixo. Isso mantém _1 e _2 SEMPRE sincronizados.
        # Com --split-files, spots com um mate só (comuns em dado de captura,
        # ex.: SRR14323300) desalinham _1/_2 e o bwa aborta lá na frente
        # ("paired reads have different names"). --skip-technical descarta
        # leituras técnicas (adaptadores/barcodes).
        if fastq-dump \
            -N "$START" \
            -X "$END" \
            --split-3 \
            --skip-technical \
            --outdir "win_${WIN_IDX}" \
            "$SOURCE"; then
            FETCH_OK=1
            break
        fi
        echo "[SAMPLE]   janela ${WIN_IDX}: tentativa ${TRY}/3 falhou; aguardando $(( TRY * 20 ))s" >&2
        sleep $(( TRY * 20 ))
    done
    if [[ "$FETCH_OK" -ne 1 ]]; then
        echo "[SAMPLE] ERRO: janela ${WIN_IDX} falhou após 3 tentativas." >&2
        exit 1
    fi

    # --split-files nomeia pelo acesso, não pelo prefixo pedido
    for MATE in 1 2; do
        # -print -quit em vez de `| head -1`: evita SIGPIPE no find sob pipefail
        FOUND=$(find "win_${WIN_IDX}" -maxdepth 1 -name "*_${MATE}.fastq" -print -quit)
        if [[ -z "$FOUND" ]]; then
            echo "[SAMPLE] ERRO: janela ${WIN_IDX} não produziu o mate ${MATE}." >&2
            exit 1
        fi
        cat "$FOUND" >> "${PREFIX}_${MATE}.fastq"
    done

    rm -rf "win_${WIN_IDX}"
done < sampling_plan.txt

OBTAINED=$(awk 'END{printf "%d", NR/4}' "${PREFIX}_1.fastq")
OBTAINED_R2=$(awk 'END{printf "%d", NR/4}' "${PREFIX}_2.fastq")
echo "[SAMPLE] Reads obtidos: R1=${OBTAINED} R2=${OBTAINED_R2}" >&2

# Guarda de sincronia: R1 e R2 têm de ter o mesmo número de reads. Com --split-3
# isto deve valer sempre; a checagem existe para um descompasso futuro falhar
# aqui, com mensagem clara, em vez de virar um erro obscuro do bwa/montador.
if [[ "$OBTAINED" -ne "$OBTAINED_R2" ]]; then
    echo "[SAMPLE] ERRO: R1 (${OBTAINED}) e R2 (${OBTAINED_R2}) dessincronizados." >&2
    exit 1
fi

if [[ "$OBTAINED" -eq 0 ]]; then
    echo "[SAMPLE] ERRO: nenhuma read extraída." >&2
    exit 1
fi
