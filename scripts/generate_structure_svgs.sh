#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────────
# generate_structure_svgs.sh
# Gera SVGs de estrutura secundária de tRNAs e rRNAs a partir da saída do MITOS2
# Utiliza RNAplot (ViennaRNA) disponível na imagem mitos2:1.0
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

MITOS_DIR="${1:?Uso: $0 <diretório_saída_mitos2>}"
OUTDIR="${2:-${MITOS_DIR}/structure_svgs}"

MITFI_DIR="${MITOS_DIR}/mitfi-global"

if [[ ! -d "$MITFI_DIR" ]]; then
    echo "ERRO: diretório mitfi-global não encontrado em ${MITOS_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}/tRNA" "${OUTDIR}/rRNA"

# ──────────────────────────────────────────────────────────────────────────────
# Função: extrai sequências + estruturas do .nc e gera SVG via RNAplot
# Filtra apenas hits primários (|1|1| no header = best hit)
# ──────────────────────────────────────────────────────────────────────────────
generate_svgs() {
    local nc_file="$1"
    local out_subdir="$2"
    local type_label="$3"

    if [[ ! -f "$nc_file" ]]; then
        echo "AVISO: ${nc_file} não encontrado, pulando ${type_label}..." >&2
        return
    fi

    local count=0
    local name="" seq="" struct=""

    while IFS= read -r line; do
        if [[ "$line" == ">"* ]]; then
            # Processa entrada anterior se existir
            if [[ -n "$name" && -n "$seq" && -n "$struct" ]]; then
                local svgfile="${out_subdir}/${name}.svg"
                echo -e ">${name}\n${seq}\n${struct}" | \
                    RNAplot --output-format=svg --auto-id 2>/dev/null
                # RNAplot cria arquivo no diretório atual, mover para destino
                local rnaplot_out
                rnaplot_out=$(ls -t *_ss.svg 2>/dev/null | head -1)
                if [[ -n "$rnaplot_out" ]]; then
                    mv "$rnaplot_out" "$svgfile"
                    count=$((count + 1))
                fi
            fi

            # Parse header: >Contig1|code|start|end|strand|1|len|score|...|anticodon|name|model|1|1|mode
            # Filtra apenas hits primários (penúltimo e antepenúltimo campo = 1|1)
            local is_primary
            is_primary=$(echo "$line" | awk -F'|' '{print $(NF-2)"|"$(NF-1)}')
            if [[ "$is_primary" != "1|1" ]]; then
                name=""
                continue
            fi

            # Extrai nome do gene (ex: trnF, rrnL)
            name=$(echo "$line" | awk -F'|' '{print $(NF-3)}')
            seq=""
            struct=""
        elif [[ -n "$name" ]]; then
            if [[ -z "$seq" ]]; then
                seq="$line"
            else
                struct="$line"
            fi
        fi
    done < "$nc_file"

    # Processa última entrada
    if [[ -n "$name" && -n "$seq" && -n "$struct" ]]; then
        local svgfile="${out_subdir}/${name}.svg"
        echo -e ">${name}\n${seq}\n${struct}" | \
            RNAplot --output-format=svg --auto-id 2>/dev/null
        local rnaplot_out
        rnaplot_out=$(ls -t *_ss.svg 2>/dev/null | head -1)
        if [[ -n "$rnaplot_out" ]]; then
            mv "$rnaplot_out" "$svgfile"
            count=$((count + 1))
        fi
    fi

    echo "${type_label}: ${count} SVGs gerados em ${out_subdir}"
}

# Diretório temporário para RNAplot (ele cria arquivos no cwd)
TMPDIR=$(mktemp -d)
cd "$TMPDIR"

generate_svgs "${MITFI_DIR}/sequence.fas-0_tRNAout.nc" "${OUTDIR}/tRNA" "tRNA"
generate_svgs "${MITFI_DIR}/sequence.fas-0_rRNAout.nc" "${OUTDIR}/rRNA" "rRNA"

# Limpa temp
rm -rf "$TMPDIR"

echo ""
echo "Estruturas secundárias geradas em: ${OUTDIR}"
ls -la "${OUTDIR}/tRNA/" "${OUTDIR}/rRNA/"
