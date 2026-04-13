#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# run_pipeline.sh — Executa o pipeline de forma resiliente a desconexão
#
# Uso:
#   ./run_pipeline.sh -profile a_leari
#   ./run_pipeline.sh -profile test
#   ./run_pipeline.sh -profile a_leari -resume -w work/a_leari/2026-04-12_14h30
#   ./run_pipeline.sh -profile a_leari --sra_max_reads 20000000
#
# O pipeline roda em background com nohup. O terminal pode ser fechado
# sem interromper a execução. O log é salvo em logs/<especie>_<timestamp>.log
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Extrair nome do perfil para o log ────────────────────────────────────────
PROFILE="unknown"
prev=""
for i in "$@"; do
    if [[ "$prev" == "-profile" ]]; then
        PROFILE="$i"
        break
    fi
    prev="$i"
done

TIMESTAMP=$(date +"%Y-%m-%d_%Hh%M")
LOGDIR="logs"
mkdir -p "$LOGDIR"
LOGFILE="${LOGDIR}/${PROFILE}_${TIMESTAMP}.log"

# ── Verificar pré-requisitos ─────────────────────────────────────────────────
if ! command -v nextflow &>/dev/null; then
    echo "ERRO: nextflow não encontrado no PATH"
    exit 1
fi

if ! docker info &>/dev/null; then
    echo "ERRO: Docker não está rodando"
    exit 1
fi

# ── Executar pipeline em background ──────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════"
echo "  mitogenome-pipeline"
echo "  Perfil:    $PROFILE"
echo "  Início:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Log:       $LOGFILE"
echo "════════════════════════════════════════════════════════════════════"
echo ""
echo "  O pipeline está rodando em background."
echo "  Você pode fechar o terminal com segurança."
echo ""
echo "  Comandos úteis:"
echo "    tail -f $LOGFILE              # acompanhar em tempo real"
echo "    grep -E 'process|ERROR' $LOGFILE   # ver progresso/erros"
echo "    cat $LOGFILE | tail -20       # ver final do log"
echo ""
echo "════════════════════════════════════════════════════════════════════"

nohup nextflow run main.nf "$@" > "$LOGFILE" 2>&1 &
PIPELINE_PID=$!

echo "  PID:       $PIPELINE_PID"
echo "  Para cancelar:  kill $PIPELINE_PID"
echo "════════════════════════════════════════════════════════════════════"

# Salvar PID para referência
echo "$PIPELINE_PID" > "${LOGDIR}/.last_pid"
