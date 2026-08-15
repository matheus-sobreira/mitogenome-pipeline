/*
 * Módulo: ASSEMBLY_VALIDATION
 * Valida a montagem circularizada mapeando a AMOSTRA PILOTO de volta contra
 * ela e medindo a razão de profundidade por janela (DEC-20).
 *
 * Roda em paralelo ao MITOS2 — não bloqueia a anotação; o veredito vai para
 * o relatório publicado e para o resumo do workflow.onComplete.
 *
 * Entradas : tuple (sample_id, assembly.fasta, pilot_R1, pilot_R2)
 *            path(script de validação)
 *            path(referência p/ comparação de comprimento — a seed serve de
 *                 placeholder quando params.pilot_reference não existe)
 *
 * Saídas   : validation_report.txt  → veredito + regiões sinalizadas
 *            depth_windows.tsv      → perfil de profundidade (para figuras)
 */

process ASSEMBLY_VALIDATION {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc/assembly_validation", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(read1), path(read2)
    path  validation_script
    path  length_ref

    output:
    path 'validation_report.txt', emit: report
    path 'depth_windows.tsv',     emit: profile
    tuple val(sample_id), path('verdict.txt'), emit: verdict

    script:
    def refArg = params.pilot_reference ? "${length_ref}" : 'NONE'
    """
    #!/bin/bash
    set -euo pipefail

    bash ${validation_script} \\
        ${assembly} ${read1} ${read2} \\
        ${params.validation_window} \\
        ${params.validation_ratio_low} \\
        ${params.validation_ratio_high} \\
        ${refArg} \\
        ${task.cpus}

    echo "=== Validação pós-montagem (${sample_id}) ==="
    cat validation_report.txt
    """
}
