/*
 * Módulo: COMPILE_SUMMARY
 * Compila todos os resultados do pipeline em uma pasta organizada
 * com os entregáveis para apresentação e submissão ao GenBank.
 *
 * Entradas : tuple (sample_id, assembly)  ← montagem circularizada
 *            tuple (sample_id, mitos_dir) ← diretório MITOS2
 *
 * Saídas   : pasta summary/ com 14 entregáveis organizados
 */

process COMPILE_SUMMARY {

    tag "${sample_id}"

    publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly), path(mitos_dir)

    output:
    path "${sample_id}/*", emit: deliverables

    script:
    def organism = params.organism ?: "Unknown organism"
    def genbank_arg = ""
    """
    # Gerar arquivos GenBank (.tbl + .fsa)
    python3 ${projectDir}/scripts/gff2genbank.py \\
        --gff   ${mitos_dir}/result.gff \\
        --fasta ${assembly} \\
        --organism "${organism}" \\
        --transl-table ${params.genetic_code} \\
        --outdir genbank_tmp/

    # Compilar pasta de entregáveis
    python3 ${projectDir}/scripts/compile_summary.py \\
        --assembly  ${assembly} \\
        --mitos-dir ${mitos_dir} \\
        --genbank-dir genbank_tmp/ \\
        --organism "${organism}" \\
        --outdir ${sample_id}
    """
}
