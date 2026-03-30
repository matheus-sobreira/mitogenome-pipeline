/*
 * Módulo: FASTQC
 * Avalia a qualidade das leituras brutas (pré-trimming).
 *
 * Entradas : tuple (sample_id, read1.fastq, read2.fastq)
 * Saídas   : relatórios HTML e ZIP
 */

process FASTQC {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("*.html"),
          path("*.zip"),
          emit: reports

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir .             \\
        ${read1}               \\
        ${read2}
    """
}
