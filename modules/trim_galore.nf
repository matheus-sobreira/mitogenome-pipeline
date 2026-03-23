/*
 * Módulo: TRIM_GALORE
 * Remove adaptadores e bases de baixa qualidade das leituras paired-end.
 *
 * Entradas : tuple (sample_id, read1.fastq, read2.fastq)
 * Saídas   : tuple (sample_id, read1_val_1.fq, read2_val_2.fq)
 *            + relatórios de trimming
 */

process TRIM_GALORE {

    tag "${sample_id}"

    publishDir "${params.outdir}/reads/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("*_val_1.fq"),
          path("*_val_2.fq"),
          emit: reads
    path "*_trimming_report.txt", emit: reports

    script:
    """
    trim_galore \\
        --paired \\
        --quality   ${params.trim_quality} \\
        --length    ${params.min_length}   \\
        --cores     ${task.cpus}           \\
        ${read1}                           \\
        ${read2}
    """
}
