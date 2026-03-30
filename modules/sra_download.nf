/*
 * Módulo: SRA_DOWNLOAD
 * Baixa reads pareados de uma entrada do NCBI SRA usando fasterq-dump.
 *
 * Entradas : acesso SRA (ex: 'SRR36152783')
 * Saídas   : tuple (sample_id, read1.fastq, read2.fastq)
 */

process SRA_DOWNLOAD {

    tag "${accession}"

    publishDir "${params.outdir}/reads/raw", mode: 'copy', pattern: '*.fastq'

    input:
    val accession

    output:
    tuple val(accession),
          path("${accession}_1.fastq"),
          path("${accession}_2.fastq"),
          emit: reads

    script:
    // -X N = baixar apenas os primeiros N spots (para testes com datasets limitados)
    def max_reads = params.sra_max_reads ? "-X ${params.sra_max_reads}" : ""

    """
    # Configuração não-interativa do SRA-Toolkit
    vdb-config --restore-defaults 2>/dev/null || true

    # Baixa e converte para FASTQ pareado
    fasterq-dump \\
        ${accession} \\
        --split-files \\
        --threads ${task.cpus} \\
        --outdir . \\
        --temp /tmp \\
        ${max_reads}
    """
}
