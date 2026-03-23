/*
 * Módulo: SRA_DOWNLOAD
 * Baixa reads pareados de uma entrada do NCBI SRA usando fasterq-dump.
 *
 * Entradas : acesso SRA (ex: 'SRR2081280')
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
    // Limite de reads: usado apenas nos testes (sra_max_reads = 500000)
    def max_reads = params.sra_max_reads ? "--maxSpotId ${params.sra_max_reads}" : ""

    """
    # Tenta restaurar configuração padrão do SRA-Toolkit
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
