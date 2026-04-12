/*
 * Módulo: SRA_DOWNLOAD
 * Baixa reads pareados de uma entrada do NCBI SRA usando fasterq-dump.
 *
 * Entradas : acesso SRA (ex: 'SRR36152783'), max_reads (0 = baixar tudo)
 * Saídas   : tuple (sample_id, read1.fastq, read2.fastq)
 */

process SRA_DOWNLOAD {

    tag "${accession}"

    // Downloads do NCBI podem falhar por instabilidade de rede
    errorStrategy 'retry'
    maxRetries    2

    // Sem publishDir — reads brutos serão apagados após uso para economizar disco

    input:
    val accession
    val max_reads

    output:
    tuple val(accession),
          path("${accession}_1.fastq"),
          path("${accession}_2.fastq"),
          emit: reads

    script:
    """
    # Habilita SDL2 — necessário no SRA-Toolkit 3.x
    mkdir -p \$HOME/.ncbi
    cat > \$HOME/.ncbi/user-settings.mkfg << 'MKFG'
/LIBS/GUID = "d3b07384-d9a3-4f3c-b7e2-4c5c1d3abe42"
/repository/remote/main/SDL.2/resolver-cgis/enabled = "true"
/repository/user/cache-disabled = "true"
MKFG

    # Etapa 1: prefetch — baixa o .sra completo localmente (mais confiável que fasterq-dump remoto)
    prefetch \\
        ${accession} \\
        --max-size 50G \\
        --output-directory .

    # Etapa 2: fasterq-dump — converte o .sra local para FASTQ
    fasterq-dump \\
        ${accession}/${accession}.sra \\
        --split-files \\
        --threads ${task.cpus} \\
        --outdir . \\
        --temp /tmp

    # Libera espaço: apaga o diretório .sra (não é mais necessário)
    rm -rf ${accession}

    # Etapa 3: se sra_max_reads definido, trunca os FASTQs para economizar disco e tempo
    if [ ${max_reads} -gt 0 ]; then
        MAX_LINES=\$(( ${max_reads} * 4 ))
        echo "Truncando FASTQs para ${max_reads} reads (\$MAX_LINES linhas)..."
        head -n \$MAX_LINES ${accession}_1.fastq > tmp_1.fastq && mv tmp_1.fastq ${accession}_1.fastq
        head -n \$MAX_LINES ${accession}_2.fastq > tmp_2.fastq && mv tmp_2.fastq ${accession}_2.fastq
    fi
    """
}
