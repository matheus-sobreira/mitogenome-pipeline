/*
 * Módulo: SRA_DOWNLOAD
 * Baixa reads pareados de uma entrada do NCBI SRA usando o SRA-Toolkit.
 *
 * Quando max_reads > 0, a extração é estratificada ao longo do run em vez de
 * pegar os primeiros N spots (ver DEC-02 no vault). O custo é zero: o prefetch
 * da etapa 1 já trouxe o .sra completo para o disco local, então extrair
 * janelas espalhadas não adiciona nenhuma requisição de rede. Sem isso, a
 * montagem continuaria alimentada pelo mesmo viés de localidade que o Pilot QC
 * passou a corrigir.
 *
 * Entradas : acesso SRA (ex: 'SRR36152783'), max_reads (0 = baixar tudo)
 * Saídas   : tuple (sample_id, read1.fastq, read2.fastq)
 */

process SRA_DOWNLOAD {

    tag "${accession}"

    // Downloads do NCBI podem falhar por instabilidade de rede
    errorStrategy 'retry'
    maxRetries    2

    // Só o plano de amostragem é publicado — reads brutos são apagados após uso
    publishDir "${params.outdir}/qc/download", mode: 'copy', overwrite: true,
               pattern: 'sampling_plan.txt'

    input:
    val accession
    val max_reads
    path sampler

    output:
    tuple val(accession),
          path("${accession}_1.fastq"),
          path("${accession}_2.fastq"),
          emit: reads
    path "sampling_plan.txt", optional: true, emit: plan
    // .sra preservado quando a rearbitragem está ativa (DEC-22): o re-sorteio
    // do pool extrai janelas novas dele SEM tocar a rede. Symlink downstream,
    // sem cópia — custo é só o disco do work dir até a limpeza.
    path "${accession}.sra", optional: true, emit: sra

    script:
    def mode     = params.sra_sampling ?: 'stratified'
    def windows  = params.sra_windows  ?: 50
    def seed_rng = params.pilot_seed   ?: 42

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

    # Etapa 2: converte o .sra local para FASTQ
    if [ ${max_reads} -gt 0 ]; then
        # Amostragem estratificada sobre o .sra local — evita gravar o dataset
        # inteiro em disco (crucial para metagenomas) sem herdar o viés de
        # pegar só o início do run.
        echo "Extraindo ${max_reads} spots (modo '${mode}', ${windows} janelas)..."
        bash ${sampler} \\
            ${accession}/${accession}.sra \\
            ${accession} \\
            ${max_reads} \\
            ${mode} \\
            ${windows} \\
            ${seed_rng}
    else
        # Sem limite: usa fasterq-dump (mais rápido, multi-thread)
        fasterq-dump \\
            ${accession}/${accession}.sra \\
            --split-files \\
            --threads ${task.cpus} \\
            --outdir . \\
            --temp . \\
            --disk-limit unlimited
    fi

    # Rearbitragem ativa: preserva o .sra para re-sorteio de pool sem rede.
    # Inativa: apaga, como antes (libera ~10 GB).
    if [ "${params.assembly_retry != false ? 'keep' : 'drop'}" = "keep" ] && [ ${max_reads} -gt 0 ]; then
        mv ${accession}/${accession}.sra ./${accession}.sra
    fi
    rm -rf ${accession}
    """
}
