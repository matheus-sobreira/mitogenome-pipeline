/*
 * Módulo: RESAMPLE_POOL
 * Re-sorteia o pool de reads a partir do .sra JÁ BAIXADO, com semente nova.
 *
 * É o mecanismo de tratamento da rearbitragem (DEC-22): quando a validação
 * pós-montagem acusa anomalia de repeat, o número de cópias escolhido pelo
 * montador foi um sorteio de caminho — e a resposta é sortear de novo com
 * outro pool e deixar o detector arbitrar. Janelas com jitter de semente
 * diferente extraem SPOTS DIFERENTES do run; nenhuma requisição de rede.
 *
 * Entradas : tuple (accession, sra), max_reads, script de amostragem, semente
 * Saídas   : tuple (accession, R1, R2) — pool novo, mesmo volume
 */

process RESAMPLE_POOL {

    tag "${accession}"

    publishDir "${params.outdir}/qc/download_retry1", mode: 'copy', overwrite: true,
               pattern: 'sampling_plan.txt'

    input:
    tuple val(accession), path(sra)
    val max_reads
    path sampler
    val retry_seed

    output:
    tuple val(accession),
          path("${accession}_1.fastq"),
          path("${accession}_2.fastq"),
          emit: reads
    path "sampling_plan.txt", optional: true, emit: plan

    script:
    def mode    = params.sra_sampling ?: 'stratified'
    def windows = params.sra_windows  ?: 50

    """
    # Habilita SDL2 — o fastq-dump toca metadados mesmo com .sra local
    mkdir -p \$HOME/.ncbi
    cat > \$HOME/.ncbi/user-settings.mkfg << 'MKFG'
/LIBS/GUID = "d3b07384-d9a3-4f3c-b7e2-4c5c1d3abe42"
/repository/remote/main/SDL.2/resolver-cgis/enabled = "true"
/repository/user/cache-disabled = "true"
MKFG

    echo "[RESAMPLE] Pool re-sorteado: ${max_reads} spots, semente ${retry_seed} (DEC-22)"
    bash ${sampler} \\
        ${sra} \\
        ${accession} \\
        ${max_reads} \\
        ${mode} \\
        ${windows} \\
        ${retry_seed}
    """
}
