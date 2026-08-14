/*
 * Módulo: SRA_PILOT_SAMPLE
 * Baixa a amostra piloto de reads que alimenta a análise do PILOT_QC.
 *
 * Amostragem estratificada por padrão (params.pilot_sampling): em vez dos N
 * primeiros spots do run, extrai N spots distribuídos em janelas equidistantes
 * com deslocamento determinístico. Ver scripts/sra_sample.sh e DEC-01 no vault.
 *
 * A análise vive em módulo separado (PILOT_QC) porque exige bwa/samtools, que
 * não estão nesta imagem — e porque o split permite iterar na análise sem
 * re-baixar a amostra (ver DEC-03).
 *
 * Entradas : val(accession)
 * Saídas   : tuple(accession, R1, R2) — amostra piloto
 *            path sampling_plan.txt   — janelas usadas, para reprodutibilidade
 */

process SRA_PILOT_SAMPLE {

    tag "${accession}"

    publishDir "${params.outdir}/qc/pilot", mode: 'copy', overwrite: true,
               pattern: 'sampling_plan.txt'

    input:
    val accession
    path sampler

    output:
    tuple val(accession),
          path("${accession}_1.fastq"),
          path("${accession}_2.fastq"),
          emit: reads
    path "sampling_plan.txt", emit: plan

    script:
    def pilot_n  = params.pilot_reads ?: 500000
    def mode     = params.pilot_sampling ?: 'stratified'
    def windows  = mode == 'dense' ? (params.pilot_dense_windows ?: 500)
                                   : (params.pilot_windows ?: 10)
    def seed_rng = params.pilot_seed ?: 42

    """
    # ── Configuração do SRA-Toolkit ──────────────────────────────────────
    mkdir -p \$HOME/.ncbi
    cat > \$HOME/.ncbi/user-settings.mkfg << 'MKFG'
/LIBS/GUID = "d3b07384-d9a3-4f3c-b7e2-4c5c1d3abe42"
/repository/remote/main/SDL.2/resolver-cgis/enabled = "true"
/repository/user/cache-disabled = "true"
MKFG

    echo "[PILOT] Amostra piloto: ${pilot_n} reads, modo '${mode}'"

    bash ${sampler} \\
        ${accession} \\
        ${accession} \\
        ${pilot_n} \\
        ${mode} \\
        ${windows} \\
        ${seed_rng}
    """
}
