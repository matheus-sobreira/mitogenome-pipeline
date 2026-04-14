/*
 * Módulo: SRA_PILOT
 * Baixa uma amostra piloto de reads e analisa a qualidade para
 * determinar automaticamente o número ideal de reads (max_reads).
 *
 * Usa fastq-dump -X, que baixa APENAS os N primeiros reads da origem,
 * sem transferir o dataset completo (diferente de fasterq-dump).
 *
 * Entradas : val(accession), path(seed)
 * Saídas   : env RECOMMENDED_READS  — número recomendado de reads
 *            path pilot_report.txt  — relatório detalhado da análise
 */

process SRA_PILOT {

    tag "${accession}"

    publishDir "${params.outdir}/qc/pilot", mode: 'copy', overwrite: true

    input:
    val accession
    path seed

    output:
    env RECOMMENDED_READS, emit: recommended_reads
    env READ_LENGTH,       emit: read_length
    path "pilot_report.txt", emit: report

    script:
    def pilot_n      = params.pilot_reads ?: 500000
    def genome_range = params.genome_range
    def target_cov   = params.target_coverage ?: 500
    def mito_frac    = params.mito_fraction ?: 0

    """
    # ── Configuração do SRA-Toolkit ──────────────────────────────────────
    mkdir -p \$HOME/.ncbi
    cat > \$HOME/.ncbi/user-settings.mkfg << 'MKFG'
/LIBS/GUID = "d3b07384-d9a3-4f3c-b7e2-4c5c1d3abe42"
/repository/remote/main/SDL.2/resolver-cgis/enabled = "true"
/repository/user/cache-disabled = "true"
MKFG

    echo "[PILOT] Baixando ${pilot_n} reads de amostra piloto..."

    # fastq-dump -X limita NA ORIGEM (baixa apenas o necessário)
    fastq-dump \\
        -X ${pilot_n} \\
        --split-files \\
        --outdir . \\
        ${accession}

    # ── Análise de qualidade ─────────────────────────────────────────────
    bash ${projectDir}/scripts/pilot_qc.sh \\
        ${accession}_1.fastq \\
        ${accession}_2.fastq \\
        ${seed} \\
        "${genome_range}" \\
        ${target_cov} \\
        ${mito_frac}

    # Captura para output env do Nextflow
    RECOMMENDED_READS=\$(cat recommended_reads.txt)
    READ_LENGTH=\$(cat read_length.txt)

    # Limpeza — pilot reads não são necessários no futuro
    rm -f ${accession}_1.fastq ${accession}_2.fastq
    """
}
