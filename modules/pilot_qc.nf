/*
 * Módulo: PILOT_QC
 * Analisa a amostra piloto e recomenda o número ideal de reads (max_reads).
 *
 * Roda na imagem pilot-qc (bwa + samtools + seqtk). A fração mitocondrial é
 * MEDIDA por mapeamento, não estimada por k-mers — ver DEC-04/DEC-05 no vault
 * e o cabeçalho de scripts/pilot_qc.sh.
 *
 * A referência do mapeamento é params.pilot_reference (mitogenoma completo)
 * quando informada; caso contrário a própria seed, com escala explícita e a
 * premissa declarada no relatório (DEC-06).
 *
 * Entradas : tuple(accession, R1, R2), path(reference)
 * Saídas   : env RECOMMENDED_READS  — número recomendado de reads
 *            env READ_LENGTH        — comprimento médio observado
 *            path pilot_report.txt  — relatório detalhado da análise
 */

process PILOT_QC {

    tag "${accession}"

    publishDir "${params.outdir}/qc/pilot", mode: 'copy', overwrite: true,
               pattern: 'pilot_report.txt'

    input:
    tuple val(accession), path(read1), path(read2)
    path reference
    path qc_script
    val  novoplasty_mem

    output:
    env RECOMMENDED_READS, emit: recommended_reads
    env READ_LENGTH,       emit: read_length
    path "pilot_report.txt", emit: report

    script:
    def genome_range = params.genome_range
    def target_cov   = params.target_coverage ?: 500
    def mito_frac    = params.mito_fraction ?: 0
    def ref_mode     = params.pilot_reference ? 'genome' : 'seed'
    def bwa_opts     = params.pilot_bwa_opts ?: ''

    """
    # Teto de segurança, reportado quando aplicado (DEC-10)
    export PILOT_MAX_READS_CAP=${params.pilot_max_reads_cap ?: 25000000}

    # Constantes do montador, medidas em 13/08/2026 (DEC-13). A memória vem
    # resolvida do main.nf para que piloto e montador usem o MESMO valor.
    export NOVOPLASTY_MEM_GB=${novoplasty_mem}
    export NOVOPLASTY_READS_PER_GB=${params.novoplasty_reads_per_gb ?: 1691018}
    export NOVOPLASTY_RETENTION=${params.novoplasty_retention ?: 0.85}
    export TRIM_LOSS=${params.trim_loss_estimate ?: 0.10}

    bash ${qc_script} \\
        ${read1} \\
        ${read2} \\
        ${reference} \\
        "${genome_range}" \\
        ${target_cov} \\
        ${mito_frac} \\
        ${ref_mode} \\
        ${task.cpus} \\
        "${bwa_opts}"

    # Captura para output env do Nextflow
    RECOMMENDED_READS=\$(cat recommended_reads.txt)
    READ_LENGTH=\$(cat read_length.txt)
    """
}
