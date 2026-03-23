#!/usr/bin/env nextflow

/*
 * Pipeline simplificado para montagem de mitogenomas
 * Fase 1 (Validação): SRA-Toolkit → FastQC → Trim Galore → NOVOPlasty
 *
 * Universidade do Estado do Rio Grande do Norte (UERN)
 * Autor: Matheus Sobreira Benevides
 */

nextflow.enable.dsl = 2

include { SRA_DOWNLOAD } from './modules/sra_download'
include { FASTQC       } from './modules/fastqc'
include { TRIM_GALORE  } from './modules/trim_galore'
include { NOVOPLASTY   } from './modules/novoplasty'

// Validação dos parâmetros obrigatórios
if (!params.sra_accession) {
    error "ERROR: Informe o acesso SRA com: --sra_accession <ACCESSION>\n" +
          "       Exemplo: --sra_accession SRR2081280"
}
if (!params.seed) {
    error "ERROR: Informe o arquivo semente com: --seed <caminho/para/seed.fasta>\n" +
          "       Veja: data/seeds/COMO_OBTER_SEMENTE.md"
}

workflow {

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║         MITOGENOME PIPELINE — Fase 1 (Validação)         ║
    ╚══════════════════════════════════════════════════════════╝
    Acesso SRA  : ${params.sra_accession}
    Semente     : ${params.seed}
    Saída       : ${params.outdir}
    Range mtDNA : ${params.genome_range} bp
    """.stripIndent()

    // Canal com o acesso SRA
    accession_ch = Channel.value(params.sra_accession)

    // Arquivo semente para o NOVOPlasty
    seed_file = file(params.seed)

    // ── Etapa 1: Download das leituras via SRA-Toolkit ──────────────
    SRA_DOWNLOAD(accession_ch)

    // ── Etapa 2: Controle de qualidade das leituras brutas ──────────
    FASTQC(SRA_DOWNLOAD.out.reads)

    // ── Etapa 3: Remoção de adaptadores e trimming ──────────────────
    TRIM_GALORE(SRA_DOWNLOAD.out.reads)

    // ── Etapa 4: Montagem do mitogenoma ─────────────────────────────
    NOVOPLASTY(TRIM_GALORE.out.reads, seed_file)
}

workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║                 PIPELINE CONCLUÍDO                       ║
    ╚══════════════════════════════════════════════════════════╝
    Status  : ${workflow.success ? 'SUCESSO ✓' : 'FALHOU ✗'}
    Duração : ${workflow.duration}
    Saída   : ${params.outdir}
    """.stripIndent()
}
