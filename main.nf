#!/usr/bin/env nextflow

/*
 * Pipeline simplificado para montagem de mitogenomas
 * Fase 1 (Validação): SRA-Toolkit → FastQC → Trim Galore → NOVOPlasty
 *
 * Universidade do Estado do Rio Grande do Norte (UERN)
 * Autor: Matheus Sobreira Benevides
 */

include { SRA_DOWNLOAD } from './modules/sra_download'
include { FASTQC       } from './modules/fastqc'
include { TRIM_GALORE  } from './modules/trim_galore'
include { NOVOPLASTY   } from './modules/novoplasty'
include { MITOS2       } from './modules/mitos2'

// Validação dos parâmetros obrigatórios
if (!params.sra_accession) {
    error "ERROR: Informe o acesso SRA com: --sra_accession <ACCESSION>\n" +
          "       Exemplo: --sra_accession SRR36152783"
}
if (!params.seed) {
    error "ERROR: Informe o arquivo semente com: --seed <caminho/para/seed.fasta>\n" +
          "       Veja: data/seeds/COMO_OBTER_SEMENTE.md"
}

workflow {

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║            MITOGENOME PIPELINE — Fase 1 + 2              ║
    ╚══════════════════════════════════════════════════════════╝
    Acesso SRA  : ${params.sra_accession}
    Semente     : ${params.seed}
    Saída       : ${params.outdir}
    Range mtDNA : ${params.genome_range} bp
    Código gen. : ${params.genetic_code}  (2=vertebrado; 5=invertebrado)
    BD MITOS2   : ${params.mitos2_db ?: '(não informado — etapa de anotação pulada)'}
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

    // ── Etapa 5: Anotação funcional (MITOS2) ────────────────────────
    // Executada apenas quando --mitos2_db for informado no perfil ou CLI.
    // O banco RefSeq89 deve ser baixado uma vez:
    //   wget https://zenodo.org/record/4284483/files/mitos2-refseq89m.db.gz
    //   gunzip mitos2-refseq89m.db.gz
    if (params.mitos2_db) {

        // Filtra apenas montagens circularizadas do NOVOPlasty
        circularized_ch = NOVOPLASTY.out.assembly
            .flatten()
            .filter { fasta -> fasta.name.startsWith('Circularized_assembly') }
            .map { fasta ->
                def sample_id = fasta.baseName.replaceAll(/^Circularized_assembly_\d+_/, '')
                tuple(sample_id, fasta)
            }

        MITOS2(circularized_ch, file(params.mitos2_db))

    } else {
        log.warn "MITOS2: parâmetro 'mitos2_db' não informado. Etapa de anotação pulada."
        log.warn "        Defina em conf/<perfil>.config ou passe --mitos2_db <caminho> na linha de comando."
    }
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
