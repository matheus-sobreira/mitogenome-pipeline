#!/usr/bin/env nextflow

/*
 * Pipeline simplificado para montagem de mitogenomas
 * Fase 1 (Validação): SRA-Toolkit → FastQC → Trim Galore → NOVOPlasty
 *
 * Universidade do Estado do Rio Grande do Norte (UERN)
 * Autor: Matheus Sobreira Benevides
 */

include { SRA_PILOT    } from './modules/sra_pilot'
include { SRA_DOWNLOAD } from './modules/sra_download'
include { FASTQC       } from './modules/fastqc'
include { TRIM_GALORE  } from './modules/trim_galore'
include { NOVOPLASTY   } from './modules/novoplasty'
include { MITOS2           } from './modules/mitos2'
include { COMPILE_SUMMARY  } from './modules/compile_summary'

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
    Max reads   : ${params.sra_max_reads != null ? params.sra_max_reads : '(automático — Pilot QC)'}
    Read length : ${params.sra_max_reads != null ? params.read_length + ' bp (config)' : '(automático — Pilot QC)'}
    Max memory  : ${params.max_memory != null ? params.max_memory + ' GB (config)' : '(automático)'}
    BD MITOS2   : ${params.mitos2_db ?: '(não informado — etapa de anotação pulada)'}
    Organismo   : ${params.organism ?: '(não informado)'}
    """.stripIndent()

    // Canal com o acesso SRA
    accession_ch = Channel.value(params.sra_accession)

    // Arquivo semente para o NOVOPlasty
    seed_file = file(params.seed)

    // ── Etapa 0: Pilot QC (automático se sra_max_reads não definido) ──
    // Baixa uma amostra pequena (500K reads) e analisa qualidade,
    // fração mitocondrial e adaptadores para recomendar o max_reads ideal.
    if (params.sra_max_reads != null) {
        log.info "Pilot QC: pulado (sra_max_reads = ${params.sra_max_reads})"
        max_reads_ch   = Channel.value(params.sra_max_reads)
        read_length_ch = Channel.value(params.read_length)
    } else {
        log.info "Pilot QC: ativado — analisando ${params.pilot_reads ?: 500000} reads piloto..."
        SRA_PILOT(accession_ch, seed_file)
        max_reads_ch   = SRA_PILOT.out.recommended_reads
        read_length_ch = SRA_PILOT.out.read_length
    }

    // ── Etapa 1: Download das leituras via SRA-Toolkit ──────────────
    SRA_DOWNLOAD(accession_ch, max_reads_ch)

    // ── Etapa 2: Controle de qualidade das leituras brutas ──────────
    FASTQC(SRA_DOWNLOAD.out.reads)

    // ── Etapa 3: Remoção de adaptadores e trimming ──────────────────
    TRIM_GALORE(SRA_DOWNLOAD.out.reads)

    // ── Etapa 4: Montagem do mitogenoma ─────────────────────────────
    NOVOPLASTY(TRIM_GALORE.out.reads, seed_file, read_length_ch)

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

        // ── Etapa 6: Compilar entregáveis ─────────────────────────────
        // Combina montagem + anotação MITOS2 em pasta organizada
        summary_ch = circularized_ch
            .join(MITOS2.out.annotation)
            .map { sample_id, assembly, mitos_dir ->
                tuple(sample_id, assembly, mitos_dir)
            }

        COMPILE_SUMMARY(summary_ch)

    } else {
        log.warn "MITOS2: parâmetro 'mitos2_db' não informado. Etapa de anotação pulada."
        log.warn "        Defina em conf/<perfil>.config ou passe --mitos2_db <caminho> na linha de comando."
    }
}

workflow.onComplete {
    // Verifica se houve montagem circularizada nos resultados
    def assemblyDir = file("${params.outdir}/assembly")
    def circularized = assemblyDir.exists() ?
        assemblyDir.listFiles()?.findAll { it.name.startsWith('Circularized_assembly') } : []
    def hasCircularized = circularized && !circularized.isEmpty()

    def status
    if (!workflow.success) {
        status = 'FALHOU ✗'
    } else if (!hasCircularized) {
        status = 'INCOMPLETO ⚠ — montagem não circularizou'
    } else {
        status = 'SUCESSO ✓'
    }

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║                 PIPELINE CONCLUÍDO                       ║
    ╚══════════════════════════════════════════════════════════╝
    Status  : ${status}
    Duração : ${workflow.duration}
    Saída   : ${params.outdir}
    """.stripIndent()

    if (!hasCircularized && workflow.success) {
        log.warn "O NOVOPlasty não conseguiu circularizar o genoma mitocondrial."
        log.warn "Possíveis causas:"
        log.warn "  • Cobertura mitocondrial insuficiente no dataset"
        log.warn "  • Seed muito divergente da espécie-alvo"
        log.warn "  • Dataset não é WGS (ex: exoma, RAD-seq, RNA-seq)"
        log.warn "Verifique o log do NOVOPlasty em: ${params.outdir}/assembly/"
    }
}
