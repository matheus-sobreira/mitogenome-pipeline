#!/usr/bin/env nextflow

/*
 * Pipeline simplificado para montagem de mitogenomas
 * Fase 1 (Validação): SRA-Toolkit → FastQC → Trim Galore → NOVOPlasty
 *
 * Universidade do Estado do Rio Grande do Norte (UERN)
 * Autor: Matheus Sobreira Benevides
 */

include { SRA_PILOT_SAMPLE    } from './modules/sra_pilot'
include { PILOT_QC            } from './modules/pilot_qc'
include { ASSEMBLY_VALIDATION } from './modules/assembly_validation'
// Cadeia de rearbitragem (DEC-22) — os mesmos processos, aliás distintos para
// publishDir/config próprios. Nextflow não tem laço; a cadeia é desenrolada.
include { RESAMPLE_POOL                                      } from './modules/resample'
include { TRIM_GALORE         as TRIM_GALORE_RETRY           } from './modules/trim_galore'
include { NOVOPLASTY          as NOVOPLASTY_RETRY            } from './modules/novoplasty'
include { ASSEMBLY_VALIDATION as ASSEMBLY_VALIDATION_RETRY   } from './modules/assembly_validation'
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

/*
 * Resolve UMA VEZ a memória que o NOVOPlasty receberá.
 *
 * Antes isso era auto-detectado dentro do próprio processo NOVOPLASTY, o que
 * bastava enquanto ninguém mais precisava do valor. Agora o PILOT_QC também
 * precisa: a janela de intake do montador (≈1,691 M reads por GB) é o que
 * limita o volume útil de download, então o piloto tem de saber com quanta
 * memória o montador vai rodar. Duas cópias da mesma regra de auto-detecção
 * divergiriam; um valor só, calculado aqui e passado aos dois, não. Ver DEC-13.
 *
 * Limitação conhecida: lê o /proc/meminfo da máquina que lança o workflow. Em
 * executor distribuído (HPC), o nó de execução pode ter memória diferente —
 * nesse caso, informar params.max_memory explicitamente.
 */
def resolveNovoplastyMemory() {
    if (params.max_memory != null) {
        return params.max_memory as int
    }
    try {
        def line = new File('/proc/meminfo').readLines().find { it.startsWith('MemTotal') }
        def totalGb = (line.replaceAll(/[^0-9]/, '') as long) / 1024.0 / 1024.0
        def resolved = Math.min((totalGb * 0.5) as int, (totalGb - 4) as int)
        return Math.max(resolved, 2)
    }
    catch (Exception e) {
        log.warn "Não foi possível auto-detectar a RAM (${e.message}); usando 4 GB. " +
                 "Informe --max_memory para controlar."
        return 4
    }
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
    Amostragem  : piloto '${params.pilot_sampling}' (${params.pilot_windows} janelas) · download '${params.sra_sampling}' (${params.sra_windows} janelas) · semente ${params.pilot_seed}
    Read length : ${params.sra_max_reads != null ? params.read_length + ' bp (config)' : '(automático — Pilot QC)'}
    Max memory  : ${params.max_memory != null ? params.max_memory + ' GB (config)' : '(automático)'}
    BD MITOS2   : ${params.mitos2_db ?: '(não informado — etapa de anotação pulada)'}
    Organismo   : ${params.organism ?: '(não informado)'}
    """.stripIndent()

    // Canal com o acesso SRA
    accession_ch = Channel.value(params.sra_accession)

    // Arquivo semente para o NOVOPlasty
    seed_file = file(params.seed)

    // Scripts auxiliares entram como INPUT dos processos, não via ${projectDir}:
    // o Nextflow monta apenas o diretório da task no container, então um caminho
    // do projeto não existe lá dentro. Como input, o script é estagiado — e seu
    // conteúdo passa a fazer parte da chave de cache, de modo que editá-lo
    // invalida o -resume em vez de reaproveitar resultado obsoleto.
    sampler_script = file("${projectDir}/scripts/sra_sample.sh")
    qc_script      = file("${projectDir}/scripts/pilot_qc.sh")

    // Um único valor de memória para PILOT_QC e NOVOPLASTY (ver acima)
    novoplasty_mem = resolveNovoplastyMemory()
    log.info "Memória do NOVOPlasty: ${novoplasty_mem} GB " +
             "${params.max_memory != null ? '(configurada)' : '(auto-detectada)'} " +
             "→ janela de intake ≈ ${Math.round(1691018 * novoplasty_mem / 1e6)}M reads"

    // ── Etapa 0: Pilot QC (automático se sra_max_reads não definido) ──
    // Amostra uma fração pequena do run (500K reads, estratificada ao longo
    // dele) e analisa qualidade, fração mitocondrial e adaptadores para
    // recomendar o max_reads ideal. Amostragem e análise são processos
    // separados: exigem imagens Docker distintas (ver DEC-03).
    // A validação pós-montagem (DEC-20) reutiliza a amostra piloto, então ela
    // é coletada MESMO quando sra_max_reads é forçado — o que se pula nesse
    // caso é só a ANÁLISE (PILOT_QC), não a coleta.
    def validationEnabled = (params.assembly_validation != false)

    if (params.sra_max_reads != null) {
        log.info "Pilot QC: análise pulada (sra_max_reads = ${params.sra_max_reads})"
        max_reads_ch   = Channel.value(params.sra_max_reads)
        read_length_ch = Channel.value(params.read_length)
        if (validationEnabled) {
            log.info "Amostra piloto: coletada mesmo assim — a validação pós-montagem a reutiliza (DEC-20)"
            SRA_PILOT_SAMPLE(accession_ch, sampler_script)
        }
    } else {
        log.info "Pilot QC: ativado — ${params.pilot_reads ?: 500000} reads piloto, " +
                 "amostragem '${params.pilot_sampling ?: 'stratified'}'"
        // Referência do mapeamento: mitogenoma completo quando informado,
        // senão a própria seed com escala explícita (DEC-06)
        pilot_ref = params.pilot_reference ? file(params.pilot_reference) : seed_file
        if (params.pilot_reference) {
            log.info "Pilot QC: fração mitocondrial medida contra ${params.pilot_reference}"
        } else {
            log.warn "Pilot QC: sem --pilot_reference; medindo contra a seed e escalando " +
                     "para ${params.genome_range} bp (assume cobertura uniforme)"
        }

        SRA_PILOT_SAMPLE(accession_ch, sampler_script)
        PILOT_QC(SRA_PILOT_SAMPLE.out.reads, pilot_ref, qc_script, novoplasty_mem)
        max_reads_ch   = PILOT_QC.out.recommended_reads
        read_length_ch = PILOT_QC.out.read_length
    }

    // ── Etapa 1: Download das leituras via SRA-Toolkit ──────────────
    SRA_DOWNLOAD(accession_ch, max_reads_ch, sampler_script)

    // ── Etapa 2: Controle de qualidade das leituras brutas ──────────
    FASTQC(SRA_DOWNLOAD.out.reads)

    // ── Etapa 3: Remoção de adaptadores e trimming ──────────────────
    TRIM_GALORE(SRA_DOWNLOAD.out.reads)

    // ── Etapa 4: Montagem do mitogenoma ─────────────────────────────
    NOVOPLASTY(TRIM_GALORE.out.reads, seed_file, read_length_ch, novoplasty_mem)

    // Filtra apenas montagens circularizadas do NOVOPlasty
    // (usado pela validação pós-montagem E pela anotação)
    circularized_ch = NOVOPLASTY.out.assembly
        .flatten()
        .filter { fasta -> fasta.name.startsWith('Circularized_assembly') }
        .map { fasta ->
            def sample_id = fasta.baseName.replaceAll(/^Circularized_assembly_\d+_/, '')
            tuple(sample_id, fasta)
        }

    // ── Etapa 4b: Validação pós-montagem (DEC-20) ───────────────────
    // Mapeia a amostra piloto de volta contra a montagem: razão de
    // profundidade << 1 = repeat superexpandido; >> 1 = repeat colapsado.
    // "Circular ≠ correto": em 14/08/2026 uma montagem de A. leari
    // circularizou com um VNTR da região controle superexpandido em ~1 kb —
    // dentro do genome_range, invisível à checagem de comprimento.
    // Canal padrão: sem validação (ou sem rearbitragem), a montagem segue
    // direto — comportamento anterior preservado.
    final_assembly_ch = circularized_ch

    if (validationEnabled) {
        validation_ref = params.pilot_reference ? file(params.pilot_reference) : seed_file
        validation_script = file("${projectDir}/scripts/assembly_validation.sh")
        ASSEMBLY_VALIDATION(
            circularized_ch.join(SRA_PILOT_SAMPLE.out.reads),
            validation_script,
            validation_ref
        )

        // ── Etapa 4c: Rearbitragem (DEC-22) ─────────────────────────
        // ALERTA de repeat = o número de cópias foi sorteio de caminho do
        // montador. Tratamento: re-sortear o pool (semente nova, .sra local,
        // zero rede), re-montar e deixar o detector arbitrar. Uma tentativa
        // desenrolada; a montagem da retentativa segue adiante qualquer que
        // seja seu veredito (última instância — o status final o reporta).
        if (params.assembly_retry != false) {
            verdict0 = ASSEMBLY_VALIDATION.out.verdict
                .map { sid, f -> tuple(sid, f.text.trim()) }

            asm_ok0 = circularized_ch.join(verdict0)
                .filter { sid, fasta, v -> v != 'ALERTA' }
                .map    { sid, fasta, v -> tuple(sid, fasta) }

            need_retry = verdict0
                .filter { sid, v -> v == 'ALERTA' }
                .map    { sid, v -> sid }

            retry_seed = (params.pilot_seed ?: 42) + (params.retry_seed_offset ?: 1000)

            RESAMPLE_POOL(
                need_retry.combine(SRA_DOWNLOAD.out.sra).map { sid, sra -> tuple(sid, sra) },
                max_reads_ch,
                sampler_script,
                retry_seed
            )
            TRIM_GALORE_RETRY(RESAMPLE_POOL.out.reads)
            NOVOPLASTY_RETRY(TRIM_GALORE_RETRY.out.reads, seed_file, read_length_ch, novoplasty_mem)

            circularized_retry_ch = NOVOPLASTY_RETRY.out.assembly
                .flatten()
                .filter { fasta -> fasta.name.startsWith('Circularized_assembly') }
                .map { fasta ->
                    def sample_id = fasta.baseName.replaceAll(/^Circularized_assembly_\d+_/, '')
                    tuple(sample_id, fasta)
                }

            ASSEMBLY_VALIDATION_RETRY(
                circularized_retry_ch.join(SRA_PILOT_SAMPLE.out.reads),
                validation_script,
                validation_ref
            )

            final_assembly_ch = asm_ok0.mix(circularized_retry_ch)
        }
    } else {
        log.warn "Validação pós-montagem: DESATIVADA (--assembly_validation false)"
    }

    // ── Etapa 5: Anotação funcional (MITOS2) ────────────────────────
    // Executada apenas quando --mitos2_db for informado no perfil ou CLI.
    // O banco RefSeq89 deve ser baixado uma vez:
    //   wget https://zenodo.org/record/4284483/files/mitos2-refseq89m.db.gz
    //   gunzip mitos2-refseq89m.db.gz
    if (params.mitos2_db) {

        // Anota a montagem ARBITRADA: a aprovada na primeira validação, ou a
        // da retentativa (DEC-22) — nunca uma reprovada que tinha substituta.
        MITOS2(final_assembly_ch, file(params.mitos2_db))

        // ── Etapa 6: Compilar entregáveis ─────────────────────────────
        // Combina montagem + anotação MITOS2 em pasta organizada
        summary_ch = final_assembly_ch
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

    // Vereditos da validação pós-montagem (DEC-20) e da rearbitragem (DEC-22)
    def valReport = file("${params.outdir}/qc/assembly_validation/validation_report.txt")
    def valLine = valReport.exists() ?
        valReport.readLines().find { it.startsWith('VEREDITO') } : null
    def retryReport = file("${params.outdir}/qc/assembly_validation_retry1/validation_report.txt")
    def retryLine = retryReport.exists() ?
        retryReport.readLines().find { it.startsWith('VEREDITO') } : null

    // O veredito EFETIVO é o da última instância que rodou
    def effectiveLine = retryLine ?: valLine

    def status
    if (!workflow.success) {
        status = 'FALHOU ✗'
    } else if (!hasCircularized) {
        status = 'INCOMPLETO ⚠ — montagem não circularizou'
    } else if (effectiveLine != null && effectiveLine.contains('ALERTA')) {
        status = (retryLine != null) ?
            'SUCESSO com ALERTA ⚠ — anomalia persistiu após rearbitragem' :
            'SUCESSO com ALERTA ⚠ — circularizou, mas a validação sinalizou anomalia'
    } else if (retryLine != null) {
        status = 'SUCESSO ✓ (via rearbitragem — 1ª montagem descartada por anomalia)'
    } else {
        status = 'SUCESSO ✓'
    }

    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║                 PIPELINE CONCLUÍDO                       ║
    ╚══════════════════════════════════════════════════════════╝
    Status    : ${status}
    Validação : ${valLine ?: '(não executada)'}${retryLine != null ? "\n    Rearbitr. : ${retryLine}" : ''}
    Duração   : ${workflow.duration}
    Saída     : ${params.outdir}
    """.stripIndent()

    if (effectiveLine != null && effectiveLine.contains('ALERTA')) {
        log.warn "A montagem circularizou mas tem região(ões) com profundidade anômala —"
        log.warn "circular ≠ correto. Detalhes em ${params.outdir}/qc/assembly_validation*/"
    }

    if (!hasCircularized && workflow.success) {
        log.warn "O NOVOPlasty não conseguiu circularizar o genoma mitocondrial."
        log.warn "Possíveis causas:"
        log.warn "  • Cobertura mitocondrial insuficiente no dataset"
        log.warn "  • Seed muito divergente da espécie-alvo"
        log.warn "  • Dataset não é WGS (ex: exoma, RAD-seq, RNA-seq)"
        log.warn "Verifique o log do NOVOPlasty em: ${params.outdir}/assembly/"
    }
}
