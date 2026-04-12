#!/usr/bin/env nextflow
// Teste isolado do módulo MITOS2

nextflow.enable.dsl = 2

include { MITOS2 } from './modules/mitos2'

workflow {
    assembly_ch = Channel.of(
        tuple('SRR28399504', file("${projectDir}/results/a_leari/assembly/Circularized_assembly_1_SRR28399504.fasta"))
    )
    MITOS2(assembly_ch, file(params.mitos2_db))
}
