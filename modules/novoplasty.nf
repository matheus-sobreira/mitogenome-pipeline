/*
 * Módulo: NOVOPLASTY
 * Monta o genoma mitocondrial usando a estratégia seed-and-extend.
 *
 * Entradas : tuple (sample_id, read1.fq, read2.fq)  ← leituras trimadas
 *            path(seed.fasta)                         ← sequência semente
 *
 * Saídas   : arquivos FASTA da montagem
 *            log de execução
 *            arquivo de configuração gerado
 *
 * Saídas possíveis do NOVOPlasty:
 *   Circularized_assembly_1_<sample>.fasta  → montagem completa circular ✓
 *   Contigs_1_<sample>.fasta                → contigs parciais (re-rodar)
 *   Option_<N>_<sample>.fasta               → montagens alternativas
 */

process NOVOPLASTY {

    tag "${sample_id}"

    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path  seed

    output:
    path "*.fasta",    emit: assembly, optional: true
    path "log_*.txt",  emit: log, optional: true
    path "config.txt", emit: config

    script:
    """
    # Gera o arquivo de configuração do NOVOPlasty dinamicamente
    cat > config.txt << EOF
Project:
-----------------------
Project name          = ${sample_id}
Type                  = mito
Genome range          = ${params.genome_range}
K-mer                 = ${params.k_mer}
Max memory            = ${params.max_memory}
Extended log          = 0
Save assembled reads  = no
Seed Input            = ${seed}
Extend seed directly  = no
Reference sequence    =
Variance detection    = no
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = ${params.read_length}
Insert size           = ${params.insert_size}
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = ${read1}
Reverse reads         = ${read2}
Store Hash            = no

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PATPRO                =
EOF

    # Executa o NOVOPlasty
    perl /opt/novoplasty/NOVOPlasty4.3.pl -c config.txt
    """
}
