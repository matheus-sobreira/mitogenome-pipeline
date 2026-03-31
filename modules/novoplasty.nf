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
    def kmers = params.novoplasty_kmers ?: "${params.k_mer}"
    def max_iter = params.novoplasty_max_iterations ?: 3

    """
    #!/bin/bash
    set -euo pipefail

    ORIGINAL_SEED="${seed}"
    CIRCULARIZED=false
    MAX_ITER=${max_iter}

    for KMER in \$(echo "${kmers}" | tr ',' ' '); do

        SEED="\$ORIGINAL_SEED"

        for ITER in \$(seq 1 \$MAX_ITER); do

            echo "============================================"
            echo "=== k-mer = \$KMER | iteração \$ITER/\$MAX_ITER ==="
            echo "=== Seed: \$SEED"
            echo "============================================"

            cat > config.txt << EOF
Project:
-----------------------
Project name          = ${sample_id}
Type                  = mito
Genome Range          = ${params.genome_range}
K-mer                 = \$KMER
Max memory            = ${params.max_memory}
Extended log          = 0
Save assembled reads  = no
Seed Input            = \$SEED
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
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =
EOF

            perl /opt/novoplasty/NOVOPlasty4.3.1.pl -c config.txt

            # Circularizou?
            if find . -maxdepth 1 -name 'Circularized_assembly_*${sample_id}.fasta' | grep -q .; then
                echo "=== SUCESSO: circularizado com k-mer \$KMER na iteração \$ITER ==="
                CIRCULARIZED=true
                break 2
            fi

            # Extrai o maior contig como nova semente
            CONTIGS="\$(find . -maxdepth 1 \\( -name 'Contigs_*${sample_id}.fasta' -o -name 'Merged_contigs_*${sample_id}.fasta' \\) | head -1)"
            if [ -z "\$CONTIGS" ]; then
                echo "=== Nenhum contig gerado, pulando para próximo k-mer ==="
                break
            fi

            awk '/^>/{if(seq && length(seq)>max){max=length(seq);hdr=h;best=seq} h=\$0;seq=""} !/^>/{seq=seq\$0} END{if(length(seq)>max){hdr=\$0;best=seq} print hdr; print best}' "\$CONTIGS" > new_seed.fasta
            NEW_SIZE=\$(grep -v '^>' new_seed.fasta | wc -c)
            echo "=== Maior contig: \$NEW_SIZE bp ==="

            SEED="new_seed.fasta"

            # Limpa outputs parciais
            find . -maxdepth 1 \\( -name 'Contigs_*${sample_id}.fasta' -o -name 'Merged_contigs_*${sample_id}.fasta' -o -name 'log_*${sample_id}.txt' \\) -delete

        done

        if [ "\$CIRCULARIZED" = true ]; then
            break
        fi

        echo "=== k-mer \$KMER esgotou \$MAX_ITER iterações sem circularizar ==="

    done

    if [ "\$CIRCULARIZED" = false ]; then
        echo "=== AVISO: nenhuma combinação de k-mer/iteração circularizou a montagem ==="
    fi

    # Libera espaço: apaga reads trimados originais (resolve symlinks)
    for f in ${read1} ${read2}; do
        real=\$(readlink -f "\$f")
        if [ -f "\$real" ]; then
            rm -f "\$real"
            echo "Limpeza: apagado \$real"
        fi
    done
    """
}
