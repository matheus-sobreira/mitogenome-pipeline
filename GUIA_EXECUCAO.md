# Guia de Execução — mitogenome-pipeline

## Pré-requisitos

- Docker (`docker info`)
- Nextflow ≥ 23.x (`nextflow -version`)
- Java 17+ (`java -version`)
- ~50 GB de disco livre

---

## 1. Preparação (apenas uma vez)

### Construir imagens Docker

```bash
for dir in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker build -t "mitogenome-pipeline/${dir}:1.0" "docker/${dir}"
done
```

### Baixar banco de dados MITOS2

```bash
mkdir -p data/databases
wget https://zenodo.org/record/4284483/files/refseq89m.tar.bz2
tar -xjf refseq89m.tar.bz2 -C data/databases/
rm refseq89m.tar.bz2
```

---

## 2. Executar o pipeline

```bash
# Arara-azul-de-lear (estudo principal)
./run_pipeline.sh -profile a_leari

# D. bifasciatus (validação — sem anotação MITOS2)
./run_pipeline.sh -profile test
```

O `run_pipeline.sh` executa em background com `nohup` — o terminal pode ser fechado sem interromper. O log fica em `logs/`.

Para acompanhar:

```bash
tail -f logs/a_leari_*.log
```

Alternativamente, para execução direta no terminal:

```bash
nextflow run main.nf -profile a_leari
```

---

## 3. Retomar execução interrompida

```bash
# Listar runs anteriores
ls work/a_leari/

# Retomar um run específico
nextflow run main.nf -profile a_leari -resume -w work/a_leari/2026-04-12_14h30
```

---

## 4. Forçar número de reads (pular Pilot QC)

```bash
nextflow run main.nf -profile a_leari --sra_max_reads 20000000
```

Se `sra_max_reads` não for definido, o Pilot QC roda automaticamente e calcula o valor ideal.

---

## 5. Fluxo do pipeline

```
SRA_PILOT → SRA_DOWNLOAD → FASTQC → TRIM_GALORE → NOVOPLASTY → MITOS2 → COMPILE_SUMMARY
```

| Etapa | Função |
|---|---|
| SRA_PILOT | Analisa 500K reads e calcula o volume ideal (opcional) |
| SRA_DOWNLOAD | Baixa e trunca reads do NCBI SRA |
| FASTQC | Controle de qualidade |
| TRIM_GALORE | Remoção de adaptadores e bases de baixa qualidade |
| NOVOPLASTY | Montagem *de novo* por seed-and-extend |
| MITOS2 | Anotação funcional (13 CDS, 22 tRNA, 2 rRNA) |
| COMPILE_SUMMARY | Compila 14 categorias de entregáveis |

---

## 6. Saídas

Os resultados ficam em `results/<espécie>/summary/`:

| Arquivo | Conteúdo |
|---|---|
| `01_genome_assembly.fasta` | Mitogenoma circularizado |
| `02_coding_genes_nt.fasta` | 13 CDS (nucleotídeos) |
| `03_coding_genes_aa.fasta` | 13 CDS (aminoácidos) |
| `04_ribosomal_genes.fasta` | 12S + 16S rRNA |
| `05_transport_genes.fasta` | 22 tRNAs |
| `06_gene_positions.tsv` | Posição dos genes |
| `07_start_stop_codons.tsv` | Start/Stop codons |
| `08_trna_anticodons.tsv` | Anticódons tRNA |
| `09_circular_map.svg/pdf` | Mapa circular |
| `10_annotation.gff` | Anotação GFF3 |
| `11_structure_svgs/` | Estruturas secundárias tRNA/rRNA (SVG) |
| `12_quality_plots/` | Plots de qualidade MITOS2 |
| `13_gene_order.txt` | Ordem gênica |
| `genbank_submission/` | `.gbk` + `.tbl` + `.fsa` para GenBank |

---

## 7. Adicionar nova espécie

1. Buscar no [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) um dataset WGS paired-end
2. Obter o gene COX1 de uma espécie próxima como semente (ver `data/seeds/COMO_OBTER_SEMENTE.md`)
3. Criar `conf/nova_especie.config`:

```groovy
params {
    sra_accession    = 'SRRXXXXXXX'
    seed             = "${projectDir}/data/seeds/especie_proxima_cox1.fasta"
    outdir           = 'results/nova_especie'
    genome_range     = '15000-18500'
    novoplasty_kmers = '39,33'
    mitos2_db        = "${projectDir}/data/databases/refseq89m"
    genetic_code     = 2        // 2 = vertebrado, 5 = invertebrado
    organism         = 'Genus species'
}
```

4. Registrar em `nextflow.config` → seção `profiles`
5. Executar: `./run_pipeline.sh -profile nova_especie`

---

## 8. Limpeza

```bash
# Cache de trabalho de uma espécie
rm -rf work/a_leari/

# Arquivos temporários do Nextflow
rm -rf .nextflow/ .nextflow.log*
```
