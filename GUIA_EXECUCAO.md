# Guia de Execução — mitogenome-pipeline

## Pré-requisitos

- Docker (`docker info`)
- Nextflow ≥ 23.x (`nextflow -version`)
- Java 17+ (`java -version`)
- ~50 GB de disco livre

---

## 1. Preparação (apenas uma vez)

### Obter as imagens Docker

As 5 imagens do pipeline (`sra-tools`, `fastqc`, `trim-galore`, `novoplasty`, `mitos2`) estão publicadas no **GitHub Container Registry** com tag `:1.0` (imutável). O Nextflow puxa automaticamente na primeira execução; não é necessário fazer build local.

Para pré-baixar (opcional, evita atraso no primeiro run):

```bash
for tool in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker pull "ghcr.io/matheus-sobreira/mitogenome-pipeline/${tool}:1.0"
done
```

Se preferir buildar localmente (por exemplo, para modificar uma etapa), buildando **e retaggeando** para o namespace que o Nextflow espera:

```bash
for tool in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker build -t "ghcr.io/matheus-sobreira/mitogenome-pipeline/${tool}:1.0" "docker/${tool}"
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

O pipeline tem três perfis disponíveis:

```bash
# Arara-azul-de-lear (espécie-alvo do TCC)
./run_pipeline.sh -profile a_leari

# Diploprion bifasciatus (espécie-controle — verifica que o pipeline reproduz uma referência conhecida; anotação MITOS2 desativada por padrão)
./run_pipeline.sh -profile test

# Arara-azul-grande (validação cruzada recíproca; usa cox1 da montagem de A. leari como semente)
./run_pipeline.sh -profile a_hyacinthinus
```

### Como funciona o `run_pipeline.sh`

O script executa o Nextflow em **background via `nohup`** — o terminal pode ser fechado, a conexão SSH pode cair, e o pipeline continua rodando. Cada execução gera um log timestamped em `logs/<perfil>_<timestamp>.log` e grava o PID em `logs/.last_pid` para facilitar acompanhamento.

Para acompanhar a execução em tempo real:

```bash
tail -f logs/a_leari_*.log
```

Para verificar se o processo continua vivo:

```bash
kill -0 $(cat logs/.last_pid) && echo "rodando" || echo "terminou"
```

Alternativamente, para execução direta no terminal (sem nohup):

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
SRA_PILOT_SAMPLE → PILOT_QC → SRA_DOWNLOAD → FASTQC → TRIM_GALORE → NOVOPLASTY → MITOS2 → COMPILE_SUMMARY
```

| Etapa | Função |
|---|---|
| SRA_PILOT_SAMPLE | Amostra 500K reads distribuídos ao longo do run (opcional) |
| PILOT_QC | Analisa a amostra e calcula o volume ideal de download |
| SRA_DOWNLOAD | Baixa do NCBI SRA e amostra o volume definido |
| FASTQC | Controle de qualidade |
| TRIM_GALORE | Remoção de adaptadores e bases de baixa qualidade |
| NOVOPLASTY | Montagem *de novo* por seed-and-extend |
| MITOS2 | Anotação funcional (13 CDS, 22 tRNA, 2 rRNA) |
| COMPILE_SUMMARY | Compila 14 categorias de entregáveis |

---

## 6. Saídas

Os resultados ficam em `results/<espécie>/summary/deliverables/`:

| Arquivo / Diretório | Conteúdo |
|---|---|
| `01_genome_assembly.fasta` | Mitogenoma circularizado |
| `02_coding_genes_nt.fasta` | 13 CDS (nucleotídeos) |
| `03_coding_genes_aa.fasta` | 13 CDS (aminoácidos) |
| `04_ribosomal_genes.fasta` | 12S + 16S rRNA |
| `05_transport_genes.fasta` | 22 tRNAs |
| `06_gene_positions.tsv` | Posição dos genes |
| `07_start_stop_codons.tsv` | Start/Stop codons |
| `08_trna_anticodons.tsv` | Anticódons tRNA |
| `09_circular_map.svg` + `.pdf` | Mapa circular do mitogenoma (Biopython) |
| `09b_genome_map_linear.png` | Mapa linear do MITOS2 |
| `10_annotation.gff` | Anotação GFF3 completa |
| `13_gene_order.txt` | Ordem gênica linear |
| `structure_svgs/` | Estruturas secundárias de tRNA/rRNA em SVG (ViennaRNA/RNAplot) |
| `quality_plots/` | Plots de qualidade da anotação MITOS2 |
| `genbank_submission/` | `.gbk` + `.tbl` + `.fsa` prontos para submissão ao GenBank |

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
    genetic_code     = 2        // Código genético (transl_table do NCBI):
                                //   2 = mitocondrial vertebrado (aves, peixes, mamíferos, répteis)
                                //   5 = mitocondrial invertebrado
                                //   ver: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
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
