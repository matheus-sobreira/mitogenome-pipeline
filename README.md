# mitogenome-pipeline

Pipeline de bioinformática para montagem e anotação automatizada de genomas mitocondriais.

Desenvolvido como Trabalho de Conclusão de Curso (TCC) — Bacharelado em Ciência da Computação, UERN.

**Autor:** Matheus Sobreira Benevides

---

## Visão Geral

Pipeline modular, containerizado (Docker) e orquestrado com Nextflow para:

1. Download de leituras do NCBI SRA (`SRA-Toolkit`)
2. Controle de qualidade (`FastQC`)
3. Remoção de adaptadores e trimming (`Trim Galore`)
4. Montagem do mitogenoma (`NOVOPlasty`)

> **Fase 2 (em desenvolvimento):** anotação funcional com MITOS2 e tRNAscan-SE.

---

## Estrutura do Projeto

```
mitogenome-pipeline/
├── main.nf                # Workflow principal (Nextflow DSL2)
├── nextflow.config        # Configuração global + Docker
├── build_images.ps1       # Script para construir imagens Docker (Windows/PowerShell)
├── modules/               # Módulos Nextflow (um por ferramenta)
├── docker/                # Dockerfiles de cada ferramenta
├── conf/
│   └── test.config        # Perfil de teste — C. elegans (SRR36152783)
└── data/seeds/            # Sequências semente para o NOVOPlasty
```

---

## Pré-requisitos

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (com WSL2 Integration ativada)
- [Nextflow](https://www.nextflow.io/) >= 23.04  
- Java >= 11

---

## Execução Rápida (Fase de Validação)

### 1. Construir as imagens Docker (apenas uma vez)

```bash
# Linux / WSL2
cd mitogenome-pipeline
docker build -t mitogenome-pipeline/sra-tools:1.0   docker/sra-tools/
docker build -t mitogenome-pipeline/fastqc:1.0       docker/fastqc/
docker build -t mitogenome-pipeline/trim-galore:1.0  docker/trim-galore/
docker build -t mitogenome-pipeline/novoplasty:1.0   docker/novoplasty/
```

### 2. Obter a semente

Siga as instruções em [`data/seeds/COMO_OBTER_SEMENTE.md`](data/seeds/COMO_OBTER_SEMENTE.md).

### 3. Rodar o pipeline de teste

```bash
nextflow run main.nf -profile test
```

Espécie de validação: *Caenorhabditis elegans* — mtDNA completo: **13.794 bp** (NC_001328.1)  
Dataset: `SRR36152783` | Limite de teste: 500.000 read pairs (~300 MB FASTQ) | Cobertura estimada: ~544x

### 4. Resultados

```
results/test_c_elegans/
├── reads/raw/          ← FASTQs brutos
├── reads/trimmed/      ← FASTQs após trimming
├── qc/fastqc_raw/      ← Relatórios HTML do FastQC
└── assembly/           ← FASTA do mitogenoma montado + log NOVOPlasty
```

---

## Parâmetros do Pipeline

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `--sra_accession` | — | Acesso SRA (obrigatório) |
| `--seed` | — | Arquivo FASTA com a semente (obrigatório) |
| `--outdir` | `results` | Diretório de saída |
| `--sra_max_reads` | `null` | Limite de read pairs para download |
| `--genome_range` | `13000-14500` | Tamanho esperado do mtDNA (bp) |
| `--trim_quality` | `20` | Qualidade Phred mínima (Trim Galore) |
| `--min_length` | `50` | Comprimento mínimo pós-trimming |

---

## Referências

- Andrews, S. (2010). FastQC.
- Dierckxsens, N. et al. (2017). NOVOPlasty. *Nucleic Acids Research*.
- Di Tommaso, P. et al. (2017). Nextflow. *Nature Biotechnology*.
- Krueger, F. (2019). Trim Galore. Babraham Bioinformatics.
- Martin, M. (2011). Cutadapt. *EMBnet.journal*.
- Merkel, D. (2014). Docker. *Linux Journal*.
