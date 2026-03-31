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

## Resultados Obtidos

### *Anodorhynchus leari* (arara-azul-de-lear) — Estudo principal

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR28399504 (HiSeq X Ten, 2×150 bp) |
| Reads utilizadas | 20 M pares (de 118,5 M totais) |
| Montagem | **Circularizada** ✓ |
| Tamanho | **16.986 bp** |
| Cobertura média | **306×** |
| Referência (*A. hyacinthinus*) | 16.999 bp (diferença < 0,1%) |
| Semente | COX1 *A. hyacinthinus* (NC_082165.1) |

### *Diploprion bifasciatus* (peixe-sabão barrado) — Validação

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR36182901 (NovaSeq X Plus, 2×151 bp) |
| Montagem | **Circularizada** ✓ |
| Referência publicada | PZ143763.1 (16.805 bp) |

---

## Estrutura do Projeto

```
mitogenome-pipeline/
├── main.nf                # Workflow principal (Nextflow DSL2)
├── nextflow.config        # Configuração global + Docker + perfis
├── modules/               # Módulos Nextflow (um por ferramenta)
├── docker/                # Dockerfiles de cada ferramenta
├── conf/
│   ├── test.config        # Perfil de teste — D. bifasciatus
│   └── a_leari.config     # Perfil principal — A. leari
├── data/seeds/            # Sequências semente (COX1) para o NOVOPlasty
├── docs/                  # Fundamentação teórica e TCC atualizado
└── results/               # Saídas do pipeline (montagens, QC, relatórios)
```

---

## Pré-requisitos

- [Docker](https://www.docker.com/) (Linux nativo ou Docker Desktop com WSL2)
- [Nextflow](https://www.nextflow.io/) >= 23.04
- Java >= 11

---

## Execução

### 1. Construir as imagens Docker (apenas uma vez)

```bash
cd mitogenome-pipeline
docker build -t mitogenome-pipeline/sra-tools:1.0   docker/sra-tools/
docker build -t mitogenome-pipeline/fastqc:1.0       docker/fastqc/
docker build -t mitogenome-pipeline/trim-galore:1.0  docker/trim-galore/
docker build -t mitogenome-pipeline/novoplasty:1.0   docker/novoplasty/
```

### 2. Rodar o pipeline

```bash
# Estudo principal — A. leari (arara-azul-de-lear)
nextflow run main.nf -profile a_leari

# Validação — D. bifasciatus (peixe-sabão barrado)
nextflow run main.nf -profile test

# Retomar execução interrompida
nextflow run main.nf -profile a_leari -resume
```

### 3. Resultados

```
results/a_leari/
├── qc/
│   ├── fastqc_raw/         ← Relatórios HTML do FastQC
│   └── trim_galore/        ← Relatórios de trimagem
└── assembly/
    ├── Circularized_assembly_1_SRR28399504.fasta   ← Mitogenoma montado
    ├── config.txt           ← Configuração do NOVOPlasty
    └── log_SRR28399504.txt  ← Log da montagem
```

### 4. Criar novo perfil para outra espécie

1. Criar `conf/nova_especie.config` (copiar `conf/a_leari.config` como base)
2. Alterar `sra_accession`, `seed`, `outdir`, `genome_range`
3. Adicionar no `nextflow.config`:
   ```groovy
   nova_especie { includeConfig 'conf/nova_especie.config' }
   ```
4. Rodar: `nextflow run main.nf -profile nova_especie`

---

## Parâmetros do Pipeline

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `sra_accession` | — | Acesso SRA (obrigatório) |
| `seed` | — | Arquivo FASTA com a semente COX1 (obrigatório) |
| `outdir` | `results` | Diretório de saída |
| `sra_max_reads` | `null` | Limite de read pairs (truncagem via `head`) |
| `genome_range` | `15000-18500` | Tamanho esperado do mtDNA (bp) |
| `novoplasty_kmers` | `'39,33'` | K-mers para NOVOPlasty (tentados em ordem) |
| `read_length` | `150` | Comprimento das reads |
| `insert_size` | `300` | Tamanho estimado do fragmento |

---

## Referências

- Andrews, S. (2010). FastQC.
- Dierckxsens, N. et al. (2017). NOVOPlasty. *Nucleic Acids Research*.
- Di Tommaso, P. et al. (2017). Nextflow. *Nature Biotechnology*.
- Krueger, F. (2019). Trim Galore. Babraham Bioinformatics.
- Martin, M. (2011). Cutadapt. *EMBnet.journal*.
- Merkel, D. (2014). Docker. *Linux Journal*.
