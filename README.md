# mitogenome-pipeline

Pipeline de bioinformática para montagem *de novo* e anotação automatizada de genomas mitocondriais.

Desenvolvido como Trabalho de Conclusão de Curso (TCC) — Bacharelado em Ciência da Computação, UERN.

**Autor:** Matheus Sobreira Benevides

---

## Visão Geral

Pipeline modular, containerizado (Docker) e orquestrado com Nextflow DSL2 para montagem e anotação de mitogenomas a partir de dados públicos do NCBI SRA. Cobre o ciclo completo: do dado bruto ao GenBank Flat File.

**Etapas:** Pilot QC → SRA Download → FastQC → Trim Galore → NOVOPlasty → MITOS2 → Compile Summary

---

## Resultados Obtidos

### *Anodorhynchus leari* (arara-azul-de-lear) — Estudo principal

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR28399504 (HiSeq X Ten, 2×150 bp) |
| Montagem | **Circularizada** — **16.986 bp** |
| Cobertura média | **306×** |
| Anotação | 13 CDS, 22 tRNA, 2 rRNA |
| Referência (*A. hyacinthinus*) | 16.999 bp (Δ 13 bp) |
| Tempo de execução | ~57 min (notebook i7, 16 GB RAM) |

### *Diploprion bifasciatus* — Validação

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR36182901 (NovaSeq X Plus, 2×151 bp) |
| Montagem | **Circularizada** — ~16.800 bp |

---

## Estrutura do Projeto

```
mitogenome-pipeline/
├── main.nf                  # Workflow principal (Nextflow DSL2)
├── nextflow.config          # Configuração global + perfis
├── run_pipeline.sh          # Execução resiliente (nohup)
├── modules/                 # Módulos Nextflow (um por etapa)
├── scripts/                 # Scripts auxiliares (Python/Bash)
├── docker/                  # Dockerfiles (5 imagens)
├── conf/                    # Perfis por espécie
├── data/
│   ├── seeds/               # Sementes COX1 para NOVOPlasty
│   └── databases/           # Banco MITOS2 (RefSeq89m)
└── docs/                    # Documentação e TCC
```

---

## Execução Rápida

```bash
# 1. Construir imagens Docker (uma vez)
for dir in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker build -t "mitogenome-pipeline/${dir}:1.0" "docker/${dir}"
done

# 2. Baixar banco MITOS2 (uma vez)
mkdir -p data/databases
wget https://zenodo.org/record/4284483/files/refseq89m.tar.bz2
tar -xjf refseq89m.tar.bz2 -C data/databases/ && rm refseq89m.tar.bz2

# 3. Executar
./run_pipeline.sh -profile a_leari
```

Ver o [Guia de Execução](GUIA_EXECUCAO.md) para detalhes sobre retomada de execuções, adição de novas espécies e limpeza.

---

## Entregáveis (14 categorias)

O módulo COMPILE_SUMMARY gera automaticamente em `results/<espécie>/summary/`:

- Genoma circularizado (FASTA)
- Genes codificadores em nucleotídeos e aminoácidos
- Genes ribossomais (12S + 16S) e transportadores (22 tRNAs)
- Tabelas de posição gênica, códons e anticódons
- Mapa circular do genoma (SVG/PDF)
- Anotação GFF3 completa
- Estruturas secundárias de tRNA/rRNA (SVG via ViennaRNA)
- GenBank Flat File (`.gbk`, `.tbl`, `.fsa`) pronto para submissão

---

## Ferramentas (versões pinadas)

| Ferramenta | Versão | Imagem Docker |
|---|---|---|
| SRA-Toolkit | 3.0.10 | `sra-tools:1.0` |
| FastQC | 0.12.1 | `fastqc:1.0` |
| Trim Galore + Cutadapt | 0.6.10 + 4.6 | `trim-galore:1.0` |
| NOVOPlasty | 4.3.1 | `novoplasty:1.0` |
| MITOS2 + ViennaRNA | 2.1.9 | `mitos2:1.0` |

---

## Referências

- Dierckxsens, N. et al. (2017). NOVOPlasty. *Nucleic Acids Research*.
- Di Tommaso, P. et al. (2017). Nextflow. *Nature Biotechnology*.
- Donath, A. et al. (2019). MITOS2. *Molecular Ecology Resources*.
- Lorenz, R. et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*.
- Merkel, D. (2014). Docker. *Linux Journal*.
