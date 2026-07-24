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

O pipeline foi aplicado a três espécies, cada uma com um papel distinto no desenho experimental.

### *Anodorhynchus leari* (arara-azul-de-lear) — Espécie-alvo do TCC

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR28399504 (HiSeq X Ten, 2×150 bp) |
| Montagem | **Circularizada** — **16.986 bp** |
| Cobertura média | **176×** |
| K-mer de circularização | 39 (1ª tentativa) |
| Anotação | 13 CDS, 22 tRNA, 2 rRNA |
| Identidade vs *A. hyacinthinus* | ~97 % nucleotídica; Δ 13 bp em tamanho total |
| Tempo de execução | ~30 min (notebook i7, 16 GB RAM) |

### *Anodorhynchus hyacinthinus* (arara-azul-grande) — Validação cruzada

Montagem *de novo* a partir de dataset metagenômico independente (SRR36400750, 2025), usando como semente o cox1 extraído da montagem de *A. leari* produzida pelo próprio pipeline. Comparação contra a referência publicada (NCBI OR209186.1, depositada por outro grupo em 2023) — dados independentes, pipelines independentes.

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR36400750 (NovaSeq X, metagenômico) |
| Montagem | **Circularizada** — **16.999 bp** |
| Cobertura média | **144×** |
| **Identidade vs OR209186.1** | **99,888 %** (19 SNPs, concentrados no D-loop) |

### *Diploprion bifasciatus* (peixe-sabão barrado) — Espécie-controle

| Métrica | Valor |
|---|---|
| Dataset SRA | SRR36182901 (NovaSeq X Plus, 2×151 bp) |
| Montagem | **Circularizada** — ~16.800 bp |
| K-mer de circularização | 33 (após fallback automático de k=39) |

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
# 1. As 5 imagens Docker já estão publicadas no GitHub Container Registry
#    com tag :1.0 (imutável); o Nextflow puxa automaticamente.
#    Para pré-baixar:
for tool in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker pull "ghcr.io/matheus-sobreira/mitogenome-pipeline/${tool}:1.0"
done

# 2. Baixar banco MITOS2 (uma vez)
mkdir -p data/databases
wget https://zenodo.org/record/4284483/files/refseq89m.tar.bz2
tar -xjf refseq89m.tar.bz2 -C data/databases/ && rm refseq89m.tar.bz2

# 3. Executar (escolha o perfil)
./run_pipeline.sh -profile a_leari            # espécie-alvo
./run_pipeline.sh -profile a_hyacinthinus     # validação cruzada
./run_pipeline.sh -profile test               # espécie-controle (D. bifasciatus)
```

Ver o [Guia de Execução](GUIA_EXECUCAO.md) para detalhes sobre retomada de execuções, adição de novas espécies e limpeza.

---

## Entregáveis

O módulo COMPILE_SUMMARY gera automaticamente em `results/<espécie>/summary/deliverables/`:

- Genoma circularizado (FASTA)
- Genes codificadores em nucleotídeos e aminoácidos
- Genes ribossomais (12S + 16S) e transportadores (22 tRNAs)
- Tabelas de posição gênica, códons start/stop e anticódons
- Mapa circular do genoma (SVG + PDF, via Biopython) e mapa linear (PNG, do MITOS2)
- Anotação GFF3 completa e ordem gênica linearizada
- Estruturas secundárias de tRNA/rRNA em SVG (via ViennaRNA/RNAplot)
- Plots de qualidade da anotação MITOS2
- GenBank Flat File (`.gbk`, `.tbl`, `.fsa`) pronto para submissão ao NCBI

Para a lista completa com nomes de arquivo, ver a [Seção 6 do Guia de Execução](GUIA_EXECUCAO.md#6-saídas).

---

## Ferramentas (versões fixadas)

| Ferramenta | Versão | Imagem Docker |
|---|---|---|
| SRA-Toolkit | 3.0.10 | `sra-tools:1.0` |
| FastQC | 0.12.1 | `fastqc:1.0` |
| Trim Galore + Cutadapt | 0.6.10 + 4.6 | `trim-galore:1.0` |
| NOVOPlasty | 4.3.1 | `novoplasty:1.0` |
| MITOS2 + ViennaRNA | 2.1.9 | `mitos2:1.0` |

---

## Licença e Uso

> ⚠️ **USO ESTRITAMENTE ACADÊMICO E NÃO COMERCIAL** — *strictly for academic, non-commercial use.*

Este pipeline orquestra ferramentas de terceiros, cada uma com sua própria licença. O montador **NOVOPlasty** é licenciado apenas para **uso não comercial** e não permite distribuição de trabalhos derivados; como o pipeline e as imagens Docker o incluem, **o conjunto destina-se exclusivamente a fins acadêmicos e de pesquisa não comerciais**. Para uso comercial, obtenha uma licença diretamente com o autor do NOVOPlasty.

Os scripts próprios (Nextflow, Dockerfiles e Python) são disponibilizados para uso acadêmico, **sem garantia**. As imagens em `ghcr.io/matheus-sobreira/...` são fornecidas "como estão", sem compromisso de atualização.

| Ferramenta | Licença |
|---|---|
| SRA-Toolkit | Domínio público (NCBI) |
| Cutadapt | MIT |
| Nextflow / Docker Engine | Apache 2.0 |
| Biopython | Biopython License (tipo BSD/MIT) |
| FastQC / Trim Galore / MITOS2 | GPL v3 |
| ViennaRNA | Licença própria (crédito obrigatório; sem redistribuição paga) |
| **NOVOPlasty** | **Licença própria — não comercial** |

Detalhes completos em [`LICENSE`](LICENSE).

---

## Referências

- Dierckxsens, N. et al. (2017). NOVOPlasty. *Nucleic Acids Research*.
- Di Tommaso, P. et al. (2017). Nextflow. *Nature Biotechnology*.
- Donath, A. et al. (2019). MITOS2. *Molecular Ecology Resources*.
- Lorenz, R. et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*.
- Merkel, D. (2014). Docker. *Linux Journal*.
