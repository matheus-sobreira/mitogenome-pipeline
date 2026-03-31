# Contexto do Projeto — TCC Matheus Sobreira Benevides

> Arquivo de referência para retomar o trabalho entre sessões.
> **NÃO versionar** (está no .gitignore).

---

## Projeto

**Título:** Desenvolvimento e Implementação de um Pipeline de Bioinformática para Montagem e Análise de Mitogenomas
**Autor:** Matheus Sobreira Benevides
**Curso:** Bacharelado em Ciência da Computação — UERN
**Orientador:** Prof. Wilfredo | **Prof. TCC:** Prof. Carlos
**Repositório:** https://github.com/matheus-sobreira/mitogenome-pipeline

---

## Stack Tecnológica

| Componente | Tecnologia |
|---|---|
| Orquestração | Nextflow DSL2 |
| Containerização | Docker |
| Controle de versão | Git + GitHub |
| Princípios | FAIR + reprodutibilidade computacional |

---

## Pipeline

### Fase 1 — Validação (implementada, aguardando execução)

| Etapa | Ferramenta | Módulo |
|---|---|---|
| 1 | SRA-Toolkit (`fasterq-dump`) | `modules/sra_download.nf` |
| 2 | FastQC | `modules/fastqc.nf` |
| 3 | Trim Galore + Cutadapt | `modules/trim_galore.nf` |
| 4 | NOVOPlasty | `modules/novoplasty.nf` |

### Fase 2 — Anotação (não iniciada)

| Etapa | Ferramenta |
|---|---|
| 5 | MITOS2 |
| 6 | tRNAscan-SE |

---

## Espécie de Validação (testes do pipeline)

- **Espécie:** *Diploprion bifasciatus* (peixe-sabão barrado)
- **mtDNA publicado:** PZ143763.1 → **16.805 bp**
- **Dataset SRA:** `SRR36182901` (NovaSeq X Plus, 2×151 bp, 12 M pares, ~1,2 GB)
- **Semente:** gene cox1 de *D. bifasciatus* (PZ143763.1, 1.560 bp)
- **Configuração:** `conf/test.config` + perfil `test` em `nextflow.config`
- **Resultado:** montagem circularizada com sucesso (k-mer 33)
  - `results/test_d_bifasciatus/assembly/Circularized_assembly_1_SRR36182901.fasta`

## Objeto Principal do TCC

- **Espécie:** *Anodorhynchus leari* (arara-azul-de-lear)
- Espécie endêmica do Brasil, ameaçada de extinção
- mtDNA montado previamente em atividade acadêmica, sem publicação formal
- **Dataset SRA:** `SRR28399504` (HiSeq X Ten, 118,5 M spots, ~11 GB comprimido)
- **Semente COX1:** *Anodorhynchus hyacinthinus* — NC_082165.1, pos. 5359–6906, 1.548 bp
  - Mesmo gênero; mtDNA de referência: 16.999 bp
- **Arquivo semente:** `data/seeds/a_hyacinthinus_cox1.fasta`
- **Configuração:** `conf/a_leari.config` + perfil `a_leari` em `nextflow.config`

---

## Estado Atual (31/03/2026)

### Concluído ✅
- [x] Estrutura completa do projeto Nextflow criada
- [x] 4 Dockerfiles: `sra-tools`, `fastqc`, `trim-galore`, `novoplasty`
- [x] 4 módulos Nextflow: `sra_download`, `fastqc`, `trim_galore`, `novoplasty`
- [x] `nextflow.config` (configuração global + Docker)
- [x] `README.md` + `.gitignore`
- [x] Repositório GitHub criado e sincronizado
- [x] Chave SSH configurada (`~/.ssh/id_ed25519`)
- [x] Revisão de código completa (30/03) — 11 correções aplicadas:
  - `sra_download.nf`: `--maxSpotId` → `-X` (flag correta do fasterq-dump)
  - `sra_download.nf`: adicionado `errorStrategy 'retry'` + `maxRetries 2`
  - `sra_download.nf`: estratégia `prefetch` + `fasterq-dump` local (mais confiável)
  - `novoplasty.nf`: `log_*.txt` → `optional: true`
  - `novoplasty.nf`: estratégia iterativa com múltiplos k-mers e re-seeding automático
  - `trim_galore.nf`: `--cores` limitado a 2 (evita estouro de CPU)
  - `trim_galore.nf`: limpeza automática de reads brutos após uso
  - `fastqc.nf`: adicionado `--outdir .` (explicitude)
  - `docker/novoplasty`: wget com `-O` para nome previsível do tarball
  - `main.nf`: removido `nextflow.enable.dsl = 2` duplicado
  - `main.nf` + `nextflow.config`: SRR2081280 → SRR36152783
  - `.gitignore`: removida regra `*.py` agressiva
- [x] Ambiente Ubuntu nativo configurado (Docker, Nextflow, Java, Git)
- [x] 4 imagens Docker construídas com sucesso
- [x] Perfil de teste reformulado: *C. elegans* → *Diploprion bifasciatus* (SRR36182901)
  - Dataset mais compacto (~1.2 GB), mtDNA publicado a partir deste dataset
  - Semente: gene cox1 de D. bifasciatus (PZ143763.1, 1.560 bp)
- [x] Sementes cox1 obtidas para 6 espécies adicionais (`data/seeds/`)
- [x] **Pipeline executado com sucesso para D. bifasciatus** — montagem circularizada:
  - `results/test_d_bifasciatus/assembly/Circularized_assembly_1_SRR36182901.fasta`
- [x] **`sra_download.nf` corrigido** — flag `-X` inválida no `fasterq-dump`:
  - Nova estratégia: conversão completa + `head -n (max_reads × 4)` para truncagem
- [x] **Semente A. leari obtida:** `data/seeds/a_hyacinthinus_cox1.fasta`
  - Fonte: NC_082165.1 (*A. hyacinthinus*), posições 5359–6906, 1.548 bp
- [x] **Configuração A. leari criada:** `conf/a_leari.config`
  - `sra_accession = 'SRR28399504'`, `sra_max_reads = 20000000`
  - `genome_range = '15500-18500'`, `novoplasty_kmers = '39,33'`
- [x] **Pipeline A. leari executado com sucesso (31/03/2026):**
  - SRA_DOWNLOAD → 11 GB .sra → 87 GB FASTQ → truncado para 20 M reads (~7,4 GB/arquivo)
  - FASTQC → relatórios em `results/a_leari/qc/fastqc_raw/`
  - TRIM_GALORE → 81,4% das reads com adaptador Illumina TruSeq detectado; 71,9% bases mantidas
  - NOVOPLASTY:
    - **Montagem circularizada com sucesso** — 1 contig
    - **Tamanho:** 16.986 bp (referência *A. hyacinthinus*: 16.999 bp — diferença de 13 bp, < 0,1%)
    - **Cobertura média:** 306×
    - **Reads alinhadas ao mitogenoma:** 34.644 / 20.285.220 (0,17%)
    - **k-mer utilizado:** 39
    - **Resultado:** `results/a_leari/assembly/Circularized_assembly_1_SRR28399504.fasta`
- [x] Gráfico de redução de volume de dados gerado: `pipeline_data_reduction.png`
- [x] Arquivo de fundamentação teórica para TCC: `docs/fundamentacao_teorica.md`

### Pendente ⏳
- [ ] **Fase 2 — Anotação funcional** (MITOS2 + tRNAscan-SE):
  - Criar Dockerfile para MITOS2 (`docker/mitos2/Dockerfile`)
  - Criar Dockerfile para tRNAscan-SE (`docker/trnascan/Dockerfile`)
  - Criar módulo Nextflow `modules/mitos2.nf`
  - Criar módulo Nextflow `modules/trnascan.nf`
  - Integrar no `main.nf` após NOVOPLASTY
  - Executar para *A. leari* — espera-se 13 proteínas, 22 tRNAs, 2 rRNAs
  - Gerar mapa circular anotado do mitogenoma
- [ ] Validar biologicamente a montagem de D. bifasciatus (BLAST contra PZ143763.1)
- [ ] Publicar mitogenoma de *A. leari* no GenBank

---

## Pendências no Documento do TCC

- ~~Contextualizar melhor o problema na Introdução~~ ✅ (Justificativa expandida em `docs/atualizacoes_tcc.md`)
- ~~Adicionar seção `4.2.1` — SRA-Toolkit~~ ✅ (adicionada em `docs/atualizacoes_tcc.md`)
- ~~Substituir "resultados esperados" por resultados reais no Capítulo 5~~ ✅ (Cap. 5 reescrito)
- [ ] Adicionar imagens das ferramentas e gráficos na Metodologia
- [ ] Usar `pipeline_data_reduction.png` e `docs/fundamentacao_teorica.md` na apresentação
- [ ] Adicionar resultados da Fase 2 (anotação) ao Cap. 5 quando disponíveis
- [ ] Padronizar para ABNT/UERN no Overleaf (LaTeX) após defesa do projeto
- **TCC versão completa atualizada:** `docs/atualizacoes_tcc.md`

---

## Ambiente de Desenvolvimento

| Item | Status |
|---|---|
| SO | Ubuntu (nativo / Linux) |
| Docker | Instalado ✅ |
| Java | 11+ ✅ |
| Git | Instalado ✅ |
| Nextflow | Instalado ✅ |

---

## Próxima Sessão — Retomar aqui

### Prioridade 1 — Fase 2: Anotação Funcional (MITOS2 + tRNAscan-SE)

**Objetivo:** Anotar o mitogenoma de *A. leari* (16.986 bp) identificando os 37 genes esperados.

**Passos:**
1. Criar `docker/mitos2/Dockerfile` — MITOS2 com dependência do banco de dados RefSeq
2. Criar `docker/trnascan/Dockerfile` — tRNAscan-SE 2.0
3. Criar `modules/mitos2.nf` — recebe FASTA da montagem, produz GFF3 + mapa circular
4. Criar `modules/trnascan.nf` — recebe FASTA, produz tabela de tRNAs
5. Integrar ambos no `main.nf` (após `NOVOPLASTY`)
6. Adicionar perfis/parâmetros no `nextflow.config` (código genético mitocondrial de vertebrados = 2)
7. Executar: `nextflow run main.nf -profile a_leari`
8. Validar: comparar genes anotados com *A. hyacinthinus* (NC_082165.1)
9. Atualizar Cap. 5 do TCC com resultados da anotação

**Resultado esperado:** 13 genes codificadores de proteínas + 22 tRNAs + 2 rRNAs + região controle (D-loop)

**Input:** `results/a_leari/assembly/Circularized_assembly_1_SRR28399504.fasta`

### Prioridade 2 — Pós-anotação
1. Validar biologicamente a montagem de *D. bifasciatus* (BLAST contra PZ143763.1)
2. Montar slides do TCC usando `docs/fundamentacao_teorica.md` e `pipeline_data_reduction.png`
3. Publicar mitogenoma no GenBank
