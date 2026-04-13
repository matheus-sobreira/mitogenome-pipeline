# Contexto do Projeto — TCC Matheus Sobreira Benevides

> Arquivo de referência para retomar o trabalho entre sessões.
> **NÃO versionar** (está no .gitignore).

---

## Projeto

**Título:** Construção de um Pipeline Reprodutível para Montagem e Anotação Automatizada de Mitogenomas: Aplicação à Arara-azul-de-lear (*Anodorhynchus leari*)
**Autor:** Matheus Sobreira Benevides
**Curso:** Bacharelado em Ciência da Computação — UERN
**Orientador:** Prof. Wilfredo | **Prof. TCC:** Prof. Carlos
**Repositório:** https://github.com/matheus-sobreira/mitogenome-pipeline
**Branch de backup Fase 1:** `fase1-montagem`

---

## Stack Tecnológica

| Componente | Tecnologia | Versão |
|---|---|---|
| Orquestração | Nextflow DSL2 | 25.10.4 |
| Containerização | Docker | — |
| Controle de versão | Git + GitHub | — |
| Princípios | FAIR + reprodutibilidade computacional | — |

---

## Pipeline — 7 Etapas

| Etapa | Ferramenta | Módulo | Status |
|---|---|---|---|
| 1 | SRA Pilot QC | `modules/sra_download.nf` + `scripts/pilot_qc.sh` | ✅ |
| 2 | SRA-Toolkit (`fasterq-dump`) | `modules/sra_download.nf` | ✅ |
| 3 | FastQC | `modules/fastqc.nf` | ✅ |
| 4 | Trim Galore + Cutadapt | `modules/trim_galore.nf` | ✅ |
| 5 | NOVOPlasty | `modules/novoplasty.nf` | ✅ |
| 6 | MITOS2 | `modules/mitos2.nf` | ✅ |
| 7 | Compile Summary | `modules/compile_summary.nf` + scripts Python | ✅ |

### Scripts Python
- `scripts/pilot_qc.sh` — análise piloto de qualidade (500K reads)
- `scripts/compile_summary.py` — compilação de 14 entregáveis
- `scripts/gff2genbank.py` — conversão GFF3 → feature table GenBank (.tbl + .fsa)
- `scripts/generate_genbank.py` — geração de GenBank Flat File (.gbk) + mapa circular (SVG/PDF)

### Dockerfiles (5 imagens, versões pinadas)
- `docker/sra-tools/Dockerfile` — SRA-Toolkit 3.0.10 (ubuntu:22.04)
- `docker/fastqc/Dockerfile` — FastQC 0.12.1 (ubuntu:22.04)
- `docker/trim-galore/Dockerfile` — Trim Galore 0.6.10 + Cutadapt 4.6 (ubuntu:22.04)
- `docker/novoplasty/Dockerfile` — NOVOPlasty 4.3.1 (ubuntu:22.04)
- `docker/mitos2/Dockerfile` — MITOS2 2.1.9 + ViennaRNA (miniconda3:24.1.2-0)

---

## Espécies

### Validação — *Diploprion bifasciatus* (peixe-sabão barrado)
- **mtDNA publicado:** PZ143763.1 → 16.805 bp
- **Dataset SRA:** `SRR36182901` (NovaSeq X Plus, 2×151 bp, 12 M pares, ~1,2 GB)
- **Semente:** gene cox1 de *D. bifasciatus* (PZ143763.1, 1.560 bp)
- **Configuração:** `conf/test.config` + `-profile test`
- **Resultado:** montagem circularizada (k-mer 33)

### Objeto Principal — Arara-azul-de-lear (*Anodorhynchus leari*)
- Espécie endêmica do Brasil (Raso da Catarina, Bahia), EN/IUCN
- **Dataset SRA:** `SRR28399504` (HiSeq X Ten, 118,5 M spots, ~11 GB)
- **Semente:** cox1 de *A. hyacinthinus* (NC_082165.1, 1.548 bp)
- **Configuração:** `conf/a_leari.config` + `-profile a_leari`
- **Resultado montagem:** 16.986 bp circularizado, 306×, k-mer 39
- **Resultado anotação (MITOS2):** 13 CDS + 22 tRNAs + 2 rRNAs + D-loop
- **Frameshift ND3:** detectado (nad3_0 + nad3_1), consistente com Psittaciformes

---

## Estado Atual (13/04/2026)

### Concluído ✅

**Infraestrutura:**
- [x] Projeto Nextflow DSL2 completo (main.nf + 5 módulos + 4 scripts)
- [x] 5 Dockerfiles com versões pinadas
- [x] Repositório GitHub sincronizado
- [x] Branch `fase1-montagem` como backup da Fase 1
- [x] workDir dinâmico por espécie/timestamp (`work/<espécie>/<timestamp>/`)
- [x] Script `run_pipeline.sh` resiliente a desconexão (nohup + logs)
- [x] `GUIA_EXECUCAO.md` com instruções completas

**Fase 1 — Montagem:**
- [x] Pipeline D. bifasciatus — montagem circularizada (k-mer 33)
- [x] Pipeline A. leari — 16.986 bp, 306×, circularizado
- [x] Pilot QC implementado (500K reads → análise → max_reads automático, cap 25M)

**Fase 2 — Anotação:**
- [x] MITOS2 integrado — 13 CDS + 22 tRNAs + 2 rRNAs + D-loop
- [x] Estruturas secundárias em SVG (RNAplot/ViennaRNA)
- [x] GenBank Flat File (.gbk) gerado automaticamente
- [x] Feature table (.tbl + .fsa) para submissão ao GenBank
- [x] Mapa circular (SVG + PDF) via Biopython GenomeDiagram

**Compilação de Entregáveis (summary/):**
- [x] 14 categorias de entregáveis gerados automaticamente
- [x] Nomes com organismo (`A_leari.gbk` em vez de `Contig1.gbk`)

**Documento TCC:**
- [x] `docs/TCC_ATUALIZADO.md` — documento completo (~1.150 linhas, 8 capítulos)
- [x] Título atualizado conforme sugestão do professor
- [x] Cap 1 reestruturado (6 seções: Contextualização → Problema → Objetivos → Justificativa → Metodologia → Estrutura)
- [x] Cap 4 reordenado (fluxo de execução antes das ferramentas)
- [x] Cap 5.7 — Análise de Desempenho com métricas reais
- [x] Cap 6 — novo capítulo dedicado à arara-azul-de-lear (7 seções)
- [x] Frameshift ND3 aprofundado com 5 referências adicionais
- [x] ~46 termos biológicos definidos para público de CC
- [x] 33 sugestões de figura/tabela posicionadas no texto
- [x] 40 referências bibliográficas (ABNT), zero órfãs
- [x] Otimização implementada: cap Pilot QC de 50M → 25M reads
- [x] Arquivos obsoletos removidos: `atualizacoes_tcc.md`, `fundamentacao_teorica.md`

### Pendente ⏳
- [ ] Executar pipeline limpo e validar COMPILE_SUMMARY end-to-end
- [ ] Criar figuras/tabelas referenciadas nas 33 sugestões do TCC
- [ ] Padronizar TCC para ABNT/UERN no Overleaf (LaTeX)
- [ ] Montar slides de defesa
- [ ] Submeter mitogenoma de A. leari ao GenBank

---

## Ambiente de Desenvolvimento

| Item | Status |
|---|---|
| SO | Ubuntu (nativo / Linux) |
| Docker | Instalado ✅ |
| Java | 11+ ✅ |
| Nextflow | 25.10.4 ✅ |
| Git | Instalado ✅ |
| Biopython | Instalado ✅ |
| matplotlib | Instalado ✅ |
| Banco MITOS2 | `data/databases/refseq89m` (gitignored) ✅ |

---

## Próxima Sessão — Retomar aqui

### Prioridade 1 — Execução limpa e validação
- Limpar work/ e results/ (se não foi feito)
- Executar: `./run_pipeline.sh -profile a_leari`
- Validar que COMPILE_SUMMARY gera summary/ correto (sem arquivos MITOS2 misturados)
- Confirmar 14 entregáveis completos

### Prioridade 2 — Figuras do TCC
- Gerar as figuras referenciadas nas 33 sugestões do documento
- Pelo menos as essenciais: fluxograma do pipeline, mapa circular OGDRAW, diagrama de fluxo de dados

### Prioridade 3 — Formatação e Defesa
- Migrar TCC para LaTeX/Overleaf (template UERN ABNT)
- Montar slides de defesa
- Submeter mitogenoma ao GenBank
