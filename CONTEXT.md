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

- **Espécie:** *Caenorhabditis elegans* (nematódeo modelo)
- **mtDNA:** NC_001328.1 → **13.794 bp**
- **Dataset SRA:** `SRR36152783`
  - Tamanho total: 1,18 GB comprimido / 4,39 Gb bases
  - Com limite de 500k reads: ~40 MB download / ~300 MB FASTQ
  - Cobertura estimada do mtDNA: ~544x
- **Semente:** gene cox1 (~1.530 bp) — ver `data/seeds/COMO_OBTER_SEMENTE.md`

## Objeto Principal do TCC

- **Espécie:** *Anodorhynchus leari* (arara-azul-de-lear)
- Espécie endêmica do Brasil, ameaçada de extinção
- mtDNA montado previamente em atividade acadêmica, sem publicação formal
- **Dataset SRA:** a definir (WGS paired-end Illumina)

---

## Estado Atual (30/03/2026)

### Concluído ✅
- [x] Estrutura completa do projeto Nextflow criada
- [x] 4 Dockerfiles: `sra-tools`, `fastqc`, `trim-galore`, `novoplasty`
- [x] 4 módulos Nextflow: `sra_download`, `fastqc`, `trim_galore`, `novoplasty`
- [x] `nextflow.config` (configuração global + Docker)
- [x] `conf/test.config` (perfil de teste com C. elegans)
- [x] `README.md` + `.gitignore`
- [x] Repositório GitHub criado e sincronizado
- [x] Chave SSH configurada (`~/.ssh/id_ed25519`)
- [x] Revisão de código completa (30/03) — 11 correções aplicadas:
  - `sra_download.nf`: `--maxSpotId` → `-X` (flag correta do fasterq-dump)
  - `sra_download.nf`: adicionado `errorStrategy 'retry'` + `maxRetries 2`
  - `novoplasty.nf`: `log_*.txt` → `optional: true`
  - `trim_galore.nf`: `--cores` limitado a 2 (evita estouro de CPU)
  - `fastqc.nf`: adicionado `--outdir .` (explicitude)
  - `docker/novoplasty`: wget com `-O` para nome previsível do tarball
  - `main.nf`: removido `nextflow.enable.dsl = 2` duplicado
  - `main.nf` + `nextflow.config`: SRR2081280 → SRR36152783
  - `.gitignore`: removida regra `*.py` agressiva

### Pendente ⏳
- [ ] Instalar Ubuntu no WSL2 (`wsl --install -d Ubuntu`)
- [ ] Habilitar integração Docker Desktop ↔ WSL2
- [ ] Instalar Nextflow dentro do Ubuntu
- [ ] Baixar semente cox1 de C. elegans (seguir `data/seeds/COMO_OBTER_SEMENTE.md`)
- [ ] Build das 4 imagens Docker
- [ ] Executar pipeline de teste (`nextflow run main.nf -profile test`)
- [ ] Validar output: `results/test_c_elegans/assembly/*.fasta`

---

## Pendências no Documento do TCC

- Contextualizar melhor o problema na Introdução *(nota do professor)*
- Adicionar seção `4.2.1` — SRA-Toolkit (atualmente sem seção própria)
- Definir e nomear a espécie de validação biológica de aves
- Substituir "resultados esperados" por resultados reais no Capítulo 5
- Adicionar imagens das ferramentas na Metodologia
- Padronizar para ABNT/UERN no Overleaf (LaTeX) após defesa do projeto

---

## Ambiente de Desenvolvimento

| Item | Status |
|---|---|
| SO | Windows 11 |
| Docker Desktop | 28.5.1 ✅ |
| Java | 11 ✅ |
| Git | 2.40.1 ✅ |
| WSL2 | Ativo, sem distro ⚠️ |
| Ubuntu (WSL2) | Não instalado ❌ |
| Nextflow | Não instalado ❌ |

---

## Próxima Sessão — Retomar aqui

1. Confirmar que Ubuntu/WSL2 está instalado e Docker integrado
2. Instalar Nextflow dentro do Ubuntu
3. Clonar o repositório no WSL2 e construir as imagens Docker
4. Baixar a semente e executar o pipeline de teste
