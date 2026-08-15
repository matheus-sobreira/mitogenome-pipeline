# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## O que é este repositório

Pipeline de bioinformática (Nextflow DSL2 + Docker) para montagem *de novo* e anotação automatizada de genomas mitocondriais a partir de dados públicos do NCBI SRA. É o produto do TCC (Bacharelado em Ciência da Computação, UERN) de Matheus Sobreira Benevides, aplicado à arara-azul-de-lear (*Anodorhynchus leari*).

Fluxo: `SRA_PILOT_SAMPLE → PILOT_QC → SRA_DOWNLOAD → FASTQC → TRIM_GALORE → NOVOPLASTY → ASSEMBLY_VALIDATION → [se ALERTA: RESAMPLE_POOL → TRIM_GALORE_RETRY → NOVOPLASTY_RETRY → ASSEMBLY_VALIDATION_RETRY] → MITOS2 → COMPILE_SUMMARY` (a anotação recebe a montagem **arbitrada**).

## Comandos

Não há build, lint ou suíte de testes automatizada neste repositório — é um pipeline científico avaliado pela execução ponta a ponta contra genomas de referência conhecidos, não por testes unitários.

```bash
# Executar (background resiliente via nohup, log em logs/<perfil>_<timestamp>.log)
./run_pipeline.sh -profile a_leari            # espécie-alvo do TCC
./run_pipeline.sh -profile a_hyacinthinus     # validação cruzada
./run_pipeline.sh -profile test               # espécie-controle (D. bifasciatus), anotação MITOS2 desativada por padrão

# Execução direta (sem nohup, útil para depurar em primeiro plano)
nextflow run main.nf -profile a_leari

# Retomar execução interrompida (usa cache do work/ daquele run)
nextflow run main.nf -profile a_leari -resume -w work/a_leari/2026-04-12_14h30

# Forçar um número de reads (pula a ANÁLISE do piloto; a AMOSTRA piloto ainda é
# coletada, porque a validação pós-montagem a reutiliza — DEC-20)
nextflow run main.nf -profile a_leari --sra_max_reads 20000000

# Acompanhar execução em andamento
tail -f logs/a_leari_*.log
kill -0 $(cat logs/.last_pid) && echo rodando || echo terminou

# Pré-baixar as imagens Docker (opcional; Nextflow puxa automaticamente)
for tool in sra-tools pilot-qc fastqc trim-galore novoplasty mitos2; do
    docker pull "ghcr.io/matheus-sobreira/mitogenome-pipeline/${tool}:1.0"
done

# Rebuildar uma imagem localmente após alterar um Dockerfile
# (precisa manter o mesmo namespace/tag para o Nextflow encontrá-la)
docker build -t "ghcr.io/matheus-sobreira/mitogenome-pipeline/<tool>:1.0" "docker/<tool>"

# Limpeza
rm -rf work/<especie>/        # cache de trabalho de uma espécie
rm -rf .nextflow/ .nextflow.log*
```

Para rodar apenas um script Python de pós-processamento isoladamente (fora do Nextflow, útil ao iterar em `scripts/`), invoque-o diretamente com `python3 scripts/<nome>.py --help` para ver os argumentos — cada um espera os mesmos caminhos que `modules/compile_summary.nf` passa (GFF do MITOS2, FASTA da montagem, organismo, `transl-table`).

## Arquitetura

**`main.nf`** é o orquestrador único. Valida `--sra_accession` e `--seed`, decide se roda o Pilot QC (etapa 0, opcional) com base em `params.sra_max_reads`, e condiciona as etapas MITOS2 + COMPILE_SUMMARY à presença de `params.mitos2_db` — sem esse parâmetro, o pipeline para na montagem (NOVOPlasty) e pula a anotação.

**`modules/*.nf`** — um processo Nextflow por etapa, cada um fixado (`container =`) a uma das 6 imagens Docker em `nextflow.config`. Pontos não óbvios:
- A etapa 0 são **dois** processos: `sra_pilot.nf` (`SRA_PILOT_SAMPLE`, imagem `sra-tools`) baixa a amostra, e `pilot_qc.nf` (`PILOT_QC`, imagem `pilot-qc`) a analisa. Estão separados porque a análise precisa de bwa/samtools e porque assim o `-resume` permite iterar na análise sem re-baixar a amostra.
- `sra_download.nf` usa `fasterq-dump` para o dataset completo (`sra_max_reads` ausente) e `scripts/sra_sample.sh` sobre o `.sra` já baixado pelo `prefetch` quando há limite — não são intercambiáveis.
- **Amostragem**: piloto e download usam `scripts/sra_sample.sh` (modos `head`/`stratified`/`dense`), mas com propósitos opostos — **o piloto é `stratified`** (a medição da fração precisa ser representativa) e **o download é `head`/contíguo** (DEC-23). Estratificar o download interage com a janela de intake do NOVOPlasty e desestabiliza a resolução de repeats em tandem (medido: 8M estratificado → VNTR superexpandido 18.058 bp; 8M contíguo → 16.986 bp correto). `stratified` no download permanece para comparação A/B — não o remova.
- `novoplasty.nf` filtra a saída do NOVOPlasty por `fasta.name.startsWith('Circularized_assembly')` em `main.nf`; montagens não circularizadas não seguem para MITOS2/COMPILE_SUMMARY (ver `workflow.onComplete` em `main.nf`, que reporta status `INCOMPLETO` nesse caso).
- `assembly_validation.nf` (imagem `pilot-qc`) mapeia a **amostra piloto** de volta contra a montagem circularizada e sinaliza regiões com razão de profundidade fora de `[validation_ratio_low, validation_ratio_high]`×mediana: razão baixa = repeat superexpandido na montagem, alta = repeat colapsado. Existe porque **circular ≠ correto**: uma montagem de *A. leari* circularizou com um VNTR de ~298 bp da região controle superexpandido em ~1 kb, dentro do `genome_range` (ver DEC-20 no vault). Por isso `SRA_PILOT_SAMPLE` roda mesmo com `--sra_max_reads` forçado; desativa-se com `--assembly_validation false`. Alerta exige ≥2 janelas consecutivas (isoladas são viés de cobertura, viram observação).
- **Rearbitragem (DEC-22)**: com `assembly_retry = true` (padrão), um ALERTA dispara re-sorteio do pool a partir do `.sra` (que o `SRA_DOWNLOAD` passa a preservar como output nesse modo) com semente `pilot_seed + retry_seed_offset`, re-montagem e re-validação — via aliases `*_RETRY` com `publishDir` próprios definidos no `nextflow.config` (sem os overrides, sobrescreveriam os artefatos da tentativa 0). O MITOS2 anota a montagem **arbitrada** (`final_assembly_ch`), nunca uma reprovada que tinha substituta. O Pilot QC também avisa ANTES da montagem quando a referência representa um tandem com unidade ≥ limite de resolução medido (insert real dos pares do piloto — DEC-21).
- `mitos2.nf` roda `runmitos.py` e depois gera SVGs de estrutura secundária de tRNA/rRNA chamando `RNAplot` diretamente no processo (parseia os arquivos `.nc` do MITFI manualmente).

**`conf/<especie>.config`** — um perfil por espécie (`sra_accession`, `seed`, `genome_range`, `novoplasty_kmers`, `genetic_code`, `organism`, etc.), registrado em `nextflow.config` → `profiles`. Adicionar uma espécie nova = criar um `.config` aqui + registrar o perfil em `nextflow.config`, sem tocar em `main.nf` ou nos módulos (ver seção 7 do [GUIA_EXECUCAO.md](GUIA_EXECUCAO.md)).

**`workDir` dinâmico**: `nextflow.config` computa `work/<especie>/<timestamp>/` a partir de `params.organism` no momento da chamada (bloco Groovy no final do arquivo) — cada execução fica isolada por espécie e por run, o que é o que torna `-resume -w work/<especie>/<timestamp>` necessário para retomar um run específico.

**`scripts/`** — pós-processamento Python chamado pelos módulos Nextflow (`compile_summary.nf` invoca os três abaixo em sequência), não uma biblioteca reutilizável isolada:
- `gff2genbank.py` — GFF3 do MITOS2 + FASTA → feature table GenBank (`.tbl`/`.fsa`) para submissão via `table2asn`. Contém correções específicas de domínio sobre a anotação bruta do MITOS2 (documentadas no docstring do arquivo): unificação do frameshift de ND3 (Mindell et al. 1998), descarte de fragmentos espúrios de baixo score, extensão 5'/3' do CDS até o start/stop em fase, e D-loop cobrindo toda a maior região não codificante. Essas regras existem porque o MITOS2 as reporta de forma inconsistente entre montagens — qualquer alteração aqui deve preservar esse comportamento para as três espécies já validadas (ver métricas no [README](README.md)).
- `generate_genbank.py` — gera o GenBank Flat File (`.gbk`) e o mapa circular (SVG/PDF via Biopython `GenomeDiagram`) a partir da mesma saída do MITOS2.
- `compile_summary.py` — compila os artefatos das etapas anteriores em `results/<especie>/summary/deliverables/` (14 categorias — ver seção 6 do [GUIA_EXECUCAO.md](GUIA_EXECUCAO.md)).
- `pilot_qc.sh` — lógica de análise do Pilot QC chamada por `modules/pilot_qc.nf`. A fração mitocondrial é **medida** por `bwa mem` + `samtools depth`, não estimada por k-mers: `fração = (profundidade_média × genome_avg) / (read_len × reads_totais)`. Usa profundidade e não contagem de reads mapeados porque numa referência curta (seed COX1 de ~1,5 kb) as bordas dominam a contagem; a média de profundidade apara `read_len` bp de cada extremidade. Sem `--pilot_reference`, o `genome_avg` na fórmula é o que escala a medição da seed para o mitogenoma inteiro — premissa de cobertura uniforme, declarada no relatório. Não há piso na fração (removido: mascarava o sinal de dataset metagenômico); os limites de número de reads existem e são reportados quando aplicados. Ver `Projeto/Decisoes/DEC-04` a `DEC-06` e `DEC-10` no vault Obsidian.
- `sra_sample.sh` — amostragem estratificada, chamada tanto por `modules/sra_pilot.nf` quanto por `modules/sra_download.nf`. Usa PRNG Park-Miller próprio em awk, não `rand()`, porque o produto precisa caber em 2^53 para ser exato e reprodutível entre mawk e gawk.

**`docker/<tool>/Dockerfile`** — 6 imagens com versões fixadas (SRA-Toolkit 3.0.10, BWA 0.7.18 + Samtools 1.19.2 + Seqtk 1.4, FastQC 0.12.1, Trim Galore 0.6.10 + Cutadapt 4.6, NOVOPlasty 4.3.1, MITOS2 2.1.9 + ViennaRNA), publicadas em `ghcr.io/matheus-sobreira/mitogenome-pipeline/<tool>:1.0`. O Nextflow puxa essas imagens pelo nome/tag definidos em `nextflow.config` → `process.withName`; um build local precisa ser retaggeado para o mesmo namespace para ser usado.

## Licenciamento

O **NOVOPlasty** (usado em `docker/novoplasty/`) é licenciado apenas para uso não comercial e não permite distribuição de trabalhos derivados. Como consequência, todo o pipeline e as imagens Docker publicadas se destinam exclusivamente a uso acadêmico/pesquisa não comercial (ver [LICENSE](LICENSE) e a tabela de licenças por ferramenta no [README](README.md)). Isso é relevante para qualquer decisão sobre redistribuição, empacotamento ou uso comercial de partes do pipeline.

## Documentação do TCC (fora da arquitetura do pipeline)

Além do código do pipeline, o repositório também abriga a monografia do TCC associada (`latex/`, gitignored — só existe localmente) e documentação de apoio em `docs/`. É uma frente de trabalho separada da engenharia do pipeline:

- **`docs/CONTEXTO_TCC.md`** (gitignored) é o documento de continuidade **autoritativo** para qualquer tarefa relacionada ao TCC (redação, banca, submissão ao GenBank, roadmap de publicação/mestrado) — leia-o antes de atuar nessa frente; ele mesmo declara precedência sobre memórias/anotações mais antigas.
- O TCC já foi **defendido e aprovado**; o trabalho atual é polimento final + preparação de artigo(s)/mestrado. O §9 do CONTEXTO_TCC.md é o núcleo estratégico (comentários integrais da banca, roadmap por eixo, estratégia de publicação).
- **Fluxo de edição**: o usuário edita e compila no **Overleaf**; a cópia local em `latex/` serve para diagnóstico/verificação — ao editar localmente, sempre avisar quais arquivos precisam ser replicados manualmente lá.
- Compilação local de teste: `cd latex && pdflatex -interaction=nonstopmode -halt-on-error main.tex` (2× ao mexer em paginação/sumário — sumário e referências cruzadas só estabilizam na segunda passada).

**Pendências do §12 do CONTEXTO_TCC.md — conferidas contra o repositório (12/08/2026):**
- ~~Gerar entregáveis GenBank (`.tbl`/`.fsa`/`.gbk`) de *A. hyacinthinus*~~ — **já concluído**, existem em `results/a_hyacinthinus/summary/deliverables/genbank_submission/` (gerados 24/07/2026, mesmo dia da última atualização do documento). Falta apenas a submissão em si via `table2asn`/BankIt, para as duas espécies.
- Slides de defesa, listados como "status desconhecido" — a pasta `Apresentação/` contém `roteiro-fala.md`/`.pdf`, `roteiro-slides.md` e `roteiro-completo.md` já finalizados; tudo indica que está concluído, não pendente.
- Envio do PDF final ao orientador — depende de sincronização Overleaf ↔ local e de revalidar a ficha catalográfica (§5 do CONTEXTO_TCC.md) caso a paginação tenha mudado; não verificável a partir do estado do repositório sozinho.
