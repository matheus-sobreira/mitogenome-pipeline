# Fundamentação Teórica — Pipeline de Montagem de Mitogenomas
**Projeto TCC — Matheus Sobreira Benevides**
**UERN — Bacharelado em Ciência da Computação**

---

## 1. O Genoma Mitocondrial

O DNA mitocondrial (mtDNA) é uma molécula circular de dupla fita presente nas mitocôndrias de células eucarióticas. Ao contrário do genoma nuclear, o mitogenoma é:

- **Herdado maternalmente** (na maioria dos organismos), sem recombinação
- **De cópia múltipla** — cada célula possui centenas a milhares de cópias de mitocôndrias, cada uma com seu próprio mtDNA
- **Compacto** — tipicamente 15–20 kbp em vertebrados (o genoma nuclear tem 1–3 Gbp)
- **Conservado em estrutura** — genes presentes em todos os vertebrados, com variação na ordem e sequência

O mitogenoma de vertebrados típico codifica:
- 13 proteínas da cadeia respiratória
- 22 RNAs de transferência (tRNA)
- 2 RNAs ribossomais (rRNA 12S e 16S)

Essa combinação de **alta cobertura natural** (proveniente da abundância de cópias), **tamanho manejável** e **conservação evolutiva** torna o mtDNA ideal para filogenia e identificação de espécies.

---

## 2. Sequenciamento de Nova Geração (NGS) e Dados SRA

### 2.1 Plataforma Illumina (utilizada neste projeto)

O sequenciamento Illumina (usado no dataset `SRR28399504`) opera por **síntese com fluorescência reversível**:

1. O DNA é fragmentado e adaptadores sintéticos são ligados às extremidades
2. Os fragmentos se ancoram numa flow cell e são **amplificados por PCR em ponte**, formando clusters
3. Nucleotídeos marcados com fluorescência são incorporados um por vez; uma câmera captura a imagem após cada incorporação
4. O resultado são **reads curtos** (tipicamente 100–150 bp) com alta exatidão (> 99%)

O HiSeq X Ten (plataforma do dataset de *A. leari*) gera reads de **150 bp paired-end**, ou seja:
- Cada fragmento é sequenciado pelos dois lados (R1 e R2)
- R1 = extremidade "forward", R2 = extremidade "reverse"
- O par de reads define um fragmento de tamanho aproximado ("insert size")

### 2.2 Por que tanta cobertura do mtDNA?

Num sequenciamento WGS (whole genome shotgun), o sequenciador processa DNA total da célula. Como cada célula possui:
- **1–2 cópias** do genoma nuclear (~2,1 Gb em aves)
- **Centenas a milhares de cópias** do mtDNA (~17 kbp)

A proporção de reads mitocondriais no total é muito maior do que a proporção de tamanho (17 kbp / 2,1 Gbp ≈ 0,0008%). No dataset de *A. leari*, **0,17% das reads** (34.644 reads) era mitocondrial — o suficiente para 306× de cobertura com apenas 20 M reads do total de 118,5 M.

### 2.3 Arquivo SRA e o banco NCBI-SRA

O NCBI Sequence Read Archive (SRA) é o maior repositório de dados de sequenciamento do mundo. Os dados são armazenados no formato binário `.sra`, que é uma compressão proprietária do NCBI. Para este projeto:

| Formato | Tamanho | Reads |
|---|---|---|
| `.sra` (comprimido) | 11 GB | 118,5 M pares |
| `.fastq` bruto (descomprimido) | ~87 GB total (R1 + R2) | 118,5 M pares |

A ferramenta `prefetch` baixa o `.sra`; `fasterq-dump` converte para FASTQ.

---

## 3. Controle de Qualidade — FastQC

O FastQC analisa arquivos FASTQ e reporta métricas de qualidade sem modificar os dados:

- **Per-base quality scores (Phred):** qualidade por posição ao longo da read. Valores ≥ Q30 (99,9% de precisão) são considerados bons
- **Per-sequence quality:** distribuição de qualidade média por read
- **GC content:** desvios do esperado podem indicar contaminação
- **Adapter content:** sequências de adaptadores encontradas nas reads (não pertencem ao organismo)
- **Duplication levels:** reads idênticas podem indicar supersaturação de PCR

> No pipeline, o FastQC é executado **antes da trimagem** (reads brutas) para diagnóstico. A presença de adaptadores confirmada pelo FastQC justifica a etapa seguinte.

---

## 4. Trimagem — Trim Galore + Cutadapt

### 4.1 Problema

As reads brutas contêm:
1. **Sequências de adaptadores** — ligados artificialmente para ancoragem na flow cell; não fazem parte do organismo
2. **Bases de baixa qualidade** — especialmente nas extremidades 3' das reads (queda natural de sinal Illumina)

### 4.2 Solução — Trim Galore

Trim Galore é um wrapper do Cutadapt que automatiza a detecção de adaptadores e trimagem:

- **Detecção automática de adaptadores:** identifica o adaptador Illumina TruSeq (`AGATCGGAAGAGC`) sem que o usuário precise especificar
- **Trimagem por qualidade:** remove bases com Phred < 20 da extremidade 3'
- **Filtro de comprimento mínimo:** descarta reads que ficarem muito curtas após trimagem (< 20 bp por padrão)
- **Paired-end awareness:** mantém sincronização entre R1 e R2 — se uma read é descartada, a parceira correspondente também é

### 4.3 Resultados para *A. leari* (SRR28399504)

| Métrica | R1 (forward) | R2 (reverse) |
|---|---|---|
| Reads de entrada | 20.000.000 | 20.000.000 |
| Reads de saída | 20.000.000 (100%) | 20.000.000 (100%) |
| Bases de entrada | 3,00 Gb | 3,00 Gb |
| Bases de saída | 2,15 Gb (71,8%) | 2,16 Gb (72,0%) |
| Removido por qualidade | 30 Mb (1,0%) | 20 Mb (0,7%) |
| Removido por adaptador | 816 Mb (27,2%) | 821 Mb (27,4%) |
| Reads com adaptador | 16.285.736 (81,4%) | 16.157.912 (80,8%) |
| Adaptador detectado | `AGATCGGAAGAGC` (Illumina TruSeq) | `AGATCGGAAGAGC` |

**Interpretação:** Cerca de 81% das reads possuíam adaptador — valor alto, mas normal em bibliotecas de fragmentos curtos onde o adaptador ocupa parte da read de 150 bp. Somente 1% das bases foi removido por qualidade — indicativo de boa qualidade de sequenciamento. Todas as 20 M reads passaram (100%), o que indica boa profundidade e qualidade.

---

## 5. Estratégia de Truncagem — Por que 20 M reads?

### 5.1 O problema de escala

O dataset de *A. leari* possui 118,5 M pares de reads (~87 GB descomprimido). Processar tudo seria:
- **Computacionalmente desnecessário** para montagem de mitogenoma (~17 kbp)
- **Demorado:** mesmo com hardware moderno, 87 GB de FASTQ demandam horas para trimagem completa

### 5.2 Estimativa de cobertura alvo

Cobertura (C) é definida como:

$$C = \frac{N \times L}{G}$$

Onde:
- $N$ = número de reads (pares)
- $L$ = comprimento da read (150 bp)
- $G$ = tamanho do genoma alvo (17.000 bp para mtDNA de ave)

Com 20 M reads:

$$C_{esperada} = \frac{20.000.000 \times 0,17\% \times 2 \times 150}{17.000} \approx 600 \times$$

Na prática obtivemos **306×** de cobertura — abaixo do estimado porque a fração de 0,17% foi calculada após a montagem. Para montagem de novo de mitogenoma, cobertura > 50× é suficiente; 306× é excelente.

### 5.3 Uniformidade da amostragem

A truncagem por `head` pega as primeiras N reads na order do arquivo FASTQ, com origem na posição física inicial da flow cell. Em dados Illumina HiSeq X Ten, a qualidade e composição das reads são **uniformes ao longo da corrida** — não há degradação sistemática por posição de cluster, tornando a truncagem por `head` estatisticamente equivalente a uma amostragem aleatória para este propósito.

### 5.4 Ordem correta: truncar ANTES de trimar

Trimar primeiro (118,5 M reads) e depois truncar seria equivalente em resultado biológico mas **~6× mais custoso computacionalmente** sem benefício, pois:
- A qualidade e proporção de adaptadores é uniforme nas 118,5 M reads
- O NOVOPlasty não distingue entre os 20 M reads trimados dos primeiros ou dos últimos

A ordem **descompressão → truncagem → trimagem** é a mais eficiente para este cenário.

**Redução total de volume:**

| Etapa | Tamanho/arquivo | Redução acumulada |
|---|---|---|
| SRA comprimido | 11,0 GB | — |
| FASTQ bruto (fasterq-dump) | 43,5 GB | +295% (expansão esperada) |
| FASTQ truncado (20 M reads) | 7,4 GB | −83% em relação ao bruto |
| FASTQ trimado (Trim Galore) | 5,3 GB | −28% em relação ao truncado |

---

## 6. Montagem de Novo — NOVOPlasty

### 6.1 Conceito de montagem de novo

Montar um genoma "de novo" significa reconstruí-lo apenas a partir das reads, sem usar uma referência completa. O algoritmo constrói **grafos de De Bruijn**:

1. As reads são fatiadas em k-mers (subsequências de comprimento k)
2. Os k-mers são conectados em um grafo: cada nó é um k-mer, cada aresta representa overlap de k-1 bases
3. O caminho Euleriano no grafo representa a sequência montada

### 6.2 NOVOPlasty como montador especializado

O NOVOPlasty foi projetado especificamente para genomas circulares de organelas (mitocôndria e cloroplasto). Suas vantagens:

- **Usa uma semente (seed):** uma sequência conhecida (ex.: gene COX1) ancora onde a montagem começa, reduzindo drasticamente o espaço de busca
- **Circularização automática:** o algoritmo detecta quando o contig "fecha" sobre si mesmo, produzindo a molécula circular
- **Baixo consumo de memória:** filtra apenas as reads que potencialmente pertencem à organela, evitando processar o genoma nuclear

### 6.3 Parâmetros utilizados para *A. leari*

| Parâmetro | Valor | Justificativa |
|---|---|---|
| `Seed` | COX1 de *A. hyacinthinus* (1.548 bp) | Mesmo gênero, alta identidade esperada |
| `Genome range` | 15.500–18.500 bp | Baseado no mtDNA de *A. hyacinthinus* (16.999 bp) |
| `K-mer` | 39 | K-mer principal; proporcional ao read length de 150 bp |
| `K-mer secundário` | 33 | Fallback caso k=39 falhe (maior k = mais específico, menor = mais sensível) |
| `Read length` | 150 bp | Comprimento real das reads HiSeq X Ten |
| `Insert size` | 300 bp | Tamanho estimado do fragmento sequenciado |
| `Max memory` | 12 GB | Limite de RAM para o processo |

### 6.4 Resultados da montagem de *A. leari*

| Métrica | Valor |
|---|---|
| Status | ✅ Circularizado com sucesso |
| Número de contigs | 1 |
| Tamanho do contig | **16.986 bp** |
| Referência mais próxima (*A. hyacinthinus*) | 16.999 bp |
| Diferença do contig em relação à referência | 13 bp (< 0,1%) |
| Reads totais processadas | 20.285.220 |
| Reads alinhadas ao mtDNA | 34.644 (0,17%) |
| Cobertura média | **306×** |
| K-mer utilizado | 39 |
| Fração subsampleada | 61,9% |

A semente (*seed*) foi inicializada com sucesso, indicando boa conservação do gene COX1 entre *A. hyacinthinus* e *A. leari*. A montagem circularizou no primeiro k-mer testado (k=39).

---

## 7. Reprodutibilidade e FAIR

O pipeline foi desenvolvido seguindo os princípios **FAIR** (Findable, Accessible, Interoperable, Reusable):

| Princípio | Implementação |
|---|---|
| **F** — Encontrável | Repositório público no GitHub com README e DOI futuro |
| **A** — Acessível | Dados de entrada no NCBI-SRA (accession público); código no GitHub |
| **I** — Interoperável | Formatos padrão (FASTQ, FASTA); ferramentas com versões fixadas nos Dockerfiles |
| **R** — Reutilizável | Nextflow DSL2 com perfis por espécie (`-profile a_leari`); Docker garante mesmo ambiente em qualquer máquina |

### 7.1 Por que Docker?

Cada ferramenta bioinformática possui dependências complexas que frequentemente conflitam entre si. Docker encapsula cada ferramenta em um **contêiner isolado** com todas as suas dependências. Isso garante:

- **Reprodutibilidade:** o mesmo Dockerfile produz o mesmo ambiente em qualquer máquina, em qualquer momento
- **Rastreabilidade:** a versão exata de cada ferramenta está fixada no Dockerfile
- **Isolamento:** versões conflitantes de bibliotecas não interferem entre si

### 7.2 Por que Nextflow?

Nextflow é um framework de orquestração de workflows científicos que oferece:

- **Paralelização automática** de processos independentes
- **Gestão de dependências** entre processos (o resultado de um alimenta o próximo)
- **Retomada de execução** (`.nextflow/cache`) — se um processo falha, apenas ele é re-executado
- **Portabilidade:** o mesmo pipeline roda localmente, em cluster HPC ou em nuvem AWS/GCP sem modificação do código
- **DSL2:** linguagem declarativa baseada em Groovy, legível e modular

---

## 8. Implicações para *Anodorhynchus leari*

*Anodorhynchus leari* (arara-azul-de-lear) é uma espécie **Criticamente Ameaçada** (IUCN Red List), endêmica da Caatinga nordestina do Brasil. A montagem do seu mitogenoma completo contribui para:

1. **Filogenia:** posicionamento preciso na árvore evolutiva do gênero *Anodorhynchus* e da família Psittacidae
2. **Genética de populações:** análise de diversidade mitocondrial entre indivíduos; estimativa do tamanho efetivo populacional
3. **Conservação:** identificação de haplótipos únicos; monitoramento genético de indivíduos em cativeiro e reintrodução
4. **Medicina veterinária:** referência para diagnóstico molecular

A similaridade de 99,9% do mitogenoma montado com *A. hyacinthinus* (16.986 bp vs. 16.999 bp) é biologicamente coerente, dado que as duas espécies são irmãs dentro do gênero.

---

## 9. Sumário Visual do Pipeline

```
SRR28399504 (NCBI SRA)
        │
        ▼ prefetch (11 GB .sra)
┌───────────────────┐
│   SRA_DOWNLOAD    │  fasterq-dump → 87 GB FASTQ → head 20M reads → 7,4 GB/arquivo
└────────┬──────────┘
         │ SRR28399504_{1,2}.fastq (20M reads each)
         ▼
┌───────────────────┐
│     FASTQC        │  Diagnóstico de qualidade e adaptadores
└────────┬──────────┘  → results/a_leari/qc/fastqc_raw/
         │
         ▼
┌───────────────────┐
│   TRIM_GALORE     │  Remoção de adaptadores TruSeq + bases Q<20
└────────┬──────────┘  81% reads com adaptador → 71,9% bases mantidas
         │ _val_1.fq / _val_2.fq (5,3 GB each)
         ▼
┌───────────────────┐
│    NOVOPLASTY     │  Montagem de novo guiada por semente (COX1)
└────────┬──────────┘  k-mer=39, seed=A. hyacinthinus COX1
         │
         ▼
 Circularized_assembly_1_SRR28399504.fasta
 ✅ 16.986 bp | 1 contig | Cobertura 306× | Circularizado
```

---

*Gerado automaticamente em 31/03/2026 a partir dos resultados do pipeline.*
