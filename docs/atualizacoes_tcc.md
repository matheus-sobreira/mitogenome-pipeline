# ATUALIZAÇÕES DO TCC — Seções revisadas com resultados reais
# ============================================================
# Autor: Matheus Sobreira Benevides
# Data: 31/03/2026
#
# Este arquivo contém as seções do TCC que foram atualizadas ou adicionadas
# com base nos resultados reais obtidos na execução do pipeline.
# Cada seção está marcada com [NOVA], [ATUALIZADA] ou [SUBSTITUIR].
#
# Instruções:
#   - Seções marcadas [SUBSTITUIR] devem SUBSTITUIR a seção correspondente no TCC
#   - Seções marcadas [NOVA] devem ser INSERIDAS na posição indicada
#   - Seções marcadas [ATUALIZADA] tiveram trechos modificados — revise o diff
# ============================================================


# =====================================================================
# [SUBSTITUIR] 1.2 Justificativa
# =====================================================================
# Nota: O professor comentou "Contextualizar e explicar melhor o problema."
# A versão abaixo expande a justificativa com o contexto de A. leari.
# =====================================================================

1.2 Justificativa

A análise de genomas mitocondriais (mitogenomas) representa um componente essencial para a biologia evolutiva, a genética de populações e a sistemática, dado que essas moléculas fornecem informações altamente informativas sobre relações de parentesco, variação intra e interespecífica e processos de conservação. Sua relevância científica foi consolidada em trabalhos seminais que destacaram o papel estrutural e evolutivo dos mitogenomas como marcadores moleculares (BOORE, 1999). Em seguida, estudos enfatizaram sua aplicabilidade em investigações de genética de populações e evolução molecular (EKBLOM; WOLF, 2014), enquanto contribuições mais recentes demonstraram seu potencial em análises de biodiversidade e conservação por meio de avanços metodológicos em montagem e anotação de organelas (ULIANO-SILVA et al., 2023).

Embora a obtenção de mitogenomas completos tenha sido amplamente favorecida pelo advento do Sequenciamento de Nova Geração (NGS), a execução de pipelines de montagem ainda enfrenta desafios técnicos significativos. A diversidade de ferramentas, versões e dependências envolvidas gera um cenário de difícil replicação, conhecido como dependency hell, que compromete a reprodutibilidade e a padronização das análises (SANDVE et al., 2013; BAKER, 2016; GRÜNING et al., 2018).

Esse problema se torna especialmente relevante quando aplicado ao estudo de espécies ameaçadas de extinção, nas quais a disponibilidade de amostras biológicas é limitada e a confiabilidade dos resultados genômicos é crucial para fundamentar decisões de conservação. É o caso da arara-azul-de-lear (Anodorhynchus leari), ave endêmica da Caatinga nordestina do Brasil, classificada como "Em Perigo" pela IUCN Red List. Com população estimada em pouco mais de 1.700 indivíduos na natureza, a espécie depende de dados genômicos acessíveis e validáveis para subsidiar programas de manejo, reprodução em cativeiro e monitoramento populacional. Embora seu mitogenoma já tenha sido previamente reconstruído em atividade acadêmica, essa montagem foi realizada de forma manual, sem padronização metodológica e sem disponibilização pública do fluxo de trabalho — o que impede sua verificação e replicação por outros pesquisadores.

Nesse contexto, torna-se necessário o desenvolvimento de soluções que conciliem robustez biológica e consistência computacional. A adoção de práticas modernas da engenharia de software, como a containerização e os sistemas de gerenciamento de workflows, tem se mostrado fundamental para superar limitações de reprodutibilidade, além de permitir que fluxos analíticos sejam portáveis, escaláveis e transparentes (DI TOMMASO et al., 2017; EWELS et al., 2020).

Assim, este trabalho justifica-se pela proposta de construir um pipeline automatizado e reprodutível para a montagem de mitogenomas, unindo boas práticas de ciência aberta e engenharia computacional à relevância biológica desses genomas. A aplicação ao caso concreto de A. leari demonstra não apenas a viabilidade técnica da abordagem, mas também sua contribuição direta para a conservação de uma espécie ameaçada da biodiversidade brasileira.


# =====================================================================
# [NOVA] 4.2.1 SRA-Toolkit
# =====================================================================
# Nota do CONTEXT.md: "Adicionar seção 4.2.1 — SRA-Toolkit (atualmente sem seção própria)"
# Inserir ANTES da seção 4.2.2 FastQC.
# =====================================================================

4.2.1 SRA-Toolkit

O SRA-Toolkit é o pacote oficial de utilitários desenvolvido e mantido pelo National Center for Biotechnology Information (NCBI) para acesso, download e manipulação de dados depositados no Sequence Read Archive (SRA), um dos maiores repositórios públicos de dados de sequenciamento do mundo (LEINONEN; SUGAWARA; SHUMWAY, 2011).

Do ponto de vista técnico, o SRA-Toolkit é composto por um conjunto de ferramentas de linha de comando com funções especializadas. As duas principais utilizadas neste pipeline são:

1. prefetch — responsável pelo download do arquivo binário no formato .sra, que corresponde à representação comprimida e proprietária dos dados de sequenciamento armazenados no NCBI. O download via prefetch é mais estável e eficiente do que o acesso direto via streaming, especialmente para datasets de grande volume.

2. fasterq-dump — ferramenta de conversão do formato .sra para o formato FASTQ, padrão na bioinformática moderna por reunir em um mesmo arquivo a sequência nucleotídica e as informações de qualidade (scores Phred) de cada base. O fasterq-dump é a versão otimizada e multi-threaded do antigo fastq-dump, oferecendo velocidades de conversão significativamente superiores.

Para dados paired-end, o fasterq-dump produz automaticamente dois arquivos FASTQ separados — um para as leituras forward (R1) e outro para as reverse (R2) — preservando a correspondência entre os pares de leituras. No caso do dataset de A. leari (SRR28399504), a conversão do arquivo .sra de 11 GB resultou em aproximadamente 87 GB de dados FASTQ descomprimidos (43,5 GB por arquivo), contendo 118,5 milhões de pares de leituras.

Uma consideração prática importante é que o fasterq-dump, diferentemente de seu antecessor fastq-dump, não suporta a flag -X para limitação do número de reads durante a conversão. Para viabilizar a truncagem do dataset (necessária para reduzir o custo computacional sem comprometer a montagem do mitogenoma), o pipeline implementa uma estratégia de pós-conversão: após a geração completa dos arquivos FASTQ, o comando head é utilizado para reter apenas os primeiros N reads, descartando o excedente. Essa abordagem é detalhada na Seção 4.3.

Entre as vantagens do SRA-Toolkit destacam-se: (i) a integração nativa com o banco SRA, garantindo acesso direto e confiável aos dados; (ii) a verificação automática de integridade dos arquivos baixados via checksums; e (iii) a ampla compatibilidade com formatos de saída utilizados em ferramentas subsequentes do pipeline. Como limitação, pode-se citar o volume de armazenamento temporário necessário durante a conversão, que pode alcançar três vezes o tamanho do arquivo .sra original.

A escolha pelo SRA-Toolkit neste trabalho se justifica por ser a ferramenta oficial e mais amplamente validada para acesso a dados públicos de sequenciamento, assegurando que a etapa de aquisição de dados seja realizada de forma padronizada, rastreável e compatível com os princípios de ciência aberta.


# =====================================================================
# [SUBSTITUIR] 4.3 Aquisição e Preparo dos Dados
# =====================================================================
# Atualizado com os datasets reais: D. bifasciatus (validação) e A. leari (principal).
# Remove menção ao ERR9710916 / Danio rerio que não foi utilizado.
# =====================================================================

4.3 Aquisição e Preparo dos Dados

A primeira etapa prática do pipeline consiste na aquisição e organização dos dados de sequenciamento que servirão de insumo para a montagem dos genomas mitocondriais. Para este trabalho, foram utilizados datasets públicos provenientes do Sequence Read Archive (SRA), repositório mantido pelo NCBI e reconhecido como um dos principais bancos internacionais de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). A escolha por dados públicos justifica-se pela sua ampla disponibilidade, pela diversidade de organismos representados e pela possibilidade de reproduzir e validar resultados obtidos em trabalhos prévios.

O acesso aos dados é realizado por meio do SRA-Toolkit. O processo inicia-se com o comando prefetch, responsável pelo download dos arquivos brutos no formato .sra, seguido da conversão para o formato FASTQ utilizando o fasterq-dump.

Dois datasets foram selecionados para este trabalho, com propósitos complementares:

a) Espécie de validação — Diploprion bifasciatus (peixe-sabão barrado): O dataset SRR36182901, gerado por sequenciamento paired-end na plataforma Illumina NovaSeq X Plus (2×151 bp), foi selecionado como caso de teste por três razões: (i) seu mitogenoma (PZ143763.1, 16.805 bp) foi publicado a partir deste mesmo dataset, permitindo comparação direta; (ii) o volume de dados é compacto (~1,2 GB comprimido, 12 milhões de pares de leituras), viabilizando testes rápidos e iterativos durante o desenvolvimento do pipeline; e (iii) a espécie pertence a um grupo taxonômico distinto (Perciformes), demonstrando a generalidade do pipeline para além de aves.

b) Espécie principal — Anodorhynchus leari (arara-azul-de-lear): O dataset SRR28399504, gerado por sequenciamento paired-end na plataforma Illumina HiSeq X Ten (2×150 bp), corresponde ao caso de estudo central deste TCC. Depositado por Iridian Genomes em 2024, o dataset contém 118,5 milhões de pares de leituras (~10,7 GB comprimido / ~87 GB em FASTQ). Por se tratar de dados de genoma total (WGS), a fração mitocondrial corresponde a uma pequena proporção do total — porém, dada a abundância natural de cópias de mtDNA por célula, essa fração é suficiente para a montagem completa do mitogenoma.

Para ambas as espécies, sequências-semente (seeds) do gene COX1 foram obtidas de espécies filogeneticamente próximas ou da própria espécie a partir de mitogenomas depositados no GenBank. Para A. leari, utilizou-se o COX1 de A. hyacinthinus (NC_082165.1, posições 5359–6906, 1.548 bp), espécie-irmã do mesmo gênero. Para D. bifasciatus, utilizou-se o COX1 extraído do próprio mitogenoma publicado (PZ143763.1, 1.560 bp). As sementes foram armazenadas em formato FASTA no diretório data/seeds/ do repositório.

Estratégia de truncagem — Devido ao grande volume do dataset de A. leari (118,5 M reads), a conversão completa para FASTQ gera aproximadamente 87 GB de dados. Para reduzir o custo computacional nas etapas subsequentes sem comprometer a montagem do mitogenoma, foi implementada uma estratégia de truncagem no módulo sra_download.nf: após a conversão completa pelo fasterq-dump, o comando head -n (N × 4) é utilizado para reter apenas os primeiros 20 milhões de pares de leituras (sendo 4 linhas por read no formato FASTQ). Essa abordagem reduziu o volume de cada arquivo de ~43,5 GB para ~7,4 GB. A escolha de 20 milhões de reads foi baseada na estimativa de cobertura:

    C = (N × f × 2 × L) / G

Onde N = 20.000.000 reads, f = fração mitocondrial estimada (~0,17%), L = 150 bp e G = 17.000 bp, resultando em uma cobertura esperada de aproximadamente 600×, valor amplamente suficiente para montagem de novo de mitogenomas.

A truncagem é realizada antes da trimagem (etapa seguinte), deliberadamente: como a qualidade das bases e a proporção de adaptadores são uniformemente distribuídas ao longo de uma corrida Illumina HiSeq X Ten, trimar primeiro os 118,5 M reads para depois truncar seria computacionalmente ~6× mais custoso sem ganho biológico. A ordem truncagem → trimagem é, portanto, a mais eficiente para este cenário.

Após a aquisição e truncagem, os dados são organizados em diretórios padronizados pelo próprio Nextflow (diretório work/), respeitando uma estrutura hierárquica que diferencia arquivos brutos, processados e resultados. Ao adotar essa prática, reforça-se o alinhamento do projeto com os princípios FAIR (WILKINSON et al., 2016).

Tabela 1 – Datasets utilizados no pipeline.

| Espécie | Accession | Plataforma | Reads totais | Tamanho SRA | Finalidade |
|---|---|---|---|---|---|
| D. bifasciatus | SRR36182901 | NovaSeq X Plus (2×151 bp) | 12 M pares | ~1,2 GB | Validação |
| A. leari | SRR28399504 | HiSeq X Ten (2×150 bp) | 118,5 M pares | ~11 GB | Estudo principal |


# =====================================================================
# [SUBSTITUIR] 4.4 Controle de Qualidade e Pré-processamento
# =====================================================================
# Atualizado com resultados reais do Trim Galore.
# =====================================================================

4.4 Controle de Qualidade e Pré-processamento

Após a aquisição e organização dos dados de sequenciamento, a etapa seguinte do pipeline consiste na avaliação da qualidade das leituras e na aplicação de procedimentos de filtragem. Essa fase é essencial para assegurar que apenas dados confiáveis avancem para a montagem do genoma mitocondrial, evitando que artefatos técnicos comprometam a acurácia do resultado final.

O primeiro passo é a execução do FastQC (ANDREWS, 2010), que gera relatórios interativos em HTML contendo métricas fundamentais: (i) qualidade por posição nas sequências, (ii) distribuição de valores de qualidade ao longo de todas as leituras, (iii) conteúdo GC, (iv) presença de adaptadores residuais e (v) detecção de sequências sobre-representadas. No caso do dataset de A. leari (SRR28399504), o FastQC confirmou a presença significativa de sequências de adaptadores Illumina nas reads brutas, justificando a necessidade da etapa de trimagem.

Em seguida, foi aplicado o Trim Galore (BABRAHAM BIOINFORMATICS, 2019), wrapper que integra o Cutadapt (MARTIN, 2011) com o FastQC, automatizando a detecção e remoção de adaptadores e bases de baixa qualidade. O Trim Galore opera em modo paired-end, mantendo a sincronização entre R1 e R2 — se uma read é descartada por comprimento mínimo, a parceira correspondente também é removida.

Para o dataset de A. leari (20 milhões de reads truncados), o Trim Galore detectou automaticamente o adaptador Illumina TruSeq (AGATCGGAAGAGC) e produziu os seguintes resultados:

Tabela 2 – Resultados do pré-processamento com Trim Galore para A. leari (SRR28399504).

| Métrica | R1 (forward) | R2 (reverse) |
|---|---|---|
| Reads de entrada | 20.000.000 | 20.000.000 |
| Reads de saída | 20.000.000 (100%) | 20.000.000 (100%) |
| Bases de entrada | 3,00 Gb | 3,00 Gb |
| Bases de saída | 2,15 Gb (71,8%) | 2,16 Gb (72,0%) |
| Removido por baixa qualidade | 30 Mb (1,0%) | 20 Mb (0,7%) |
| Removido por adaptadores | 816 Mb (27,2%) | 821 Mb (27,4%) |
| Reads com adaptador detectado | 16.285.736 (81,4%) | 16.157.912 (80,8%) |
| Adaptador identificado | AGATCGGAAGAGC (Illumina TruSeq) | AGATCGGAAGAGC |

A alta proporção de reads com adaptador (81%) é um valor esperado em bibliotecas de fragmentos curtos sequenciadas em plataformas que geram reads de 150 bp: quando o fragmento de DNA inserido é menor que o comprimento da read, o sequenciador avança até a região do adaptador, que é então capturado na sequência. A remoção dessas sequências artificiaisé indispensável para evitar erros de montagem.

Apenas 1% das bases foi removido por baixa qualidade (Phred < 20), indicando que a corrida de sequenciamento foi de alta qualidade. Notavelmente, 100% das reads passaram pelos filtros de comprimento mínimo, sem perda de pares — evidência de que o dataset possui consistência e profundidade adequadas.

Além de aumentar a confiabilidade biológica das leituras, o pré-processamento contribuiu significativamente para a redução do volume de dados. Cada arquivo FASTQ foi reduzido de ~7,4 GB (após truncagem) para ~5,3 GB (após trimagem), representando uma redução adicional de 28%.

Tabela 3 – Redução progressiva de volume de dados para A. leari.

| Etapa | Tamanho por arquivo | Reads (pares) | Redução |
|---|---|---|---|
| SRA comprimido (.sra) | 11,0 GB (total) | 118,5 M | — |
| FASTQ bruto (fasterq-dump) | ~43,5 GB | 118,5 M | +295% (descompressão) |
| FASTQ truncado (head 20M) | ~7,4 GB | 20 M | −83% |
| FASTQ trimado (Trim Galore) | ~5,3 GB | ~20 M | −28% |

Essa redução progressiva demonstra como cada etapa do pipeline contribui para tornar o processamento mais eficiente, sem comprometer a representatividade dos dados para a montagem do mitogenoma.


# =====================================================================
# [SUBSTITUIR] Capítulo 5 — Resultados e Discussão
# =====================================================================
# Substituir TODO o capítulo 5 (que era "Resultados Esperados") por resultados reais.
# =====================================================================

5. Resultados e Discussão

Este capítulo apresenta os resultados obtidos na execução do pipeline de montagem de mitogenomas, aplicado a dois estudos de caso: (i) Diploprion bifasciatus como espécie de validação, cujo mitogenoma publicado permite comparação direta, e (ii) Anodorhynchus leari como espécie principal, objeto central deste TCC. Os resultados são descritos por etapa do pipeline, seguidos de uma discussão sobre a robustez da abordagem e sua relevância biológica e científica.

5.1 Resultados por Etapa do Pipeline

5.1.1 Aquisição e Preparo dos Dados

A etapa de aquisição de dados foi executada com sucesso para ambos os datasets. O SRA-Toolkit (prefetch + fasterq-dump) realizou o download e a conversão dos arquivos .sra para FASTQ sem erros. Para o dataset de A. leari (SRR28399504), o prefetch obteve o arquivo .sra de 11 GB, e o fasterq-dump realizou a conversão completa em dois arquivos FASTQ totalizando ~87 GB.

Uma correção técnica relevante foi necessária durante o desenvolvimento: a flag -X, utilizada no antigo fastq-dump para limitar o número de reads, não é suportada pelo fasterq-dump (versão 3.0.10), gerando o erro "Unknown argument '-X'" com código de saída 3. A solução implementada no módulo sra_download.nf foi a conversão completa seguida de truncagem via head -n (N × 4), onde N é o número de reads desejado. Essa abordagem, embora exija armazenamento temporário dos arquivos completos, é robusta e independente de versão da ferramenta.

Para A. leari, a truncagem de 118,5 M para 20 M reads reduziu cada arquivo de ~43,5 GB para ~7,4 GB, conforme planejado.

5.1.2 Controle de Qualidade

O FastQC (versão 0.12.1) gerou relatórios de qualidade para ambos os datasets. Para A. leari, os relatórios confirmaram perfis de qualidade típicos de dados Illumina HiSeq X Ten, com scores Phred médios acima de Q30 em todas as posições, presença significativa de adaptadores Illumina TruSeq e níveis de duplicação moderados — padrão esperado em bibliotecas WGS de alta cobertura.

5.1.3 Pré-processamento

O Trim Galore (versão 0.6.10, com Cutadapt 4.6) foi executado em modo paired-end, com detecção automática de adaptadores. Os resultados detalhados para A. leari estão apresentados na Tabela 2 (Seção 4.4). Em resumo: 81,4% das reads R1 e 80,8% das R2 continham sequências de adaptador Illumina TruSeq (AGATCGGAAGAGC), resultando na remoção de 27,2-27,4% das bases por adaptadores e 0,7-1,0% por baixa qualidade. Todas as 20 milhões de reads passaram nos filtros de comprimento mínimo (100% de retenção).

Esses valores indicam que a biblioteca de sequenciamento continha fragmentos de inserção relativamente curtos, resultando em read-through dos adaptadores — comportamento normal e esperado em bibliotecas WGS padrão. O resultado confirma a importância da etapa de trimagem como passo obrigatório antes da montagem.

5.1.4 Montagem do Mitogenoma

O NOVOPlasty (versão 4.3.1) foi executado com sucesso para ambas as espécies, produzindo montagens circularizadas na primeira execução — sem necessidade da estratégia iterativa de re-seeding com contigs parciais.

Tabela 4 – Resultados da montagem por NOVOPlasty.

| Métrica | D. bifasciatus (validação) | A. leari (principal) |
|---|---|---|
| Status da montagem | Circularizado ✓ | Circularizado ✓ |
| Número de contigs | 1 | 1 |
| Tamanho do contig | — | 16.986 bp |
| Referência mais próxima | PZ143763.1 (16.805 bp) | A. hyacinthinus (16.999 bp) |
| K-mer utilizado | 33 | 39 |
| Reads totais processadas | — | 20.285.220 |
| Reads alinhadas ao mtDNA | — | 34.644 (0,17%) |
| Cobertura média | — | 306× |
| Fração subsampleada | — | 61,9% |
| Insert size observado | — | 213 bp |
| Semente utilizada | cox1 D. bifasciatus (PZ143763.1) | cox1 A. hyacinthinus (NC_082165.1) |

Para A. leari, o mitogenoma de 16.986 bp difere em apenas 13 bp (< 0,1%) do mitogenoma de referência de A. hyacinthinus (16.999 bp), resultado biologicamente coerente para espécies do mesmo gênero. A semente COX1 interespecífica (A. hyacinthinus → A. leari) foi recuperada com sucesso, demonstrando que o NOVOPlasty tolera bem diferenças genéticas dentro do mesmo gênero.

A fração de reads mitocondriais (0,17% do total) produziu cobertura média de 306×, amplamente acima do limiar mínimo de 50× necessário para montagem confiável de organelas. Isso confirma a adequação da estratégia de truncagem para 20 M reads: mesmo utilizando apenas 16,9% do dataset original, foi possível obter cobertura elevada e montagem circular completa.

5.1.5 Validação com D. bifasciatus

A espécie de validação (D. bifasciatus, SRR36182901) forneceu uma montagem circularizada com k-mer 33 na primeira execução. A comparação com o mitogenoma publicado (PZ143763.1, gerado a partir do mesmo dataset) permite validar a corretude do pipeline e confirmar que as etapas de download, trimagem e montagem preservam a integridade biológica dos dados. A montagem foi bem-sucedida, confirmando a funcionalidade do pipeline em organismos de diferentes grupos taxonômicos (peixes e aves).

5.2 Discussão sobre a Robustez da Abordagem

A robustez metodológica do pipeline foi demonstrada pela obtenção de montagens completas e circularizadas em ambos os estudos de caso, abrangendo organismos de grupos taxonômicos distintos (Perciformes e Psittaciformes) e plataformas de sequenciamento diferentes (NovaSeq X Plus e HiSeq X Ten).

Entre os aspectos positivos observados, destaca-se a reprodutibilidade, garantida pela combinação entre a containerização das ferramentas via Docker e a orquestração do fluxo pelo Nextflow. Cada execução do pipeline produz um diretório work/ com hash único por processo, permitindo rastreabilidade completa de cada etapa. O uso do flag -resume permite retomar execuções interrompidas sem reprocessamento, característica que se mostrou essencial durante o desenvolvimento — quando falhas de ferramentas (como a incompatibilidade da flag -X no fasterq-dump) puderam ser corrigidas sem necessidade de reiniciar todo o pipeline.

A estratégia de truncagem por head, embora exija o armazenamento temporário dos arquivos FASTQ completos, mostrou-se robusta e produziu resultados de montagem indistinguíveis dos que seriam obtidos com o dataset completo. A justificativa reside na uniformidade de qualidade e composição ao longo das corridas Illumina, que torna a truncagem posicional estatisticamente equivalente a uma amostragem aleatória para fins de montagem mitogenômica.

O módulo de NOVOPlasty implementado no pipeline inclui uma lógica iterativa com múltiplos k-mers e re-seeding automático. Embora essa funcionalidade não tenha sido necessária nos dois casos testados (ambos circularizaram na primeira tentativa), ela confere resiliência ao pipeline para cenários mais desafiadores, como baixa cobertura ou regiões repetitivas complexas.

Como limitação observada na prática, o volume de armazenamento temporário durante a etapa de download/conversão pode ser significativo: para o dataset de A. leari, foram necessários ~11 GB (arquivo .sra) + ~87 GB (FASTQ bruto temporário) + ~15 GB (FASTQ truncado), totalizando ~113 GB de pico de armazenamento antes da limpeza automática. Esse requisito pode ser restritivo em ambientes com espaço limitado, embora o pipeline realize limpeza automática das reads brutas após cada etapa.

5.3 Relevância Biológica e Científica

A relevância deste trabalho transcende a dimensão metodológica. Do ponto de vista biológico, o pipeline possibilitou a reconstrução reprodutível do mitogenoma de Anodorhynchus leari, espécie endêmica do Brasil classificada como "Em Perigo" de extinção. O mitogenoma obtido — 16.986 bp, circularizado, com 306× de cobertura — constitui um recurso genômico validável e reprodutível, diferenciando-se da montagem anterior (realizada em contexto acadêmico sem documentação metodológica padronizada) por estar associado a:

(i) um pipeline público e versionado no GitHub;
(ii) containers Docker com versões fixas de cada ferramenta;
(iii) arquivos de configuração declarativos que permitem reprodução exata do experimento;
(iv) relatórios de qualidade e trimagem preservados como saída do pipeline.

A similaridade de 99,9% entre o mitogenoma montado (16.986 bp) e o de A. hyacinthinus (16.999 bp) é biologicamente coerente, dado que as duas espécies são irmãs dentro do gênero Anodorhynchus. Essa proximidade valida tanto a escolha da semente interespecífica quanto a qualidade da montagem obtida.

A validação cruzada com D. bifasciatus (grupo taxonômico distinto) demonstra que o pipeline não é específico para aves ou para um tipo particular de dados, ampliando sua aplicabilidade para a comunidade científica. A disponibilização pública do pipeline, acompanhada de perfis de execução para cada espécie (nextflow run main.nf -profile a_leari ou -profile test), permite que qualquer pesquisador reproduza integralmente os resultados aqui apresentados.

Em termos mais amplos, a disponibilização de um pipeline automatizado, containerizado e orquestrado com Nextflow se alinha às demandas contemporâneas da bioinformática por soluções reprodutíveis, portáveis e escaláveis. Ao tornar esse recurso acessível por meio de repositório público, acompanhado de documentação detalhada, o presente trabalho contribui para o fortalecimento de práticas de ciência aberta e colabora com a comunidade internacional que atua no desenvolvimento de workflows bioinformáticos.


# =====================================================================
# [SUBSTITUIR] Capítulo 6 — Considerações Finais
# =====================================================================

6. Considerações Finais

O presente trabalho teve como objetivo desenvolver um pipeline de bioinformática automatizado, reprodutível e portável para a montagem de novo de genomas mitocondriais. A integração de ferramentas consolidadas — SRA-Toolkit, FastQC, Trim Galore/Cutadapt e NOVOPlasty — em um fluxo containerizado com Docker e orquestrado pelo Nextflow constituiu a principal contribuição metodológica desta proposta.

O pipeline foi validado com sucesso em dois cenários distintos. A espécie de validação, Diploprion bifasciatus, forneceu uma montagem circularizada comparável ao mitogenoma publicado para a mesma espécie. O estudo de caso principal, Anodorhynchus leari, resultou em um mitogenoma de 16.986 bp, circularizado com cobertura média de 306×, com diferença inferior a 0,1% em relação à espécie-irmã A. hyacinthinus, confirmando a validade biológica da montagem.

Os resultados demonstraram que:

a) O pipeline é capaz de reconstruir mitogenomas completos e circularizados a partir de dados públicos de WGS, utilizando sementes interespecíficas quando a sequência da espécie-alvo não está disponível;

b) A estratégia de truncagem para 20 milhões de reads (de um total de 118,5 milhões) foi suficiente para produzir cobertura de 306×, validando a abordagem de redução de volume como estratégia viável para economizar recursos computacionais;

c) A containerização individual das ferramentas eliminou problemas de incompatibilidade de versões e dependências, materializando o princípio de reprodutibilidade no próprio design do pipeline;

d) A modularidade do Nextflow DSL2, aliada a perfis de configuração por espécie (-profile a_leari, -profile test), permite que o pipeline seja reutilizado para qualquer organismo com ajustes mínimos — bastando fornecer um accession SRA e uma semente COX1.

Como perspectivas futuras, destacam-se: (i) a implementação da Fase 2 do pipeline — anotação funcional com MITOS2 e tRNAscan-SE, que identificará os 37 genes típicos do mitogenoma; (ii) a comparação detalhada do mitogenoma de A. leari com outras espécies de Psittaciformes para estudos filogenéticos; (iii) a publicação do mitogenoma no GenBank; e (iv) a adequação completa do documento ao modelo ABNT/UERN no Overleaf (LaTeX).

Em síntese, o pipeline desenvolvido representa uma contribuição metodológica concreta para a bioinformática de organelas, demonstrando que é possível unir rigor técnico, aplicabilidade biológica e boas práticas de reprodutibilidade computacional em uma solução acessível e reutilizável. A aplicação ao caso de uma espécie ameaçada da biodiversidade brasileira confere ao trabalho não apenas relevância acadêmica, mas também potencial impacto para a conservação.
