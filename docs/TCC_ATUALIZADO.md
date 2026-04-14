# TCC — Construção de um Pipeline Reprodutível para Montagem e Anotação Automatizada de Mitogenomas: Aplicação à Arara-azul-de-lear (*Anodorhynchus leari*)

**Autor:** Matheus Sobreira Benevides
**Curso:** Bacharelado em Ciência da Computação — UERN
**Orientador:** Prof. Wilfredo
**Prof. TCC:** Prof. Carlos e Prof. Bruno

---

## Capítulo 1 — Introdução

### 1.1 Contextualização e Tema

Os genomas mitocondriais — também denominados **mitogenomas** — constituem o conjunto completo de material genético contido na mitocôndria, **organela** (estrutura subcelular com função especializada) responsável pela produção de energia nas **células eucarióticas** (células que possuem núcleo delimitado por membrana — categoria que inclui animais, plantas e fungos, em oposição às bactérias). Essas moléculas de DNA ocupam uma posição central na biologia evolutiva, genética de populações e sistemática, uma vez que oferecem informações valiosas sobre parentesco, variação intra e interespecífica e processos de conservação. O DNA mitocondrial, por suas características estruturais e funcionais — como a **herança predominantemente materna** (transmissão exclusiva pela linhagem feminina, sem contribuição paterna), a ausência de **recombinação** (processo em que **cromossomos**, as estruturas lineares de DNA no núcleo celular, trocam segmentos durante a reprodução sexual — ausente no mtDNA, cada molécula é uma cópia integral da mãe) e a elevada taxa evolutiva em comparação ao DNA nuclear —, tornou-se um marcador molecular amplamente utilizado em diferentes contextos biológicos (BOORE, 1999; EKBLOM; WOLF, 2014).

O avanço das tecnologias de sequenciamento de nova geração (*Next-Generation Sequencing* — NGS) ampliou o acesso a grandes volumes de dados, favorecendo a obtenção de mitogenomas completos e de alta qualidade. No entanto, a reconstrução dessas moléculas a partir de dados brutos permanece um desafio técnico, especialmente no caso de **leituras curtas** (*reads* — os fragmentos individuais de sequência de DNA produzidos pelo sequenciador, tipicamente com 100 a 300 pares de bases), que exigem algoritmos sofisticados para superar problemas como regiões repetitivas e possíveis falhas de **montagem** — o processo computacional de reconstruir a sequência original do genoma a partir desses fragmentos sobrepostos. Nesse cenário, ferramentas especializadas, como o NOVOPlasty, têm desempenhado papel fundamental ao implementar estratégias de *seed-and-extend* capazes de recuperar organelas de forma eficiente a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017).

Embora a diversidade de softwares tenha contribuído para avanços significativos na montagem de organelas, a ausência de padronização metodológica ainda compromete a reprodutibilidade das análises. A multiplicidade de versões, dependências e configurações de ambiente gera obstáculos práticos para replicar experimentos, fenômeno amplamente conhecido como *dependency hell* (GRÜNING et al., 2018). Esse problema se insere no contexto mais amplo da crise de reprodutibilidade científica (BAKER, 2016), que tem motivado a adoção, na bioinformática, de práticas da engenharia de software, como a conteinerização de ambientes computacionais com Docker (MERKEL, 2014) e a orquestração de fluxos com sistemas de gerenciamento de workflows, a exemplo do Nextflow (DI TOMMASO et al., 2017).

Além da montagem propriamente dita, a anotação funcional dos genes mitocondriais constitui uma etapa igualmente relevante, que permite identificar os genes codificadores de proteínas, RNAs transportadores e RNAs ribossômicos presentes na molécula reconstruída. Ferramentas como o MITOS2 (DONATH et al., 2019) automatizam essa tarefa, mas sua execução — normalmente restrita a interfaces web — permanece desvinculada dos pipelines de montagem, criando uma descontinuidade no fluxo analítico.

Soma-se a esses desafios computacionais uma barreira de infraestrutura frequentemente subestimada. Datasets de sequenciamento de genoma total alcançam dezenas a centenas de gigabytes — o dataset da arara-azul-de-lear (*Anodorhynchus leari*) utilizado neste trabalho, por exemplo, totaliza 87 GB em formato FASTQ descomprimido. Processar esse volume na íntegra exigiria centenas de gigabytes de armazenamento temporário e memória RAM proporcional, recursos indisponíveis em computadores pessoais típicos de estudantes e pequenos laboratórios. Sem estratégias inteligentes de redução de dados, a montagem de mitogenomas a partir de dados públicos ficaria restrita a grupos com acesso a servidores dedicados.

### 1.2 Problema de Pesquisa

Diante do cenário apresentado — marcado pela fragmentação das ferramentas de montagem e anotação, pela ausência de padronização metodológica que compromete a reprodutibilidade, pela desconexão entre etapas analíticas que obriga o uso de interfaces web manuais, e pela barreira computacional imposta pelo volume dos dados de sequenciamento —, este trabalho busca responder à seguinte questão:

**De que maneira é possível integrar, em um único fluxo de trabalho automatizado, as etapas de aquisição, controle de qualidade, montagem *de novo*, anotação funcional e preparação para submissão ao GenBank de genomas mitocondriais, garantindo reprodutibilidade, portabilidade e viabilidade computacional em hardware doméstico?**

### 1.3 Objetivos

#### 1.3.1 Objetivo Geral

Desenvolver e validar um pipeline de bioinformática automatizado, reprodutível e portável para a montagem *de novo* e anotação funcional de genomas mitocondriais, utilizando tecnologias de conteinerização (Docker) e sistemas de gerenciamento de workflows (Nextflow).

#### 1.3.2 Objetivos Específicos

a) Containerizar, por meio de Dockerfiles com versões pinadas, cada uma das ferramentas bioinformáticas necessárias ao pipeline (SRA-Toolkit, FastQC, Trim Galore, NOVOPlasty e MITOS2).

b) Desenvolver o script de orquestração do pipeline utilizando Nextflow DSL2, definindo os processos, suas entradas, saídas e interdependências em módulos independentes.

c) Integrar os containers Docker ao workflow, garantindo que cada etapa da análise seja executada em seu próprio ambiente isolado e pré-configurado.

d) Implementar um módulo de análise piloto de qualidade (Pilot QC) capaz de determinar automaticamente o volume ideal de dados a serem processados, com base na qualidade das leituras e na fração mitocondrial estimada.

e) Integrar a anotação funcional automatizada (MITOS2) ao pipeline, incluindo a geração de estruturas secundárias de tRNAs e rRNAs via ViennaRNA.

f) Automatizar a geração de entregáveis em formato GenBank Flat File, com tratamento adequado de particularidades gênicas como o *frameshift* do ND3 em aves.

g) Validar o pipeline automatizado por meio de estudos de caso com dados públicos de sequenciamento de pelo menos duas espécies de grupos **taxonômicos** distintos — isto é, de posições diferentes na classificação biológica hierárquica (espécie, gênero, família, ordem, etc.) —, comparando os genomas montados e anotados com resultados de referência.

h) Disponibilizar o pipeline completo e sua documentação em repositório público no GitHub, assegurando sua acessibilidade e reprodutibilidade pela comunidade científica.

### 1.4 Justificativa

A análise de mitogenomas representa um componente essencial para a biologia evolutiva, a genética de populações e a sistemática, dado que essas moléculas fornecem informações relevantes sobre relações de parentesco, variação intra e interespecífica e processos de conservação. Sua relevância científica foi consolidada em trabalhos seminais que destacaram o papel estrutural e evolutivo dos mitogenomas como marcadores moleculares (BOORE, 1999), em estudos que enfatizaram sua aplicabilidade em investigações de genética de populações e evolução molecular (EKBLOM; WOLF, 2014) e em contribuições mais recentes que demonstraram seu potencial em análises de biodiversidade e conservação por meio de avanços metodológicos em montagem e anotação de organelas (ULIANO-SILVA et al., 2023).

Embora a obtenção de mitogenomas completos tenha sido amplamente favorecida pelo advento do NGS, a execução de pipelines de montagem ainda enfrenta desafios técnicos significativos. A diversidade de ferramentas, versões e dependências envolvidas gera um cenário de difícil replicação, conhecido como *dependency hell*, que compromete a reprodutibilidade e a padronização das análises (SANDVE et al., 2013; BAKER, 2016; GRÜNING et al., 2018).

Esse problema se torna especialmente relevante quando aplicado ao estudo de espécies ameaçadas de extinção, nas quais a disponibilidade de amostras biológicas é limitada e a confiabilidade dos resultados genômicos é crucial para fundamentar decisões de conservação. É o caso da **arara-azul-de-lear** (*Anodorhynchus leari*), ave endêmica da Caatinga nordestina do Brasil, classificada como "Em Perigo" pela IUCN Red List. Com população estimada em pouco mais de 1.700 indivíduos na natureza (ICMBio, 2022), a espécie depende de dados genômicos acessíveis e validáveis para subsidiar programas de manejo, reprodução em cativeiro e monitoramento populacional. Embora seu mitogenoma já tenha sido previamente reconstruído em atividade acadêmica, essa montagem foi realizada de forma manual, sem padronização metodológica e sem disponibilização pública do fluxo de trabalho — o que impede sua verificação e replicação por outros pesquisadores.

Um desafio adicional reside na completude dos pipelines existentes. A maioria das ferramentas e workflows publicados concentra-se exclusivamente na etapa de montagem, sem integrar a anotação funcional nem a preparação dos resultados para submissão a bancos de dados públicos como o GenBank. Essa lacuna obriga o pesquisador a recorrer a interfaces web externas (e.g., MITOS2 via UseGalaxy) e a conversões manuais de formato, introduzindo etapas não rastreáveis que comprometem a reprodutibilidade do ciclo completo de análise.

Nesse contexto, torna-se necessário o desenvolvimento de soluções que conciliem robustez biológica, consistência computacional e acessibilidade de hardware. A adoção de práticas modernas da engenharia de software, como a containerização e os sistemas de gerenciamento de workflows, tem se mostrado fundamental para superar limitações de reprodutibilidade, além de permitir que fluxos analíticos sejam portáveis, escaláveis e transparentes (DI TOMMASO et al., 2017; EWELS et al., 2020).

Assim, este trabalho justifica-se pela proposta de construir um pipeline automatizado e reprodutível que cubra o ciclo completo — desde a aquisição dos dados brutos até a geração de arquivos prontos para submissão ao GenBank (banco de dados público mantido pelo NCBI para depósito de sequências genéticas) —, unindo boas práticas de ciência aberta e engenharia computacional à relevância biológica dos mitogenomas. A aplicação ao caso concreto da arara-azul-de-lear demonstra não apenas a viabilidade técnica da abordagem, mas também sua contribuição direta para a conservação de uma espécie ameaçada da biodiversidade brasileira.

### 1.5 Metodologia

O presente trabalho adota uma abordagem de pesquisa aplicada, de natureza experimental, com desenvolvimento iterativo e incremental de software. O pipeline foi concebido a partir da experiência adquirida na disciplina Tópicos Especiais em Bioinformática A, ministrada pelo Prof. Jorge Estefano de Santana Souza, na qual um fluxo análogo de montagem mitogenômica foi utilizado de forma manual e não padronizada. O presente TCC formaliza, automatiza e expande substancialmente esse fluxo, integrando-o em um pipeline reprodutível implementado em Nextflow DSL2 e containerizado com Docker, seguindo boas práticas de engenharia de software científico (WILSON et al., 2014).

A validação foi realizada por meio de dois estudos de caso com dados públicos do NCBI Sequence Read Archive (SRA): (i) *Diploprion bifasciatus* (Perciformes — peixe-sabão barrado), como espécie-controle para verificação funcional do pipeline; e (ii) arara-azul-de-lear (*Anodorhynchus leari*, Psittaciformes — aves), como espécie-alvo para a aplicação final. Em ambos os casos, os mitogenomas foram montados *de novo* a partir de leituras Illumina de genoma total e comparados com referências disponíveis no GenBank.

O pipeline integra sete etapas sequenciais automatizadas: (1) análise piloto de qualidade para dimensionamento do volume de dados; (2) download e truncamento dos dados do SRA; (3) controle de qualidade com FastQC; (4) remoção de adaptadores e bases de baixa qualidade com Trim Galore; (5) montagem *de novo* com NOVOPlasty; (6) anotação funcional com MITOS2; e (7) compilação automática de 14 categorias de entregáveis, incluindo arquivos em formato GenBank Flat File prontos para submissão. Os detalhes de cada etapa, incluindo ferramentas, parâmetros e versões utilizadas, são apresentados no Capítulo 4.

### 1.6 Estrutura do Trabalho

O presente trabalho está organizado em sete capítulos. Este primeiro capítulo apresentou a contextualização do tema, o problema de pesquisa, os objetivos, a justificativa e um breve resumo metodológico.

O **Capítulo 2** reúne o referencial teórico, abordando os fundamentos da montagem *de novo* de genomas, as especificidades do DNA mitocondrial, a crise de reprodutibilidade em ciência computacional, a tecnologia de conteinerização com Docker e os sistemas de gerenciamento de workflows científicos.

O **Capítulo 3** analisa os trabalhos relacionados, revisando as principais ferramentas disponíveis para montagem de organelas — MitoHiFi, MITObim, GetOrganelle, MToolBox e NOVOPlasty — e apresenta uma análise comparativa entre elas.

O **Capítulo 4** detalha a metodologia, descrevendo o processo de desenvolvimento de software adotado, a visão geral e o fluxo de execução do pipeline, sua arquitetura modular, as ferramentas e tecnologias utilizadas, e cada uma das sete etapas da análise.

O **Capítulo 5** apresenta os resultados obtidos nas execuções do pipeline e a discussão, incluindo as métricas de montagem e anotação para ambas as espécies, a análise de desempenho computacional e uma avaliação da robustez da abordagem.

O **Capítulo 6** é dedicado integralmente à arara-azul-de-lear, apresentando a caracterização detalhada de seu mitogenoma, a análise dos genes codificadores de proteínas, RNAs transportadores e ribossomais, a comparação com a espécie congênere *Anodorhynchus hyacinthinus* e as implicações dos resultados para a conservação da espécie.

O **Capítulo 7** encerra o trabalho com as considerações finais, sintetizando as contribuições alcançadas e indicando perspectivas para trabalhos futuros.

---

## Capítulo 2 — Referencial Teórico

O presente capítulo reúne os fundamentos conceituais necessários para compreender a proposta metodológica deste trabalho. São abordados, inicialmente, os princípios da montagem *de novo* de genomas e as especificidades do DNA mitocondrial, que justificam sua relevância como objeto de estudo. Em seguida, discute-se a questão da reprodutibilidade em ciência computacional e as soluções contemporâneas oferecidas pela engenharia de software, com destaque para a tecnologia de conteinerização. Por fim, apresentam-se os sistemas de gerenciamento de workflows científicos, cuja integração com containers constitui a base para a construção de pipelines modernos, portáveis e reprodutíveis.

### 2.1 Montagem *De Novo* de Genomas e a Análise Mitocondrial

Em bioinformática, a **montagem de genomas** (*genome assembly*) consiste no processo computacional de reconstruir a sequência original de DNA a partir de fragmentos obtidos por tecnologias de Sequenciamento de Nova Geração (NGS). O conceito é análogo à reconstrução de um texto completo a partir de milhões de trechos sobrepostos — cada trecho sendo uma **leitura** (*read*), um fragmento curto de sequência gerado pelo sequenciador. O método mais amplamente utilizado para gerar esses fragmentos é o *shotgun sequencing*, no qual o DNA genômico é clivado aleatoriamente em milhões de pequenos pedaços que, posteriormente, são sequenciados. Essa estratégia garante alta **cobertura** do genoma — isto é, cada posição do DNA é representada por múltiplas leituras, tipicamente dezenas a centenas de vezes —, mas resulta em dados fragmentados que precisam ser cuidadosamente alinhados e sobrepostos para permitir a reconstrução da sequência completa (EKBLOM; WOLF, 2014). A complexidade computacional dessa tarefa motivou o desenvolvimento de algoritmos especializados, baseados em técnicas como grafos de De Bruijn — estrutura em que cada nó representa uma subsequência de $k$ nucleótidos consecutivos (denominada **k-mer**), e cada aresta representa sobreposições entre k-mers adjacentes — e estratégias de sobreposição-layout-consenso (MILLER; KOREN; SUTTON, 2010).

> **Figura 1** — Diagrama ilustrativo do processo de *shotgun sequencing*: clivagem aleatória do DNA genômico em fragmentos, sequenciamento em leituras curtas e sobreposição para reconstrução da sequência original.

Conforme ilustrado na Figura 1, o processo envolve três estágios principais: a fragmentação aleatória do DNA genômico, o sequenciamento paralelo dos milhões de fragmentos resultantes e a reconstrução computacional da sequência original por meio da identificação de sobreposições entre as leituras. É essa última etapa — a montagem — que constitui o desafio computacional central abordado neste trabalho.

Na montagem genômica, fragmentos sobrepostos são inicialmente agrupados em **contigs** (*contiguous sequences* — sequências contínuas reconstruídas a partir do alinhamento de leituras sobrepostas, análogas a blocos de texto já montados). Quando existe informação adicional que permite inferir a ordem relativa entre diferentes contigs — como dados de **leituras pareadas** (*paired-end reads*: pares de leituras gerados a partir das duas extremidades do mesmo fragmento de DNA, com distância conhecida entre elas) —, esses podem ser organizados em estruturas mais amplas chamadas **scaffolds**. Assim, enquanto os contigs representam sequências consolidadas diretamente a partir das leituras, os scaffolds constituem hipóteses mais abrangentes sobre a organização genômica, servindo como passo intermediário até a obtenção de montagens completas e anotadas.

As tecnologias de NGS diferem substancialmente no comprimento das leituras produzidas. Sequenciadores Illumina, por exemplo, geram leituras curtas de alta acurácia (50 a 300 pares de bases), que exigem algoritmos sofisticados para lidar com regiões repetitivas. Em contrapartida, tecnologias como PacBio e Oxford Nanopore produzem leituras longas, que podem alcançar milhares de pares de bases em uma única sequência, favorecendo a continuidade da montagem, ainda que a custos mais elevados e com menor rendimento global (EKBLOM; WOLF, 2014).

> **Figura 2** — Comparativo das principais tecnologias de sequenciamento: Illumina (leituras curtas, alta acurácia), PacBio HiFi (leituras longas, alta fidelidade) e Oxford Nanopore (leituras ultra-longas, menor acurácia).

Como evidenciado na Figura 2, a escolha da tecnologia de sequenciamento envolve um compromisso entre comprimento das leituras, acurácia e custo. O presente trabalho utiliza dados Illumina (leituras curtas de alta acurácia), que constituem a grande maioria dos datasets públicos disponíveis e são a entrada padrão das ferramentas de montagem de organelas mais utilizadas.

O DNA mitocondrial (mtDNA) reúne um conjunto de particularidades que o tornam um alvo privilegiado para a montagem *de novo*. Antes de descrevê-las, convém definir conceitos fundamentais da biologia molecular que serão recorrentes ao longo deste trabalho.

Um **gene** é um trecho de DNA que contém as instruções para a síntese de uma molécula funcional — geralmente uma **proteína**, macromolécula que executa a maioria das funções celulares (catálise enzimática, transporte de substâncias, estrutura celular). O DNA é composto por uma sequência de **nucleotídeos**, unidades químicas identificadas pelas letras A (adenina), T (timina), G (guanina) e C (citosina). Dois nucleotídeos complementares em fitas opostas (A com T, G com C) formam um **par de bases** (bp, do inglês *base pair*), a unidade de medida padrão para comprimento de sequências genômicas. A informação genética contida no DNA é lida em blocos de três nucleotídeos consecutivos chamados **códons**: cada códon especifica um **aminoácido**, os “blocos de construção” das proteínas (existem 20 aminoácidos padrão, como fenilalanina, valina, isoleucina, entre outros). Um gene que codifica uma proteína possui um **códon de início** (tipicamente ATG, que sinaliza onde a tradução deve começar) e um **códon de término** ou **códon de parada** (como TAA, TAG ou TGA, que sinaliza onde a tradução deve parar). O fluxo da informação genética segue o chamado **dogma central da biologia molecular**: o DNA é transcrito em **RNA mensageiro** (mRNA), que por sua vez é traduzido em proteína pelos **ribossomos** — complexos macromoleculares compostos por proteínas e **RNA ribossomal** (rRNA). Nesse processo de tradução, cada códon do mRNA é reconhecido por uma molécula de **RNA transportador** (tRNA), que carrega o aminoácido correspondente e o posiciona no ribossomo para ser incorporado à cadeia proteica em crescimento. A região do tRNA que reconhece o códon é chamada **anticódon** — uma sequência de três nucleotídeos complementar ao códon do mRNA (ALBERTS et al., 2022).

Em metazoários (animais multicelulares), o mtDNA corresponde a uma molécula **circular** de fita dupla — diferentemente dos **cromossomos** nucleares (estruturas lineares de DNA compactado localizadas no núcleo da célula, que carregam a maior parte da informação genética), que são lineares, o DNA mitocondrial forma um anel contínuo sem extremidades —, com tamanho geralmente entre 12 e 22 mil pares de bases, codificando 37 genes distribuídos em 13 proteínas associadas à **fosforilação oxidativa** (o conjunto de reações bioquímicas que a mitocôndria utiliza para produzir ATP, a molécula de energia da célula), 22 tRNAs e 2 rRNAs (ANDERSON et al., 1981; ASAKAWA et al., 1995; BOORE, 1999). Em sua revisão seminal, Boore (1999) observou que, embora a economia extrema de tamanho seja aparente em muitos casos, essa visão — até então quase axiomática — é refutada por mitogenomas com grandes quantidades de sequência não codificante descobertos em alguns artrópodes, moluscos e nematódeos, demonstrando que as pressões seletivas sobre o tamanho do mtDNA variam substancialmente entre linhagens. O sequenciamento completo do genoma mitocondrial humano por Anderson et al. (1981) representou um marco histórico, estabelecendo o modelo de referência para a organização gênica mitocondrial dos metazoários.

> **Figura 3** — Representação circular esquemática de um genoma mitocondrial típico de vertebrado, mostrando a distribuição dos 37 genes (13 CDS, 22 tRNAs, 2 rRNAs), a região controle (D-loop) e as fitas pesada (H) e leve (L).

Na Figura 3, observa-se a distribuição dos 37 genes ao longo das duas fitas da molécula circular: a fita pesada (H, externa) codifica a maioria dos genes, enquanto a fita leve (L, interna) contém o gene ND6 e oito tRNAs. A região controle (D-loop), que não codifica proteínas, é o único segmento não codificante de tamanho significativo, responsável pela regulação da replicação e transcrição do mtDNA.

A ausência de **íntrons** (segmentos não codificantes que, nos genes nucleares, intercalam as regiões codificantes e precisam ser removidos antes da tradução), a estrutura compacta e a elevada taxa de cópias por célula tornam o mtDNA especialmente acessível, facilitando sua recuperação mesmo em experimentos voltados ao genoma nuclear, nos quais ele aparece como subproduto. Além disso, características funcionais como a herança predominantemente materna, a ausência geral de recombinação e a taxa evolutiva superior à do DNA nuclear ampliam sua utilidade em estudos biológicos. Essas propriedades explicam sua ampla aplicação em análises de genética de populações, **filogeografia** (estudo da distribuição geográfica de linhagens genéticas), sistemática e conservação da biodiversidade.


Para explorar essas vantagens, foram desenvolvidos algoritmos específicos de montagem. Uma das abordagens mais empregadas é a *seed-and-extend* (semente-e-extensão), em que uma sequência inicial conhecida (*seed* — semente) — que pode ser um gene ou até mesmo um fragmento de organela de espécie próxima — serve como ponto de partida para estender iterativamente a montagem até que a molécula circular seja reconstruída. Quando as extremidades da sequência em construção se sobrepõem, confirmando que o genoma forma um anel contínuo, diz-se que a montagem foi **circularizada** — indicação de que o genoma mitocondrial completo foi recuperado com sucesso. O NOVOPlasty é um dos principais programas baseados nessa estratégia, destacando-se pela eficiência na recuperação de genomas mitocondriais e cloroplastidiais a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017). Avanços recentes, como o pipeline MitoHiFi, exploram leituras longas de alta fidelidade, evidenciando a constante evolução das ferramentas voltadas à montagem de organelas (ULIANO-SILVA et al., 2023).

> **Figura 4** — Diagrama ilustrativo da estratégia *seed-and-extend*: (a) sequência semente; (b) extensão iterativa a partir de leituras sobrepostas; (c) circularização quando as extremidades se encontram.

A Figura 4 evidencia os três estágios do algoritmo: (a) a seleção de uma sequência semente conhecida, como o gene *cox1* de uma espécie filogeneticamente próxima; (b) a extensão iterativa, em que o algoritmo busca leituras que se sobreponham às extremidades da sequência em construção; e (c) a circularização, momento em que as duas extremidades da sequência se encontram, confirmando a recuperação completa do genoma mitocondrial.

Complementarmente à montagem, a **anotação funcional** dos genes mitocondriais permite identificar e classificar os elementos codificantes da molécula reconstruída — ou seja, determinar quais genes estão presentes, suas posições exatas e suas funções biológicas. Para metazoários, ferramentas como o MITOS2 combinam alinhamentos de homologia (comparações com genes conhecidos em bancos de dados) com modelos de covariância (Infernal) para predizer genes codificadores de proteínas, tRNAs e rRNAs, produzindo anotações em formatos padronizados como **GFF3** (*General Feature Format version 3* — formato tabular que especifica a localização e o tipo de cada gene no genoma) (DONATH et al., 2019). A integração da etapa de anotação em pipelines automatizados de montagem permanece, contudo, uma lacuna na maioria das soluções disponíveis na literatura. A ampla adoção do MITOS2 pela comunidade científica — empregado como ferramenta primária de anotação em estudos recentes abrangendo vertebrados (ZHAN et al., 2024; KUNDU et al., 2024), invertebrados terrestres (SHI et al., 2024; TAO et al., 2024) e invertebrados marinhos (WEI et al., 2024; ALBOASUD; JEONG; LEE, 2024) — atesta a confiabilidade de seus resultados e reforça a pertinência de sua integração em fluxos automatizados.

### 2.2 Reprodutibilidade em Ciência Computacional

A reprodutibilidade é um pilar do método científico, consistindo na capacidade de um pesquisador independente replicar um experimento, com os mesmos materiais e métodos, e obter resultados consistentes. Em pesquisas computacionais, isso se traduz na possibilidade de executar o mesmo código e software sobre os mesmos dados e chegar ao mesmo resultado (PENG, 2011; SANDVE et al., 2013).

Contudo, a ciência moderna enfrenta o que tem sido amplamente denominado como uma "crise de reprodutibilidade". Em uma pesquisa abrangente com 1.576 pesquisadores de diversas áreas, Baker (2016) revelou que 90% dos entrevistados reconhecem a existência dessa crise — sendo que 52% a consideram significativa e 38% a classificam como leve. Uma barreira crítica identificada por Peng (2011) é que, em muitos casos, o código computacional utilizado nas análises simplesmente não está mais disponível: softwares interativos frequentemente não registram as ações dos usuários de forma concreta, e mesmo quando scripts são utilizados, a combinação de múltiplos pacotes raramente é preservada de modo reprodutível. No domínio da bioinformática, essa crise é particularmente acentuada devido à alta complexidade dos pipelines de análise. Um único pipeline pode envolver dezenas de ferramentas de software distintas, cada uma desenvolvida por diferentes grupos, em diferentes linguagens de programação e com um conjunto único e, por vezes, conflitante de dependências de bibliotecas e do próprio sistema operacional. A tarefa de instalar e configurar manualmente esse ecossistema de software para replicar uma análise é extremamente suscetível a erros, uma condição frequentemente descrita como dependency hell (inferno de dependências) (GRÜNING et al., 2018).

Pequenas variações na versão de uma ferramenta ou de uma biblioteca subjacente podem levar a resultados drasticamente diferentes, comprometendo a validade e a confiabilidade da pesquisa. Essa barreira técnica não apenas dificulta a verificação dos resultados, mas também impede a reutilização e a adaptação de métodos computacionais, contrariando os princípios FAIR (*Findable, Accessible, Interoperable, and Reusable*), que visam maximizar o valor dos dados e das análises científicas ao enfatizar, em particular, a capacidade de máquinas encontrarem e utilizarem dados automaticamente (SANDVE et al., 2013; WILKINSON et al., 2016; GRÜNING et al., 2018; BRITISH ECOLOGICAL SOCIETY, 2017).

> **Figura 5** — Diagrama comparativo do "inferno de dependências": à esquerda, instalação manual com conflitos entre versões; à direita, pipeline containerizado com Docker, cada ferramenta em seu container isolado.

O painel esquerdo da Figura 5 ilustra o cenário típico de instalação manual, em que diferentes ferramentas exigem versões conflitantes das mesmas bibliotecas, tornando a configuração do ambiente um obstáculo frequentemente intransponível para pesquisadores sem formação em TI. O painel direito mostra a solução containerizada, na qual cada ferramenta é encapsulada em seu próprio ambiente isolado, eliminando conflitos de dependências e garantindo execução idêntica em qualquer máquina.

### 2.3 Tecnologia de Conteinerização com Docker

Como solução para o problema do dependency hell e para garantir a reprodutibilidade computacional, a comunidade científica tem adotado práticas e tecnologias da engenharia de software, com destaque para a virtualização em nível de sistema operacional, também conhecida como containerização. A plataforma Docker se estabeleceu como a principal ferramenta desse paradigma, permitindo o empacotamento de uma aplicação e de todo o seu ambiente de execução - incluindo bibliotecas, códigos-fonte e dependências - em uma unidade padronizada e isolada chamada container (MERKEL, 2014; GRÜNING et al., 2018).

O funcionamento do Docker baseia-se em dois conceitos centrais: a imagem e o container. A imagem é um template estático e imutável, construído a partir de um arquivo de texto com instruções sequenciais chamado Dockerfile. Esse arquivo serve como uma "receita" que define o sistema operacional base, as dependências a serem instaladas e os comandos a serem executados, criando um ambiente de software completo e autocontido. O container, por sua vez, é uma instância executável e isolada de uma imagem. Ao contrário de máquinas virtuais tradicionais, os containers compartilham o kernel do sistema operacional hospedeiro, o que os torna leves e eficientes para iniciar (BRITISH ECOLOGICAL SOCIETY, 2017).

Essa arquitetura garante que um software containerizado se comporte de maneira idêntica em qualquer ambiente que possua o Docker instalado - seja em computadores pessoais, servidores de produção ou ambientes de nuvem - resolvendo de forma eficaz o problema da inconsistência de ambientes e sendo um passo fundamental para a criação de pipelines científicos portáveis e reprodutíveis (DI TOMMASO et al., 2017).

> **Figura 6** — Comparação entre containers Docker e máquinas virtuais (VM): Dockerfile (receita textual), imagem (template imutável) e container (instância em execução), evidenciando a leveza dos containers.

A diferença fundamental evidenciada na Figura 6 reside na ausência da camada de sistema operacional convidado nos containers: enquanto cada máquina virtual embarca um sistema operacional completo (Guest OS), os containers compartilham o kernel do hospedeiro, resultando em imagens na ordem de megabytes (versus gigabytes para VMs) e tempos de inicialização de segundos (versus minutos). Essa eficiência torna a containerização particularmente adequada para pipelines de bioinformática, nos quais dezenas de ferramentas distintas precisam ser executadas sequencialmente.

### 2.4 Sistemas de Gerenciamento de Workflows Científicos

A containerização resolve a questão do ambiente de software, mas não a orquestração das múltiplas etapas que compõem uma análise bioinformática. A execução manual e sequencial de cada container é ineficiente e propensa a erros. Para gerenciar essa complexidade, foram desenvolvidos Sistemas de Gerenciamento de Workflows (Workflow Management Systems), como o Nextflow (DI TOMMASO et al., 2017) e o Snakemake (KÖSTER; RAHMANN, 2012). Essas plataformas permitem que os pesquisadores definam pipelines computacionais complexos por meio de linguagens de alto nível.

Nesses sistemas, um workflow é decomposto em uma série de processos ou regras interconectadas. Cada processo define uma tarefa computacional discreta, especificando seus comandos, suas entradas (arquivos ou variáveis) e suas saídas. O gerenciador de workflow é então responsável por monitorar o estado dos dados e executar os processos na ordem correta, automatizando o fluxo de informação entre as etapas.

A principal vantagem dessas ferramentas reside na abstração da camada de execução. Elas integram-se nativamente com tecnologias como Docker, permitindo que cada processo seja executado em seu próprio container pré-configurado. Além disso, são compatíveis com diversos ambientes de execução, desde computadores pessoais até clusters de computação de alto desempenho (HPC) e plataformas de nuvem. Isso confere ao pipeline características essenciais para a ciência moderna: portabilidade, escalabilidade e retomada automática (resumability) em caso de falhas, evitando o reprocessamento de etapas já concluídas.

A combinação de um gerenciador de workflows com a containerização é hoje considerada o padrão-ouro para o desenvolvimento de pipelines de bioinformática, sendo a base de iniciativas comunitárias como o nf-core, que busca padronizar e compartilhar pipelines que seguem as melhores práticas de reprodutibilidade (EWELS et al., 2020).

O referencial teórico expôs os conceitos fundamentais que dão suporte à proposta deste trabalho, abordando desde os princípios de montagem de novo de genomas até as soluções contemporâneas de reprodutibilidade, como a containerização e os sistemas de gerenciamento de workflows. Esses conceitos fornecem a base conceitual sobre a qual se apoiam os trabalhos relacionados, discutidos no Capítulo 3, que permitem situar a presente proposta no estado da arte da bioinformática aplicada à montagem de mitogenomas.

---

## Capítulo 3 — Trabalhos Relacionados

Nesta seção, são apresentados e analisados trabalhos recentes da literatura que, assim como este projeto, propõem soluções computacionais para a montagem de genomas de organelas. O objetivo desta análise comparativa é contextualizar a presente proposta no estado da arte da área, identificando as abordagens metodológicas existentes e destacando os diferenciais e a contribuição original do pipeline a ser desenvolvido neste trabalho.

### 3.1 MitoHiFi: Montagem com Leituras Longas

Um trabalho de destaque na área é o MitoHiFi, um pipeline desenvolvido para a montagem de genomas mitocondriais a partir de dados de **leituras longas** de alta fidelidade (PacBio HiFi — tecnologia que produz leituras de milhares de pares de bases com alta precisão, em contraste com as leituras curtas de 100–300 bp geradas pela plataforma Illumina). Proposto por Uliano-Silva et al. (2023), o objetivo principal era criar um método que aproveitasse a extensão e a precisão deste tipo de dado para gerar mitogenomas completos e corretamente circularizados com maior eficiência. A metodologia do MitoHiFi integra ferramentas como o hifiasm-meta para montagem e o GetOrganelle para a etapa de finalização, sendo orquestrada por meio de scripts em Shell e utilizando **Conda** (gerenciador de pacotes e ambientes, alternativa ao Docker para gestão de dependências de software) para o gerenciamento das dependências.


A abordagem do MitoHiFi é extremamente atual por focar em uma tecnologia de sequenciamento de ponta. Seus resultados são expressivos: o pipeline foi utilizado para montar 374 mitogenomas (368 Metazoa e 6 Fungi) em projetos como o Darwin Tree of Life, o Vertebrate Genomes Project e o Aquatic Symbiosis Genome Project, além de identificar automaticamente variantes decorrentes de heteroplasmia e inserções nucleares de sequências mitocondriais (NUMTs) (ULIANO-SILVA et al., 2023). No entanto, sua implementação via scripts e Conda, embora eficaz para seu propósito, pode apresentar limitações de portabilidade e reprodutibilidade quando comparada a soluções baseadas em containers e sistemas de workflow. O diferencial do presente TCC reside precisamente neste ponto: enquanto o MitoHiFi inova na aplicação de um novo tipo de dado biológico, nossa proposta foca em uma arquitetura de software mais robusta, utilizando Docker e Nextflow, para garantir que o pipeline de montagem (baseado em leituras curtas) seja universalmente reprodutível e portável.

### 3.2 MITObim: Estratégia de Baiting and Iterative Mapping

O MITObim, desenvolvido por Hahn et al. (2013), foi um dos primeiros métodos amplamente utilizados para montagem de organelas a partir de dados de NGS. Sua estratégia de *baiting and iterative mapping* (captura e mapeamento iterativo) consiste em iniciar o processo com uma sequência guia (um fragmento do genoma alvo ou de uma espécie próxima) e, a partir dela, recuperar leituras homólogas (isto é, leituras suficientemente similares à sequência guia para corresponderem ao genoma mitocondrial), refinando o alinhamento de forma iterativa até a montagem completa do mitogenoma.


Em testes com dados reais de NGS, o MITObim demonstrou capacidade de recuperar mitogenomas completos com acurácia superior a 99,5% em menos de 24 horas utilizando um computador desktop padrão (HAHN; BACHMANN; CHEVREUX, 2013). Apesar desses resultados e de sua relevância histórica como pioneiro na abordagem, o MITObim apresenta limitações em termos de eficiência computacional e escalabilidade, especialmente quando aplicado a conjuntos de dados maiores ou de maior complexidade, vindo a ser gradualmente substituído por ferramentas mais modernas, como o NOVOPlasty e o GetOrganelle.

### 3.3 GetOrganelle: Montagem Baseada em Grafos

O GetOrganelle, proposto por Jin et al. (2020), é uma das ferramentas mais recentes e populares para montagem de organelas. Ele utiliza **grafos de montagem** (estruturas matemáticas em que os nós representam fragmentos de sequência e as arestas representam sobreposições entre eles) derivados de dados de NGS para isolar e reconstruir genomas de organelas, oferecendo resultados robustos tanto para mitocôndrias quanto para cloroplastos (organelas fotossintetizantes de plantas). Em uma avaliação sistemática com 50 datasets publicados de plantas, o GetOrganelle foi capaz de remontar plastomas circulares em 47 dos 50 casos — uma taxa de sucesso significativamente superior à do NOVOPlasty. Além disso, a avaliação por mapeamento de reads demonstrou que os plastomas montados pelo GetOrganelle são geralmente mais acurácios que os publicados originalmente ou remontados pelo NOVOPlasty (JIN et al., 2020). A ferramenta tem se destacado também por sua capacidade de gerar todas as configurações possíveis quando os genomas apresentam configurações flip-flop ou outros isômeros mediados por repetições, oferecendo resultados consistentes mesmo em organismos com maior complexidade genômica.


No entanto, assim como outras ferramentas tradicionais, o GetOrganelle depende de um ambiente de software adequado para ser executado, geralmente configurado manualmente pelo usuário ou via gerenciadores como Conda, o que ainda pode limitar sua reprodutibilidade plena em diferentes contextos computacionais. Cabe notar que, embora o GetOrganelle apresente taxas de sucesso e acurácia superiores ao NOVOPlasty em plastomas de plantas, a literatura recente demonstra que o NOVOPlasty permanece como a ferramenta predominante para montagem de mitogenomas de metazoários — sendo adotado em estudos recentes de larga escala com 213 mitogenomas de lagartos (ZHAN et al., 2024) e em trabalhos individuais de peixes, arraias, moluscos e invertebrados (KUNDU et al., 2024; GUERREIRO et al., 2025; SHI et al., 2024; TAO et al., 2024). Essa preferência pela comunidade reforça a escolha do NOVOPlasty como montador neste pipeline.

### 3.4 MToolBox: Foco em Genomas Mitocondriais Humanos

O MToolBox, apresentado por Calabrese et al. (2014), foi uma das primeiras pipelines integradas voltadas especificamente para o estudo de genomas mitocondriais humanos. Ele combina montagem, anotação e análise funcional em um único fluxo, permitindo identificar **variantes mitocondriais** (alterações pontuais na sequência de DNA que distinguem indivíduos ou populações) relevantes para estudos biomédicos e populacionais.

Apesar de sua contribuição para o campo, o MToolBox é fortemente orientado ao contexto humano e não é tão amplamente aplicável a outros organismos. Além disso, por ter sido desenvolvido antes da adoção massiva de containers e workflows, carece de recursos de portabilidade e reprodutibilidade presentes em soluções mais modernas.


### 3.5 Análise Comparativa

A Tabela a seguir sintetiza as características das ferramentas e pipelines analisados, confrontando-as com a proposta deste trabalho. A comparação considera critérios relevantes para a reprodutibilidade, portabilidade e completude do fluxo analítico.

| Critério | MitoHiFi | MITObim | GetOrganelle | MToolBox | NOVOPlasty | **Este trabalho** |
|---|---|---|---|---|---|---|
| **Tipo de dado** | Long reads (HiFi) | Short reads | Short/Long reads | Short reads | Short reads | Short reads |
| **Montagem** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Anotação funcional** | — | — | — | ✓ (humano) | — | ✓ (metazoários) |
| **Visualização (mapa circular)** | — | — | — | — | — | ✓ (Biopython) |
| **GenBank Flat File** | — | — | — | — | — | ✓ |
| **Containerização** | Parcial (Conda) | — | Conda | — | — | Docker (5 imagens) |
| **Workflow manager** | — | — | — | — | — | Nextflow DSL2 |
| **Pilot QC automático** | — | — | — | — | — | ✓ |
| **Estrutura secundária RNA** | — | — | — | — | — | ✓ (ViennaRNA) |
| **Execução em laptop** | Requer cluster/HPC | Servidor recomendado | Servidor recomendado | Desktop/servidor | Desktop | ✓ (notebook 16 GB) |
| **Código aberto** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

A análise comparativa evidencia que, embora cada ferramenta contribua com inovações relevantes para domínios específicos, nenhuma oferece um fluxo integrado que cubra todas as etapas — desde a aquisição dos dados brutos até a geração de arquivos prontos para submissão ao GenBank — em um ambiente completamente containerizado e orquestrado. Adicionalmente, nenhuma das soluções analisadas incorpora mecanismos de redução automática do volume de dados (como o Pilot QC e a extração seletiva via `fastq-dump -X`), o que as torna dependentes de hardware com grande capacidade de armazenamento e memória. A proposta deste trabalho preenche essas lacunas ao combinar montagem, anotação funcional, geração de entregáveis, práticas de reprodutibilidade e acessibilidade computacional em um único pipeline modular, executável em computadores pessoais.

Para complementar a análise tabular, a Figura 7 apresenta uma visualização em gráfico radar que sintetiza o desempenho relativo de cada ferramenta em cinco eixos qualitativos, pontuados em escala de 1 a 5 conforme os critérios da Tabela 3. Essa pontuação foi atribuída com base nas informações publicadas nos artigos originais de cada ferramenta e na análise direta da documentação e dos repositórios oficiais.

**Tabela 3 — Critérios de pontuação do gráfico radar comparativo (escala 1–5).**

| Eixo | Critério de pontuação |
|---|---|
| Reprodutibilidade | 1 = instalação manual; 2 = Conda; 3 = Conda + versões pinadas; 4 = Docker ou Singularity; 5 = Docker + workflow manager com versões pinadas |
| Completude do fluxo | 1 = apenas montagem; 2 = montagem + anotação parcial; 3 = montagem + anotação completa; 4 = inclui visualização ou conversão de formatos; 5 = fluxo completo (aquisição → montagem → anotação → entregáveis GenBank) |
| Acessibilidade computacional | 1 = requer cluster/HPC; 2 = servidor dedicado recomendado; 3 = desktop padrão; 4 = desktop com dados reduzidos manualmente; 5 = notebook com redução automática de dados |
| Automação | 1 = configuração manual extensiva; 2 = scripts isolados; 3 = pipeline com parâmetros manuais; 4 = pipeline com defaults razoáveis; 5 = execução com um único comando + decisões automáticas (Pilot QC, fallback k-mer) |
| Generalidade taxonômica | 1 = uma espécie (e.g., humano); 2 = um grupo taxonômico restrito; 3 = metazoários; 4 = metazoários + outros eucariotos; 5 = qualquer organismo com genoma mitocondrial |

> **Figura 7** — Gráfico radar comparando MitoHiFi, MITObim, GetOrganelle, MToolBox, NOVOPlasty e o pipeline proposto em 5 eixos: completude do fluxo, reprodutibilidade, acessibilidade computacional, automação e generalidade taxonômica.

Conforme ilustrado na Figura 7, o pipeline proposto neste trabalho (representado em vermelho) domina quatro dos cinco eixos avaliados — reprodutibilidade, completude do fluxo, automação e acessibilidade computacional —, atingindo a pontuação máxima (5) em cada um deles, enquanto as demais ferramentas concentram-se predominantemente na faixa de 1 a 3. O competidor mais próximo varia conforme o eixo considerado: o MitoHiFi alcança 4 em automação, o MToolBox atinge 3 em completude, e o GetOrganelle e o NOVOPlasty obtêm 4 em generalidade taxonômica — justamente o único eixo em que este trabalho não lidera, com pontuação 3 (escopo restrito a metazoários). Essa diferença, contudo, decorre de uma escolha deliberada de projeto: ao concentrar o pipeline em genomas mitocondriais de metazoários, foi possível otimizar cada etapa — do Pilot QC à geração do GenBank Flat File — para esse grupo taxonômico, em vez de buscar uma generalidade que diluísse a profundidade do fluxo. O padrão visual resultante é claro: enquanto as demais ferramentas formam polígonos compactos e irregulares no interior do gráfico, o pipeline proposto projeta-se de forma ampla em quase todos os eixos, evidenciando graficamente sua proposta de integração, reprodutibilidade e acessibilidade.

A análise dos trabalhos relacionados permitiu identificar diferentes abordagens para a montagem de organelas, incluindo pipelines baseados em leituras longas, estratégias iterativas e montadores especializados como o NOVOPlasty. Embora cada solução apresente contribuições relevantes, permanecem lacunas quanto à padronização, portabilidade, completude do fluxo analítico e reprodutibilidade. Nesse contexto, o Capítulo 4 descreve a metodologia adotada neste trabalho, detalhando a concepção do pipeline proposto e as ferramentas que o compõem.

---

## Capítulo 4 — Metodologia

Nesta seção, é apresentado o plano metodológico para o desenvolvimento do pipeline de montagem de mitogenomas. A concepção deste pipeline teve como ponto de partida o fluxo de análise mitogenômica utilizado na disciplina Tópicos Especiais em Bioinformática A, ministrada pelo Prof. Jorge Estefano de Santana Souza na Universidade Federal do Rio Grande do Norte (UFRN). Naquela disciplina, as etapas de download, controle de qualidade, montagem e anotação eram executadas manualmente, sem padronização de versões ou ambiente computacional. O presente trabalho parte dessa experiência para projetar uma solução automatizada e reprodutível, formalizando cada etapa em módulos containerizados e orquestrados. A abordagem combina etapas de análise bioinformática com práticas consolidadas de engenharia de software, seguindo um processo de desenvolvimento sistemático descrito na seção 4.1, antes de detalhar a arquitetura do pipeline e suas etapas constituintes.

### 4.1 Processo de Desenvolvimento de Software

O desenvolvimento do pipeline seguiu uma abordagem de **desenvolvimento iterativo e incremental**, metodologia consolidada na engenharia de software em que o sistema é construído em ciclos sucessivos, cada um produzindo uma versão funcional e testável do produto (LARMAN; BASILI, 2003). Essa abordagem se contrapõe ao modelo cascata (*waterfall*), no qual todas as fases — levantamento de requisitos, projeto, implementação e teste — são executadas sequencialmente e o produto só é validado ao final. No contexto de software científico, em que os requisitos frequentemente evoluem à medida que os resultados são analisados, o modelo iterativo é amplamente recomendado (WILSON et al., 2014).

O ciclo de desenvolvimento foi organizado em duas fases principais, cada uma constituindo uma iteração completa com entrega funcional:

- **Fase 1 — Montagem**: implementação dos módulos de download (SRA-Toolkit), controle de qualidade (FastQC, Trim Galore) e montagem (NOVOPlasty), com containerização em Docker e orquestração via Nextflow. Ao final desta fase, o pipeline já era capaz de produzir mitogenomas completos e circularizados a partir de um acesso SRA.
- **Fase 2 — Anotação e entregáveis**: adição dos módulos de análise piloto de qualidade (Pilot QC), anotação funcional (MITOS2), geração de estruturas secundárias de RNA (ViennaRNA/RNAplot), produção de arquivos para submissão ao GenBank (GenBank Flat File) e compilação automática de entregáveis.

Cada fase produziu um pipeline funcional e independente. A Fase 1 foi preservada como branch separada no repositório Git (`fase1-montagem`), permitindo uso e validação independente por terceiros — uma prática de versionamento semântico que garante rastreabilidade e facilita a manutenção (SOMMERVILLE, 2016).

> **Figura 8** — Diagrama de fases do desenvolvimento do pipeline: Fase 1 (Montagem) e Fase 2 (Anotação), com indicação dos ciclos iterativos e das espécies utilizadas em cada fase.

Conforme representado na Figura 8, cada fase constitui um ciclo completo de desenvolvimento–teste–correção–validação. A separação permite que a Fase 1 seja utilizada independentemente por pesquisadores interessados apenas na montagem, enquanto a Fase 2 adiciona as funcionalidades de anotação e geração de entregáveis sobre a base já validada.

#### 4.1.1 Validação por Espécie-Controle

Um elemento central da estratégia de verificação e validação (V&V) foi o emprego de uma **espécie-controle** com mitogenoma previamente depositado em banco de dados público. Antes de aplicar o pipeline à espécie-alvo (a arara-azul-de-lear), cuja sequência mitocondrial era inédita, o sistema foi integralmente validado com *Diploprion bifasciatus*, cujo mitogenoma de referência (PZ143763.1) está depositado no GenBank. Essa estratégia permitiu comparar a montagem obtida pelo pipeline com a sequência de referência conhecida, verificando a correção funcional do sistema em condições controladas antes de aplicá-lo a dados sem referência prévia.

No desenvolvimento de software científico, a validação com dados de referência conhecida é considerada uma das práticas mais importantes para garantir a confiabilidade dos resultados computacionais (WILSON et al., 2014). Essa abordagem é análoga ao conceito de *oracle testing* na engenharia de software, em que uma fonte externa de verdade é utilizada para validar as saídas do sistema sob teste (SOMMERVILLE, 2016).

#### 4.1.2 Tratamento de Erros e Correções Iterativas

O modelo iterativo também se manifestou na resolução de problemas técnicos encontrados durante a execução do pipeline. Exemplos concretos incluem:

1. **Incompatibilidade do SRA-Toolkit**: a descoberta de que o `fasterq-dump` do SRA-Toolkit 3.x não implementa o flag `-X` para limitação de leituras na origem. Inicialmente, a solução adotada foi a truncagem pós-conversão via `head`, mas a execução com o dataset metagenômico de *A. hyacinthinus* (145,8M spots, ~124 GB de FASTQ) revelou que essa abordagem é inviável para grandes volumes. A correção definitiva substituiu o `fasterq-dump` pelo `fastq-dump -X` quando `sra_max_reads` está definido, extraindo diretamente apenas os primeiros *N* spots sem gravar o dataset completo (Seção 5.9);
2. **Flag do RNAplot**: a incompatibilidade do flag `--output-format=svg` na versão empacotada do ViennaRNA, corrigida para a sintaxe `-f svg` aceita pela versão instalada;
3. **Frameshift do gene ND3**: a detecção, durante a análise dos resultados do MITOS2, de que o gene ND3 apresenta um frameshift biologicamente real e documentado em aves (MINDELL; SORENSON; DIMCHEFF, 1998), levando à implementação de uma lógica de unificação (`join`) nos scripts de conversão para GenBank.

Cada correção foi implementada, testada e integrada antes do avanço para a próxima etapa, em conformidade com o princípio de retroalimentação rápida (*fail fast, fix immediately*) característico de metodologias ágeis (SOMMERVILLE, 2016).

#### 4.1.3 Aderência às Boas Práticas de Computação Científica

A adoção das práticas descritas acima situa este trabalho no contexto das boas práticas para computação científica sistematizadas por Wilson et al. (2014), que recomendam:

1. **Controle de versão** para todo código-fonte — implementado via Git com histórico completo de commits descritivos e branches para marcos significativos;
2. **Automação de tarefas repetitivas** — implementada via Nextflow, que orquestra todas as etapas sem intervenção manual;
3. **Testes com dados conhecidos** — implementados pela validação com espécie-controle (*D. bifasciatus*);
4. **Registro de dependências e ambientes** — implementado via Dockerfiles com versões explicitamente fixadas para todas as ferramentas;
5. **Documentação do propósito e das decisões de design** — implementada no repositório (README, guia de execução) e neste documento.

Essa abordagem metodológica não é apenas uma decisão de conveniência, mas uma necessidade reconhecida na literatura: estudos demonstram que a falta de práticas formais de engenharia de software em projetos científicos é uma das principais causas de irreprodutibilidade computacional (BAKER, 2016; PENG, 2011).

### 4.2 Visão Geral do Pipeline

O pipeline proposto neste trabalho foi concebido para realizar a montagem e anotação de genomas mitocondriais a partir de dados de sequenciamento de nova geração (NGS), integrando ferramentas bioinformáticas consolidadas em um ambiente de execução automatizado, reprodutível e portável. A arquitetura segue uma estrutura modular, em que cada etapa corresponde a uma fase distinta do processamento dos dados, sendo encapsulada em containers Docker para garantir consistência na execução (MERKEL, 2014) e integrada por meio do gerenciador de workflows Nextflow (DI TOMMASO et al., 2017).

De forma geral, o pipeline é composto por sete macroetapas principais, organizadas em uma cadeia de redução progressiva de dados que transforma dezenas de gigabytes de leituras brutas em um mitogenoma anotado de ~17 KB — viabilizando a execução em computadores pessoais sem necessidade de servidores. A primeira é a **análise piloto de qualidade (Pilot QC)**, um estágio opcional que analisa automaticamente uma subamostra dos dados para determinar o volume ideal de leituras a serem processadas, evitando downloads desnecessários em datasets grandes. A segunda etapa é a **aquisição e preparo dos dados**, responsável por obter datasets públicos em repositórios como o Sequence Read Archive (SRA) do NCBI (LEINONEN; SUGAWARA; SHUMWAY, 2011). Em seguida, realiza-se o **controle de qualidade e pré-processamento**, que envolve a avaliação das leituras por meio do FastQC (ANDREWS, 2010) e a remoção de adaptadores e bases de baixa qualidade utilizando o Trim Galore (BABRAHAM BIOINFORMATICS, 2019), que por sua vez utiliza o Cutadapt (MARTIN, 2011) internamente.

A quarta etapa corresponde à **montagem do mitogenoma**, realizada com o NOVOPlasty, ferramenta especializada que utiliza a estratégia seed-and-extend para reconstruir a molécula circular completa a partir de leituras curtas (DIERCKXSENS; MARDULYN; SMITS, 2017). A quinta etapa, de **anotação funcional**, identifica os genes do genoma montado, utilizando o MITOS2 (DONATH et al., 2019) para a predição de genes codificadores de proteínas, rRNAs e tRNAs, com geração de estruturas secundárias previstas via ViennaRNA (LORENZ et al., 2011). A sexta etapa consiste na **compilação dos entregáveis**, em que um script automatizado reúne todos os produtos da análise em uma pasta organizada, incluindo arquivos no formato GenBank Flat File para submissão ao NCBI. Por fim, a sétima etapa concentra-se na **automação e reprodutibilidade**, em que cada componente do pipeline é containerizado com Docker (MERKEL, 2014) e orquestrado via Nextflow (DI TOMMASO et al., 2017).

> **Figura 9** — Fluxograma do pipeline com as 7 etapas: Pilot QC, SRA Download, FastQC, Trim Galore, NOVOPlasty, MITOS2 e Compile Summary. Cores distinguem Fase 1 (montagem, azul) e Fase 2 (anotação, verde).


A visão geral do pipeline, representada em formato de fluxograma, ilustra o encadeamento das etapas e destacando a integração entre análise bioinformática e práticas de engenharia de software. Deve-se observar que a lógica condicional do pipeline permite que a etapa de anotação (MITOS2) seja executada apenas quando o banco de dados de referência estiver disponível, tornando a Fase 1 (montagem) independente da Fase 2 (anotação). Adicionalmente, o workflow aplica um **filtro de circularização**: apenas montagens cujo nome de arquivo inicia com `Circularized_assembly` são encaminhadas ao MITOS2 e à compilação de entregáveis. Quando o NOVOPlasty não obtém circularização (produzindo apenas contigs parciais), a anotação é automaticamente omitida, evitando análises sobre montagens incompletas que poderiam gerar resultados biologicamente incorretos. Essa decisão de design garante que os entregáveis finais representem exclusivamente genomas mitocondriais completos e circulares.

### 4.3 Fluxo de Execução do Pipeline

Esta seção descreve o fluxo completo de execução do pipeline, desde os pré-requisitos até a obtenção dos resultados finais, detalhando as decisões automáticas tomadas pelo sistema e o caminho percorrido pelos dados em cada etapa.

#### 4.3.1 Pré-requisitos e Preparação do Ambiente

Antes da primeira execução, três preparações são necessárias:

1. **Construção das imagens Docker** — o repositório inclui cinco Dockerfiles (SRA-Toolkit, FastQC, Trim Galore, NOVOPlasty e MITOS2), cada um com versões pinadas de todas as dependências. A construção é realizada por um script automatizado (`build_images.ps1` / equivalente Bash), que gera as imagens localmente sem necessidade de registro externo;
2. **Download do banco de dados MITOS2** — o banco RefSeq89m (~20 MB compactado) é obtido do Zenodo e extraído em `data/databases/`. Esse passo é necessário apenas uma vez, pois o banco é reutilizado entre execuções;
3. **Configuração do perfil da espécie** — cada espécie requer um arquivo de configuração em `conf/` definindo o acesso SRA (`sra_accession`), a semente (`seed`), o intervalo de tamanho esperado (`genome_range`), os k-mers a serem tentados (`novoplasty_kmers`), o código genético (`genetic_code`) e o nome científico (`organism`).

#### 4.3.2 Comando de Execução e Resiliência

A execução é iniciada por um único comando:

```
nextflow run main.nf -profile <espécie>
```

Onde `<espécie>` corresponde ao nome do perfil definido em `nextflow.config` (e.g., `a_leari`, `test`). Para execuções de longa duração — típicas em bioinformática, onde o download de dados do SRA pode levar dezenas de minutos e a conversão dos FASTQs pode exceder uma hora —, o pipeline inclui um script wrapper (`run_pipeline.sh`) que encapsula a execução com `nohup`, permitindo que o terminal seja fechado ou a conexão SSH seja interrompida sem afetar o pipeline. Os logs são automaticamente direcionados para `logs/<perfil>_<timestamp>.log`, possibilitando acompanhamento assíncrono via `tail -f`.

Para retomar uma execução interrompida, o Nextflow oferece o mecanismo `-resume`, que reutiliza os resultados de etapas já concluídas armazenados no diretório `work/`. O pipeline organiza esse diretório por espécie e timestamp (`work/<espécie>/<timestamp>/`), permitindo identificar e gerenciar execuções individuais. Para retomar uma execução específica, basta indicar o diretório correspondente:

```
nextflow run main.nf -profile <espécie> -resume -w work/<espécie>/<timestamp>
```

#### 4.3.3 Fluxo de Dados e Decisões Automáticas

A Figura a seguir ilustra o caminho dos dados através do pipeline, incluindo os pontos de decisão automática.

> **Figura 10** — Fluxo sequencial de dados no pipeline, mostrando os formatos e volumes em cada transição, desde o arquivo .sra (11 GB) até os 14 entregáveis finais (1,6 MB). Destacam-se a redução de 87 GB para 17 KB e os dois pontos de decisão automática.

Como evidenciado na Figura 10, o pipeline implementa uma estratégia de redução progressiva: o dado original de 87 GB (FASTQ bruto) é reduzido em mais de cinco milhões de vezes até a montagem final de 17 KB. Essa redução não é uniforme — o maior salto ocorre na etapa de montagem, quando o NOVOPlasty seleciona apenas leituras mitocondriais dentre milhões de leituras genômicas.

O pipeline incorpora dois pontos de decisão automática que não requerem intervenção do usuário:

**Decisão 1 — Pilot QC (volume de dados).** Quando o parâmetro `sra_max_reads` não é definido no perfil da espécie, o pipeline ativa automaticamente o módulo Pilot QC (seção 4.7), que baixa uma subamostra de 500.000 leituras, analisa sua qualidade e calcula o número ideal de leituras para a execução principal. Quando o parâmetro é explicitamente definido (e.g., via `--sra_max_reads 20000000`), o Pilot QC é pulado e o valor fornecido é utilizado diretamente. Essa lógica condicional é implementada no workflow principal (`main.nf`) por meio de operadores Nextflow de branching.

**Decisão 2 — Filtro de Circularização (qualidade da montagem).** Após a montagem pelo NOVOPlasty, apenas arquivos cujo nome corresponde ao padrão `Circularized_assembly_*` são encaminhados ao MITOS2. Montagens parciais (contigs não circularizados) são automaticamente descartados do fluxo de anotação, evitando que genomas incompletos recebam anotações potencialmente incorretas (e.g., genes truncados nas extremidades de contigs lineares). O pipeline reporta essas montagens parciais no log para que o pesquisador possa investigar manualmente a causa da não-circularização (cobertura insuficiente, semente inadequada, ou regiões repetitivas complexas).

> **Figura 11** — Diagrama de atividade mostrando os dois pontos de decisão automática do pipeline: (1) Pilot QC — determina o volume de dados; (2) Filtro de Circularização — condiciona a anotação à obtenção de genoma circularizado.

O diagrama da Figura 11 explicita que ambas as decisões são binárias e automáticas, dispensando intervenção do usuário. Na Decisão 1, o Pilot QC avalia a qualidade das leituras e estima a fração mitocondrial para determinar o volume ideal; na Decisão 2, o filtro de circularização verifica se a montagem foi circularizada — apenas montagens completas prosseguem para anotação.

#### 4.3.4 Extensibilidade para Novas Espécies

A adição de uma nova espécie ao pipeline requer apenas a criação de um arquivo de configuração e a obtenção da semente correspondente, sem modificação de código-fonte. Os passos são:

1. Identificar no NCBI SRA um dataset de sequenciamento de genoma total para a espécie-alvo;
2. Localizar o mitogenoma de uma espécie **congênere** — isto é, pertencente ao mesmo gênero taxonômico — (ou filogeneticamente próxima) no NCBI e extrair o gene cox1 como semente;
3. Criar um arquivo `conf/<especie>.config` definindo os parâmetros específicos (acesso SRA, semente, genome_range, k-mers, código genético e nome científico);
4. Registrar o novo perfil em `nextflow.config`;
5. Executar: `nextflow run main.nf -profile <especie>`.

Essa arquitetura extensível é demonstrada pelas sete sementes incluídas no repositório (`data/seeds/`), que cobrem espécies de aves, peixes, insetos e nematóides, atestando a generalidade do pipeline para diferentes grupos taxonômicos.


### 4.4 Arquitetura Modular e Scripts do Pipeline

A organização do código segue a convenção de projetos Nextflow DSL2, em que cada etapa do pipeline corresponde a um **módulo** independente (arquivo `.nf` no diretório `modules/`), e os scripts auxiliares residem no diretório `scripts/`. Essa separação permite que cada componente seja desenvolvido, testado e mantido de forma isolada. A Tabela 2 apresenta a relação completa dos módulos, scripts e suas responsabilidades.

> **[SUGESTÃO DE TABELA — Tabela 2]**: Inserir tabela com a arquitetura modular:

| Arquivo | Tipo | Função |
|---|---|---|
| `main.nf` | Workflow | Orquestração geral: define a ordem de execução, lógica condicional (Pilot QC, MITOS2) e passagem de dados entre módulos |
| `modules/sra_download.nf` | Módulo | Download do SRA via `prefetch` + conversão condicional (`fastq-dump -X` ou `fasterq-dump`) |
| `modules/fastqc.nf` | Módulo | Controle de qualidade das leituras brutas via FastQC |
| `modules/trim_galore.nf` | Módulo | Remoção de adaptadores e trimming de qualidade via Trim Galore |
| `modules/novoplasty.nf` | Módulo | Montagem iterativa com múltiplos k-mers e re-seeding automático |
| `modules/mitos2.nf` | Módulo | Anotação funcional via MITOS2 + geração de SVGs com RNAplot |
| `scripts/pilot_qc.sh` | Script Bash | Análise piloto: calcula Q30%, taxa de adaptador, fração mitocondrial e recomenda `max_reads` |
| `scripts/compile_summary.py` | Script Python | Compila 14 categorias de entregáveis a partir das saídas dos módulos |
| `scripts/gff2genbank.py` | Script Python | Converte GFF3 → feature table (`.tbl`) + FASTA formatado (`.fsa`) para GenBank, com tratamento do frameshift ND3 |
| `scripts/generate_genbank.py` | Script Python | Gera GenBank Flat File (`.gbk`) e mapa circular (SVG/PDF) via Biopython |
| `nextflow.config` | Configuração | Parâmetros globais, perfis por espécie, configuração de containers Docker |
| `conf/*.config` | Configuração | Perfis específicos por espécie (acesso SRA, semente, k-mers, banco MITOS2) |

Essa arquitetura modular permite que novos módulos sejam adicionados (e.g., análise filogenética, BLAST) sem alterar os existentes — basta criar o arquivo `.nf`, o Dockerfile correspondente e registrar a chamada no `main.nf`.

### 4.5 Ferramentas e Tecnologias Utilizadas

A implementação de um pipeline bioinformático demanda não apenas a escolha de algoritmos adequados, mas também a definição criteriosa de ferramentas que sejam estáveis, bem documentadas e amplamente utilizadas pela comunidade científica. Nesse contexto, a seleção de softwares para este trabalho foi orientada por três princípios fundamentais: (i) funcionalidade biológica, garantindo que cada ferramenta seja apropriada para a etapa que executa, desde a aquisição dos dados até a anotação final do genoma; (ii) confiabilidade e validação prévia, privilegiando programas consolidados na literatura e em uso corrente em pesquisas genômicas; e (iii) compatibilidade com práticas modernas de engenharia de software, assegurando que cada aplicação possa ser encapsulada em containers e integrada a sistemas de gerenciamento de workflows.

Um aspecto central da engenharia deste pipeline é o **pinamento de versões**. Todas as ferramentas foram encapsuladas em containers Docker com versões explicitamente fixadas, tanto para as imagens base quanto para os softwares instalados. Essa decisão, embora impeça o uso automático de atualizações, garante que o pipeline produza resultados idênticos independentemente do momento ou do ambiente em que for executado, eliminando uma das principais fontes de irreprodutibilidade em bioinformática. A Tabela 1 apresenta as versões de cada componente utilizado.

> **[SUGESTÃO DE TABELA — Tabela 1]**: Inserir tabela com as versões pinadas:

| Componente | Versão | Imagem base |
|---|---|---|
| SRA-Toolkit | 3.0.10 | ubuntu:22.04 |
| FastQC | 0.12.1 | ubuntu:22.04 |
| Trim Galore | 0.6.10 | ubuntu:22.04 |
| Cutadapt | 4.6 | (dentro do Trim Galore) |
| NOVOPlasty | 4.3.1 | ubuntu:22.04 |
| MITOS2 | 2.1.9 | continuumio/miniconda3:24.1.2-0 |
| ViennaRNA | (via conda) | (dentro do MITOS2) |
| Nextflow | DSL2 (25.x) | - |

Além de atender a esses critérios, as ferramentas escolhidas refletem o estado da arte em bioinformática aplicada ao estudo de organelas, cobrindo desde utilitários básicos, como o SRA-Toolkit para acesso a bancos de dados, até plataformas mais complexas, como o NOVOPlasty para montagem especializada de mitogenomas. Da mesma forma, tecnologias transversais como Docker e Nextflow foram incorporadas não pelo papel biológico que desempenham, mas pela sua relevância em assegurar portabilidade, escalabilidade e reprodutibilidade, elementos indispensáveis para a robustez metodológica do pipeline.

Assim, a seguir neste capítulo é descrito de forma integrada as principais ferramentas e tecnologias que compõem o pipeline de montagem de genomas mitocondriais.

#### 4.5.1 SRA-Toolkit

O SRA-Toolkit é o pacote oficial de utilitários do NCBI para manipulação de dados depositados no Sequence Read Archive (SRA), o principal repositório público de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). No pipeline desenvolvido, dois utilitários do SRA-Toolkit são empregados em sequência: o `prefetch`, que realiza o download eficiente do arquivo `.sra` compactado, e um dos conversores para formato FASTQ paired-end — `fasterq-dump` ou `fastq-dump` — selecionado automaticamente conforme as condições de execução.

Uma decisão de implementação relevante diz respeito à estratégia de redução de dados. Em datasets de grande volume — como o utilizado para a arara-azul-de-lear, com 118,5 milhões de pares de leituras (SRR28399504) —, o processamento integral seria desnecessário e computacionalmente oneroso, dado que o genoma mitocondrial de aves tem cerca de 17 kb (BOORE, 1999) e apenas uma fração diminuta das leituras (~0,17%) corresponde ao mtDNA. O pipeline implementa duas estratégias condicionais de redução:

- **Quando `sra_max_reads` está definido**: utiliza o `fastq-dump` com a flag `-X`, que extrai diretamente apenas os primeiros *N* spots do arquivo `.sra` sem necessidade de conversão integral. Embora o `fastq-dump` seja single-thread e mais lento que o `fasterq-dump`, a redução drástica do volume processado compensa amplamente — por exemplo, extrair 25M spots de um dataset de 145,8M processa apenas 17% dos dados, com uso de disco proporcional.
- **Quando `sra_max_reads` não está definido**: utiliza o `fasterq-dump`, que é multi-threaded e até 10 vezes mais rápido que o `fastq-dump` em datasets grandes (NCBI, 2023), convertendo o arquivo integralmente para FASTQ.

Essa abordagem condicional foi adotada após a execução com o dataset metagenômico de *A. hyacinthinus* (145,8M spots) revelar que a estratégia original — conversão integral pelo `fasterq-dump` seguida de truncamento por `head` — era inviável para grandes volumes, por gerar ~124 GB de arquivos FASTQ intermediários (Seção 5.9).

#### 4.5.2 FastQC

O FastQC é uma das ferramentas mais utilizadas na etapa de controle de qualidade de dados de sequenciamento de nova geração (NGS). Desenvolvido pelo Babraham Institute, tornou-se referência para avaliação inicial de leituras em formato FASTQ, sendo amplamente empregado em pipelines de bioinformática devido à sua rapidez, robustez e facilidade de interpretação (ANDREWS, 2010).

Do ponto de vista técnico, o FastQC processa os arquivos de entrada em formato FASTQ e gera como saída relatórios em HTML interativo e arquivos compactados (.zip), contendo gráficos e estatísticas sobre diferentes aspectos das sequências. Entre os módulos mais relevantes estão:

- Per base sequence quality, que avalia a distribuição da qualidade das bases em cada posição do read;
- Per sequence quality scores, que identifica leituras com baixa qualidade geral;
- Per base GC content, que compara a distribuição de GC observada com a esperada;
- Overrepresented sequences, que detecta possíveis contaminantes ou adaptadores;
- Sequence duplication levels, que informa o grau de redundância entre leituras;
- Adapter content, que indica a presença de adaptadores de sequenciamento não removidos.

Entre as principais vantagens do FastQC estão: (i) a rapidez de execução mesmo em grandes conjuntos de dados, (ii) a geração de relatórios de fácil interpretação, e (iii) a capacidade de identificar múltiplos problemas em uma única execução. No entanto, o programa possui limitações: por ser uma análise pré-alinhamento, não oferece informações sobre a qualidade do mapeamento ou sobre a distribuição de leituras no genoma de referência. Além disso, os resultados devem ser interpretados com cautela, já que pequenas variações em métricas, como o conteúdo GC, podem refletir características biológicas legítimas e não necessariamente artefatos técnicos.

A escolha pelo FastQC neste trabalho se justifica por sua ampla aceitação na comunidade científica como ferramenta de triagem inicial de qualidade. Sua aplicação permite identificar rapidamente problemas técnicos, orientar a necessidade de trimming ou filtragem, e garantir que apenas leituras de qualidade adequada sejam utilizadas nas etapas de montagem do mitogenoma. Dessa forma, o FastQC atua como um ponto de controle essencial para assegurar a confiabilidade das análises subsequentes.

#### 4.5.3 Cutadapt e Trim Galore

Após a avaliação inicial da qualidade das leituras, pode ser necessário realizar etapas de filtragem para remover **adaptadores** residuais e **bases de baixa qualidade**. Para compreender essas operações, é necessário esclarecer dois conceitos fundamentais.

**Adaptadores** são sequências sintéticas curtas (tipicamente 30–60 nucleótidos) que são quimicamente ligadas às extremidades dos fragmentos de DNA durante a preparação da biblioteca de sequenciamento. Funcionam como “ancoradouros” que permitem ao sequenciador reconhecer e processar cada fragmento. Após o sequenciamento, essas sequências artificiais não fazem parte do genoma do organismo e precisam ser removidas antes da montagem, pois sua presença geraria sobreposições falsas entre leituras, comprometendo a reconstrução do genoma.

**Bases de baixa qualidade** referem-se a nucleótidos cuja identificação pelo sequenciador apresenta baixa confiança. A qualidade de cada base é expressa pela **pontuação Phred** ($Q$), uma escala logarítmica em que $Q = -10 \log_{10}(P)$, sendo $P$ a probabilidade de erro (EWING et al., 1998). Assim, $Q = 20$ indica 1% de probabilidade de erro (99% de acerto), enquanto $Q = 30$ indica 0,1% de erro (99,9% de acerto). Bases com $Q < 20$ são geralmente consideradas de baixa qualidade e removidas no processo de **trimming** (corte) — a operação de remoção systemática de adaptadores e bases de baixa qualidade das extremidades das leituras.

Nesse contexto, duas ferramentas complementares são utilizadas no pipeline: Cutadapt e Trim Galore.

O Cutadapt é um programa desenvolvido por Martin (2011) especificamente para identificar e remover sequências de adaptadores em dados de alto rendimento. As leituras geradas por tecnologias como Illumina frequentemente contêm fragmentos de adaptadores ligados durante a preparação da biblioteca, que, se não forem removidos, podem comprometer a montagem ou o mapeamento subsequente. O Cutadapt utiliza algoritmos de busca eficiente para detectar esses fragmentos nas extremidades das leituras e removê-los, além de permitir o corte de bases de baixa qualidade de acordo com pontuações Phred. Como entrada, o programa recebe arquivos **FASTQ** — formato de texto que armazena tanto a sequência de cada leitura quanto a qualidade Phred de cada base individual, representada por caracteres ASCII —, e como saída gera novos FASTQ já filtrados, prontos para análises subsequentes.

O Trim Galore é um wrapper que automatiza a execução do Cutadapt em conjunto com o FastQC, simplificando o processo de trimming e tornando-o mais acessível ao usuário (BABRAHAM BIOINFORMATICS, 2019). Com ele, é possível remover adaptadores mesmo quando a sequência exata não é previamente conhecida, o que é particularmente útil em experimentos com bibliotecas de origem diversa. Além disso, o Trim Galore aplica rotinas de filtragem por qualidade e integra relatórios do FastQC para avaliar o impacto do corte.

Entre as principais vantagens do Cutadapt e do Trim Galore estão: (i) a remoção precisa de adaptadores, prevenindo erros de montagem, (ii) a flexibilidade de configuração de parâmetros, e (iii) a integração automatizada que reduz a intervenção manual. Suas limitações incluem a necessidade de escolher cuidadosamente os limiares de qualidade, já que cortes excessivos podem resultar em perda de informação biológica relevante, enquanto cortes insuficientes podem deixar contaminantes residuais.

A escolha por essas ferramentas neste trabalho se justifica pelo seu amplo uso na comunidade científica e pela confiabilidade de seus resultados. A integração de Cutadapt e Trim Galore garante que as leituras que avançam para a montagem apresentem qualidade adequada, reduzindo o risco de artefatos e aumentando a robustez dos resultados obtidos pelo pipeline.

#### 4.5.4 NOVOPlasty

O NOVOPlasty é um montador especializado desenvolvido para reconstrução de genomas de organelas, como mitocôndrias e cloroplastos, a partir de dados de sequenciamento de genoma total. Publicado por Dierckxsens, Mardulyn e Smits (2017), tornou-se uma das ferramentas mais amplamente utilizadas para esse fim devido à sua eficiência e facilidade de uso.

Do ponto de vista técnico, o NOVOPlasty baseia-se na estratégia seed-and-extend. O processo de montagem é iniciado a partir de uma sequência fornecida pelo usuário — denominada seed — que pode ser um fragmento do genoma da organela em estudo, uma sequência de uma espécie próxima ou até mesmo uma sequência completa de organela de outro organismo. A partir dessa seed, o software realiza uma extensão bidirecional, identificando leituras sobrepostas e construindo gradualmente a sequência circular completa do genoma. Para acelerar o processo, o programa utiliza tabelas de hash, que permitem buscas rápidas por sobreposições entre as leituras.

O NOVOPlasty recebe como entrada arquivos FASTQ contendo leituras pareadas, bem como um arquivo FASTA com a sequência seed. A saída consiste em um arquivo FASTA representando o genoma circularizado, acompanhado de arquivos de log e, em alguns casos, múltiplas opções de montagem (quando diferentes caminhos são possíveis no grafo de montagem).

Entre as principais vantagens do NOVOPlasty destacam-se: (i) sua alta eficiência na recuperação de genomas completos de organelas mesmo a partir de dados de genoma total, (ii) a baixa exigência de parâmetros complexos, o que o torna acessível a usuários com diferentes níveis de experiência, e (iii) a capacidade de lidar com diferentes tipos de sementes, permitindo flexibilidade em estudos de espécies com genomas pouco conhecidos. No entanto, como limitações, pode-se destacar: (i) a dependência da qualidade e da escolha adequada da seed, que influencia diretamente o sucesso da montagem, e (ii) o risco de circularizações incorretas em casos de regiões altamente repetitivas ou coberturas desbalanceadas.

No pipeline implementado, o módulo NOVOPlasty incorpora uma lógica iterativa com re-seeding automático: caso a primeira execução não resulte em circularização, o maior contig obtido é automaticamente utilizado como nova semente em uma segunda tentativa, e múltiplos valores de k-mer são experimentados em sequência. Essa estratégia aumenta a taxa de sucesso da montagem sem requerer intervenção manual do usuário.

#### 4.5.5 MITOS2

O MITOS2 é uma ferramenta desenvolvida para a anotação funcional de genomas mitocondriais em metazoários. Apresentada por Donath et al. (2019), constitui a evolução da primeira versão do MITOS, oferecendo maior precisão na predição de genes, integração com bases de dados atualizadas e maior flexibilidade de execução.

Tecnicamente, o MITOS2 combina métodos de homologia (comparação com sequências conhecidas) e algoritmos de predição *ab initio* (a partir da própria sequência, sem referência externa) para identificar as principais classes de genes presentes em genomas mitocondriais. O sistema utiliza alinhamentos com modelos de substituição específicos para proteínas codificadas por mitocôndrias, além de perfis de RNA para identificar rRNAs e tRNAs por meio de modelos de covariância (Infernal). As entradas consistem em um arquivo FASTA contendo a sequência do genoma mitocondrial a ser anotado, e as saídas incluem: (i) arquivos tabulares com as coordenadas dos genes anotados, (ii) representações gráficas do genoma mostrando a organização gênica, e (iii) arquivos de anotação em formatos padrão, como GFF3 e BED (formato tabular simplificado com coordenadas genômicas). É relevante notar que, conforme demonstrado por Donath et al. (2019), as anotações do RefSeq contêm erros — em particular, atribuições implausíveis de posições de códons de início e parada, além de uma fração substancial de anotações incompletas que identificam apenas fragmentos de genes codificadores de proteínas. O MITOS2 demonstrou que uma anotação totalmente automática com acurácia muito alta é possível, representando uma melhoria drástica em relação à determinação manual de fronteiras gênicas.

A versão empregada neste trabalho (v2.1.9) foi instalada via Conda em um container Docker baseado em Miniconda 24.1.2-0, utilizando o banco de dados RefSeq89m obtido do Zenodo (DONATH et al., 2019). Essa abordagem garante que a ferramenta, normalmente executada via interface web, funcione de forma autônoma em linha de comando e completamente integrada ao workflow Nextflow.

Como extensão funcional, o pipeline integra automaticamente a geração de representações visuais da estrutura secundária dos tRNAs e rRNAs preditos. Para isso, o módulo MITOS2 extrai as sequências e estruturas dot-bracket dos arquivos de saída e as processa com o RNAplot, componente do pacote ViennaRNA (LORENZ et al., 2011), gerando arquivos SVG individuais para cada RNA não codificante. Essa funcionalidade, originalmente disponível apenas na versão web do MITOS2 (UseGalaxy), foi reproduzida localmente no pipeline, garantindo que os diagramas de estrutura secundária sejam gerados sem dependência de serviços externos.

Entre as principais vantagens do MITOS2 destacam-se: (i) a automação completa da etapa de anotação, reduzindo a necessidade de curadoria manual extensiva; (ii) a compatibilidade com diferentes grupos taxonômicos de metazoários; e (iii) a produção de relatórios gráficos de fácil interpretação. Em estudos recentes de genômica mitocondrial comparativa, a anotação do MITOS2 é frequentemente complementada pela verificação de tRNAs com tRNAscan-SE, ferramenta especializada que utiliza modelos de covariância otimizados para a predição de RNAs transportadores (ZHAN et al., 2024; KUNDU et al., 2024; SHI et al., 2024). Essa combinação MITOS2 + tRNAscan-SE tem se consolidado como padrão metodológico na área. No presente pipeline, o MITOS2 é utilizado como ferramenta única de anotação, dado que seus próprios modelos de covariância (via Infernal) já incorporam predição de tRNAs com alta sensibilidade; a integração opcional de tRNAscan-SE como módulo de verificação adicional constitui uma perspectiva futura. Por outro lado, como limitações, pode apresentar divergências em regiões com estruturas gênicas incomuns — como o caso do gene ND3 em aves, que apresenta um **frameshift insertion** (inserção de um nucleotídeo que desloca o quadro de leitura do código genético, alterando a sequência de aminoácidos traduzida daquele ponto em diante), documentado por Mindell et al. (1998), levando o MITOS2 a reportar o gene como dois **ORFs** (*Open Reading Frames* — quadros abertos de leitura, trechos de DNA que potencialmente codificam proteínas) distintos (nad3_0 e nad3_1) em vez de um único CDS (*Coding Sequence* — sequência codificante) com `join()`. Esse tipo de situação requer validação complementar e tratamento específico na conversão para formatos de submissão ao GenBank.

#### 4.5.6 Docker

No pipeline proposto, cada ferramenta foi encapsulada em um container Docker, assegurando que todas as dependências fossem executadas em ambientes controlados e consistentes (MERKEL, 2014). O pipeline utiliza cinco imagens Docker distintas, cada uma correspondendo a um estágio do processamento:

| Imagem | Base | Ferramenta | Dockerfile |
|---|---|---|---|
| `mitogenome-pipeline/sra-tools:1.0` | ubuntu:22.04 | SRA-Toolkit 3.0.10 | `docker/sra-tools/Dockerfile` |
| `mitogenome-pipeline/fastqc:1.0` | ubuntu:22.04 | FastQC 0.12.1 | `docker/fastqc/Dockerfile` |
| `mitogenome-pipeline/trim-galore:1.0` | ubuntu:22.04 | Trim Galore 0.6.10 + Cutadapt 4.6 | `docker/trim-galore/Dockerfile` |
| `mitogenome-pipeline/novoplasty:1.0` | ubuntu:22.04 | NOVOPlasty 4.3.1 | `docker/novoplasty/Dockerfile` |
| `mitogenome-pipeline/mitos2:1.0` | miniconda3:24.1.2-0 | MITOS2 2.1.9 + ViennaRNA | `docker/mitos2/Dockerfile` |

Essa abordagem eliminou a necessidade de configurações manuais complexas e garantiu que softwares como FastQC, Cutadapt, NOVOPlasty e MITOS2 operassem de forma uniforme em diferentes sistemas operacionais. Com isso, a montagem e a anotação de mitogenomas tornam-se replicáveis em qualquer ambiente computacional que disponha do Docker instalado, ampliando a portabilidade do pipeline.

#### 4.5.7 Nextflow

A orquestração do pipeline foi implementada no Nextflow DSL2, que estruturou os processos em módulos independentes (um arquivo `.nf` por etapa), conectando automaticamente suas entradas e saídas (DI TOMMASO et al., 2017). Essa organização permitiu automatizar a lógica iterativa da montagem, como a possibilidade de reaproveitar contigs parciais do NOVOPlasty como novas seeds em execuções subsequentes. Além disso, o Nextflow possibilitou a paralelização de tarefas, a retomada de execuções interrompidas (via flag `-resume`) e a integração direta com containers Docker. Dessa forma, o pipeline desenvolvido combina simplicidade de uso com escalabilidade e reprodutibilidade, adequando-se tanto a computadores locais quanto a servidores de alto desempenho.

O sistema de perfis do Nextflow foi utilizado para parametrizar execuções distintas sem alterar o código-fonte. Cada espécie estudada possui um arquivo de configuração próprio (e.g., `conf/a_leari.config`) contendo os parâmetros específicos — acesso SRA, semente, tamanho esperado do genoma, código genético — que são carregados automaticamente ao invocar o perfil correspondente (`nextflow run main.nf -profile a_leari`).

#### 4.5.8 GitHub

O GitHub é uma plataforma de hospedagem de código baseada no sistema de controle de versão Git, amplamente utilizada para o desenvolvimento colaborativo de software e projetos científicos. Para pipelines de bioinformática, a disponibilização em um repositório GitHub garante não apenas transparência, mas também acessibilidade e manutenção a longo prazo.

Neste trabalho, o GitHub foi utilizado como repositório público para disponibilizar o pipeline desenvolvido, assegurando sua acessibilidade, versionamento e documentação adequada. Essa prática se alinha às diretrizes de ciência aberta e reforça o compromisso com a reprodutibilidade e a transparência metodológica. O repositório está disponível em https://github.com/matheus-sobreira/mitogenome-pipeline e inclui: o código completo do workflow Nextflow, todos os Dockerfiles, scripts auxiliares, configurações de perfil por espécie, dados de semente (sequências cox1), documentação e guia de execução.

Adicionalmente, o sistema de branches do Git foi utilizado para preservar estados estáveis do pipeline: a branch `fase1-montagem` contém o estado validado da Fase 1 (apenas montagem, sem anotação), enquanto a branch `main` contém o pipeline completo com todas as etapas integradas.

### 4.6 Aquisição e Preparo dos Dados

A primeira etapa prática do pipeline consiste na aquisição e organização dos dados de sequenciamento que servirão de insumo para a montagem dos genomas mitocondriais. Para este trabalho, foram utilizados datasets públicos provenientes do Sequence Read Archive (SRA), repositório mantido pelo NCBI e reconhecido como um dos principais bancos internacionais de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). A escolha por dados públicos justifica-se pela sua ampla disponibilidade, pela diversidade de organismos representados e pela possibilidade de reproduzir e validar resultados obtidos em trabalhos prévios.

O acesso aos dados é realizado por meio do SRA-Toolkit v3.0.10 (encapsulado em container Docker), por meio de dois comandos sequenciais: `prefetch`, responsável pelo download do arquivo `.sra` compactado, e `fastq-dump -X` ou `fasterq-dump`, que convertem o arquivo para o formato FASTQ paired-end — o primeiro quando há limite de leituras definido, o segundo quando se deseja processar o dataset integralmente. Implementou-se também a estratégia de retries automáticos (até 2 tentativas) para lidar com instabilidades de rede, comuns em downloads de datasets grandes do NCBI.

Como estudo de caso principal, foi utilizado o dataset SRR28399504 referente ao genoma da arara-azul-de-lear (*Anodorhynchus leari*), gerado por sequenciamento paired-end na plataforma Illumina HiSeq X Ten (IRIDIAN GENOMES, 2024), com leituras de 150 pares de bases e um total de 118,5 milhões de pares. Para validação do pipeline, foi utilizado o dataset SRR36182901 referente a *Diploprion bifasciatus* (peixe-sabão barrado), gerado em plataforma NovaSeq X Plus com leituras de 151 bp.

### 4.7 Análise Piloto de Qualidade (Pilot QC)

Uma inovação incorporada ao pipeline é o módulo de análise piloto de qualidade, que soluciona um problema prático recorrente em montagens de organelas: a determinação do volume adequado de dados a serem processados. Em datasets de sequenciamento de genoma total, a proporção de leituras mitocondriais é tipicamente inferior a 1% (DIERCKXSENS; MARDULYN; SMITS, 2017), o que significa que processar integralmente um dataset de 100 GB resultaria em centenas de gigabytes de dados nucleares irrelevantes para a montagem do mitogenoma. Além da ineficiência computacional, esse volume torna a análise impraticável em computadores pessoais: um notebook com 256 GB de SSD e 16 GB de RAM não comportaria os arquivos temporários gerados pela conversão e processamento de 87 GB de FASTQs brutos. O Pilot QC resolve esse problema logo na primeira etapa, reduzindo drasticamente o volume de dados antes que qualquer processamento pesado ocorra.

O módulo Pilot QC opera da seguinte forma: quando o parâmetro `sra_max_reads` não é definido pelo usuário, o pipeline baixa automaticamente uma subamostra de 500.000 leituras (aproximadamente 500 MB) utilizando o `fastq-dump -X`, que limita a extração na origem sem necessidade de transferir o dataset completo. Essa mesma estratégia de extração limitada por `fastq-dump -X` é empregada pelo módulo `SRA_DOWNLOAD` quando `sra_max_reads` está definido, garantindo consistência e eficiência de armazenamento em todo o pipeline. Essa subamostra é então submetida ao script `scripts/pilot_qc.sh`, um analisador desenvolvido inteiramente em Bash que opera diretamente sobre os arquivos FASTQ sem dependências externas além de `awk` e `grep`. O script extrai três métricas-chave:

1. **Proporção de bases com qualidade ≥ Q30**, indicando a qualidade geral das leituras;
2. **Proporção de leituras com adaptador Illumina TruSeq** residual, estimando a perda efetiva no trimming;
3. **Fração mitocondrial estimada**, calculada por correspondência de k-mers da sequência semente contra as leituras piloto.

A partir dessas métricas, o script calcula o número total de leituras necessário para atingir uma cobertura-alvo de 500× no mitogenoma, aplicando fatores de correção por qualidade e por presença de adaptador, com uma margem de segurança de 1,5×. O resultado é limitado ao intervalo de 5 a 50 milhões de leituras. Esse valor é então passado automaticamente ao módulo de download principal, que extrai diretamente apenas o número calculado de reads via `fastq-dump -X`.

A fórmula de cálculo é:

$$R_{mito} = \frac{C_{alvo} \times G}{L \times 2}$$

$$R_{total} = \frac{R_{mito} \times F_{Q30} \times F_{adapter}}{f_{mito}} \times 1.5$$

Onde $C_{alvo}$ é a cobertura desejada, $G$ é o tamanho médio esperado do genoma, $L$ é o comprimento médio das leituras, $f_{mito}$ é a fração mitocondrial estimada, e $F_{Q30}$ e $F_{adapter}$ são fatores de correção baseados na qualidade e presença de adaptadores.

Caso o usuário opte por definir manualmente o `sra_max_reads` via parâmetro de linha de comando, a etapa Pilot QC é automaticamente pulada, preservando a flexibilidade para cenários em que o pesquisador já possui conhecimento prévio sobre seus dados.

> **Figura 12** — Fluxograma de decisão do Pilot QC: análise de Q30%, taxa de adaptadores e fração mitocondrial estimada para cálculo automático do número ideal de leituras, com caps entre 5M e 25M.

A Figura 12 evidencia que o módulo opera exclusivamente sobre uma subamostra de 500K leituras, evitando o download completo do dataset antes de determinar seus parâmetros de qualidade.

### 4.8 Controle de Qualidade e Pré-processamento

Após a aquisição e organização dos dados de sequenciamento, a etapa seguinte do pipeline consiste na avaliação da qualidade das leituras e na aplicação de procedimentos de filtragem. Essa fase é essencial para assegurar que apenas dados confiáveis avancem para a montagem do genoma mitocondrial, evitando que artefatos técnicos comprometam a acurácia do resultado final.

O primeiro passo é a execução do FastQC v0.12.1, que gera relatórios interativos em HTML contendo métricas fundamentais para a avaliação da qualidade das leituras. Em seguida, é aplicado o Trim Galore v0.6.10 (que integra o Cutadapt v4.6), com parâmetros de qualidade Phred mínimo de 20 e comprimento mínimo de leitura de 50 bp para remoção de adaptadores e bases de baixa qualidade. No módulo implementado, o número de threads de CPU foi deliberadamente limitado a 2 para evitar ultrapassar a alocação de CPU no container, dado que o Trim Galore gera internamente threads adicionais para compressão (pigz) e para o próprio Cutadapt.

Além de aumentar a confiabilidade biológica das leituras, o pré-processamento contribui significativamente para a redução do volume de dados. Na execução para a arara-azul-de-lear, por exemplo, o Trim Galore detectou adaptadores Illumina TruSeq em 81,4% das leituras, e após o tricming, 71,9% das bases foram mantidas — uma redução de volume de ~28%.

Cabe observar que o pipeline não executa uma segunda análise FastQC após a trimagem. Embora essa prática seja adotada por alguns protocolos de QC genômico, ela foi deliberadamente omitida neste contexto por três razões: (i) o Trim Galore já produz um relatório detalhado com estatísticas de trimagem; (ii) a validação definitiva da qualidade dos dados ocorre *a posteriori*, pela circularização bem-sucedida do mitogenoma e pela anotação íntegra dos 37 genes — evidências funcionais mais robustas que métricas de QC intermediárias; e (iii) a omissão desta etapa reduz o tempo de execução total sem comprometer a confiabilidade dos resultados.

### 4.9 Montagem do Mitogenoma

A montagem do genoma mitocondrial constitui a etapa central do pipeline, pois é nela que se obtém a sequência circular completa que serve de base para a anotação funcional e análises comparativas. Para essa finalidade é utilizado o NOVOPlasty v4.3.1, um dos montadores mais consolidados para organelas (DIERCKXSENS; MARDULYN; SMITS, 2017).

No pipeline implementado, o módulo NOVOPlasty incorpora diversas melhorias em relação ao uso padrão da ferramenta:

1. **Múltiplos k-mers**: o parâmetro `novoplasty_kmers` permite definir uma lista de valores de k-mer a serem tentados em sequência (e.g., `'39,33'`). O pipeline experimenta cada um até obter circularização;
2. **Re-seeding automático**: para cada k-mer, caso a primeira execução não circularize, o maior contig obtido é automaticamente extraído e utilizado como nova semente em até 5 iterações (configurável);
3. **Liberação de espaço**: após a montagem, os arquivos FASTQ trimados são automaticamente deletados para economizar espaço em disco, comportamento crucial para execuções com dados volumosos.

A seleção dos valores de k-mer segue as recomendações do próprio NOVOPlasty (DIERCKXSENS; MARDULYN; SMITS, 2017), que aceita valores ímpares no intervalo de 21 a 39. Valores maiores de k-mer (como 39) favorecem a **especificidade** da montagem — cada k-mer é mais longo e, portanto, menos ambíguo, reduzindo o risco de quimeras e extensões incorretas. Por outro lado, valores menores (como 33) aumentam a **sensibilidade**, permitindo extensões em regiões de menor cobertura ou maior divergência entre a semente e o genoma-alvo. Assim, a estratégia adotada neste pipeline inicia pelo valor máximo (39), que produz montagens mais confiáveis quando a cobertura e a qualidade dos dados são adequadas, e recorre ao valor menor (33) apenas se a circularização não for alcançada na primeira tentativa. Para a arara-azul-de-lear, a circularização ocorreu na primeira tentativa com k-mer 39, indicando que a cobertura de 176× e a qualidade do HiSeq X Ten foram suficientes para o valor mais restritivo. Já para *D. bifasciatus*, a circularização ocorreu com k-mer 33, possivelmente refletindo diferenças na cobertura ou na composição nucleotídica do dataset.

> **Figura 13** — Diagrama da lógica iterativa de k-mers: tentativa com k=39 (alta especificidade) seguida de *fallback* para k=33 (maior sensibilidade) quando a circularização não é alcançada.

A Figura 13 mostra que o pipeline tenta inicialmente o k-mer mais alto (k=39), que oferece maior especificidade na resolução de regiões repetitivas. Caso a circularização não seja obtida, o pipeline reduz automaticamente para k=33.

**Obtenção da semente (seed).** O NOVOPlasty requer uma sequência inicial conhecida para ancorar a extensão. A estratégia adotada neste trabalho consiste em: (i) identificar no NCBI uma espécie congênere ou filogeneticamente próxima que possua mitogenoma completo depositado; (ii) extrair o gene cox1 dessa referência, por ser o marcador mitocondrial mais conservado e amplamente disponível; e (iii) salvar a região em formato FASTA. Para a arara-azul-de-lear, foi utilizado o gene cox1 de *A. hyacinthinus* (NC_082165.1, posições 5359–6906, 1.548 bp), espécie congênere cuja proximidade filogenética garante homologia suficiente para o seed-and-extend funcionar eficazmente. Para *D. bifasciatus*, foi utilizado o próprio cox1 da referência publicada (PZ143763.1, 1.560 bp). O repositório do pipeline inclui um guia detalhado para obtenção de sementes (`data/seeds/COMO_OBTER_SEMENTE.md`) e exemplos para sete espécies de diferentes grupos taxonômicos.

**Seleção do intervalo de tamanho (`genome_range`).** O NOVOPlasty utiliza um intervalo esperado de tamanho do genoma para descartar extensões espúrias. A determinação desse intervalo segue a mesma lógica da semente: busca-se no NCBI o tamanho do mitogenoma de referência da espécie congênere mais próxima e aplica-se uma margem de ±1.500–2.000 bp para acomodar variação interespecífica. Para a arara-azul-de-lear, a referência *A. hyacinthinus* (16.999 bp) resultou no intervalo 15.500–18.500 bp; para *D. bifasciatus*, a referência PZ143763.1 (16.805 bp) resultou em 15.000–18.500 bp.

**Parâmetro `insert_size`.** O tamanho do inserto da biblioteca é parametrizado como 300 bp (valor padrão conservador para preparações Illumina). O NOVOPlasty utiliza esse valor como estimativa inicial, mas refina-o automaticamente durante a montagem — no caso da arara-azul-de-lear, o inserto real medido foi de 213 bp, sem impacto na qualidade da circularização.

### 4.10 Anotação Funcional

Concluída a etapa de montagem, procede-se à anotação funcional do genoma mitocondrial reconstruído, utilizando o MITOS2 v2.1.9 com o banco de dados RefSeq89m. A anotação é executada automaticamente pelo pipeline quando o parâmetro `mitos2_db` está configurado, sendo condicional para permitir que a Fase 1 (montagem) opere independentemente.

O módulo MITOS2 recebe como entrada o arquivo FASTA da montagem circularizada e produz como saída: arquivo GFF3 com coordenadas dos genes, tabelas BED, sequências proteicas preditas (FAA), mapa linear do genoma (PNG) e plots de qualidade (PDF). Adicionalmente, o pipeline gera automaticamente diagramas SVG da estrutura secundária prevista para cada tRNA e rRNA identificado, utilizando o programa RNAplot do pacote ViennaRNA.

O código genético é parametrizado (`genetic_code = 2` para vertebrados mitocondriais, `genetic_code = 5` para invertebrados), permitindo que o mesmo pipeline seja utilizado para diferentes grupos taxonômicos sem modificação de código.

### 4.11 Compilação dos Entregáveis

Uma funcionalidade diferenciada do pipeline é a compilação automática de todos os produtos da análise em uma pasta organizada de entregáveis (`summary/`), pronta para apresentação acadêmica ou submissão a bancos de dados. Essa etapa é executada pelo módulo `COMPILE_SUMMARY`, que integra três scripts Python:

1. **`gff2genbank.py`** — Converte a anotação GFF3 do MITOS2 para o formato feature table (`.tbl`) do GenBank, com tratamento especial para o frameshift do ND3 em aves (uso de `join()` e `/exception=ribosomal slippage`);
2. **`generate_genbank.py`** — Gera o GenBank Flat File (`.gbk`) e o mapa circular do genoma (SVG e PDF) via Biopython GenomeDiagram (COCK et al., 2009);
3. **`compile_summary.py`** — Reúne e organiza todos os arquivos em uma estrutura numerada.

Os entregáveis gerados são:

| # | Arquivo | Conteúdo |
|---|---|---|
| 01 | `01_genome_assembly.fasta` | Mitogenoma circularizado completo |
| 02 | `02_coding_genes_nt.fasta` | 13 CDS em nucleotídeos |
| 03 | `03_coding_genes_aa.fasta` | 13 CDS em aminoácidos |
| 04 | `04_ribosomal_genes.fasta` | 12S + 16S rRNA |
| 05 | `05_transport_genes.fasta` | 22 tRNAs |
| 06 | `06_gene_positions.tsv` | Posição física de todos os genes |
| 07 | `07_start_stop_codons.tsv` | Start/Stop codons dos CDS |
| 08 | `08_trna_anticodons.tsv` | Anticódons dos tRNAs |
| 09 | `09_circular_map.svg/pdf` | Mapa circular colorido por grupo funcional |
| 10 | `10_annotation.gff` | Anotação GFF3 completa |
| 11 | `structure_svgs/` | 22 SVGs de tRNA + 2 SVGs de rRNA |
| 12 | `genbank_submission/` | `.gbk` + `.tbl` + `.fsa` para GenBank |
| 13 | `13_gene_order.txt` | Ordem gênica linear |
| 14 | `quality_plots/` | Plots de qualidade do MITOS2 |

Todos os scripts são genéricos — recebem o nome da espécie como parâmetro (`--organism`) e derivam automaticamente o identificador de sequência (e.g., "Anodorhynchus leari" → seqid `A_leari`), tornando o pipeline reutilizável para qualquer táxon.

### 4.12 Automação e Reprodutibilidade

A automação do pipeline e a padronização de sua execução representam os principais diferenciais desta proposta em relação a abordagens convencionais de montagem de genomas mitocondriais. Mais do que reunir ferramentas em sequência, busca-se disponibilizar um fluxo de trabalho que possa ser aplicado por diferentes grupos de pesquisa em distintos ambientes computacionais, assegurando consistência analítica, transparência metodológica e facilidade de reutilização.

Para atingir esse objetivo, cada etapa do pipeline foi encapsulada em containers Docker, garantindo que as ferramentas operem em ambientes controlados e reproduzíveis, independentemente do sistema em que forem executadas (MERKEL, 2014). A orquestração pelo Nextflow automatiza a execução dos processos, gerencia dependências e possibilita paralelização quando aplicável. O sistema ainda permite a retomada da análise a partir do ponto de falha (via `-resume`), evitando reprocessamento desnecessário e aumentando a eficiência computacional (DI TOMMASO et al., 2017).

Todo o código-fonte do workflow está disponibilizado em repositório público no GitHub (https://github.com/matheus-sobreira/mitogenome-pipeline), acompanhado de documentação detalhada, guia de execução e exemplos de uso. Essa prática de ciência aberta não apenas reforça a transparência metodológica, mas também incentiva a adaptação e o aprimoramento do pipeline por outros pesquisadores, em consonância com iniciativas comunitárias como o nf-core (EWELS et al., 2020).

> **Figura 14** — Três camadas de reprodutibilidade do pipeline: (1) código-fonte no GitHub, (2) ambiente isolado via Docker com versões pinadas, (3) orquestração via Nextflow DSL2.

As três camadas representadas na Figura 14 são complementares e funcionam de forma concêntrica: a camada externa (código-fonte no GitHub) garante a acessibilidade e o versionamento; a camada intermediária (Docker com versões pinadas) garante a consistência do ambiente de execução; e a camada interna (Nextflow DSL2) garante a orquestração determinística das etapas.


---

## Capítulo 5 — Resultados e Discussão

Este capítulo apresenta os resultados obtidos com a aplicação do pipeline desenvolvido, organizados segundo as etapas do fluxo de trabalho. Diferentemente da versão anterior deste trabalho, que apresentava resultados esperados de natureza projetiva, os dados aqui descritos são produtos reais de execuções completas do pipeline, validando sua funcionalidade e robustez.

### 5.1 Execuções Realizadas

O pipeline foi executado para três espécies, com objetivos complementares:

| | *Anodorhynchus leari* | *Diploprion bifasciatus* | *Anodorhynchus hyacinthinus* |
|---|---|---|---|
| **Objetivo** | Objeto principal do TCC | Validação do pipeline | Validação cruzada |
| **Dataset SRA** | SRR28399504 | SRR36182901 | SRR36400750 |
| **Plataforma** | Illumina HiSeq X Ten | NovaSeq X Plus | NovaSeq X |
| **Tipo de dados** | WGS | WGS | Metagenômico |
| **Leituras** | 118,5M × 150 bp | 12,0M × 151 bp | 145,8M × 150 bp |
| **Semente** | cox1 *A. hyacinthinus* (NC_082165.1, 1.548 bp) | cox1 *D. bifasciatus* (PZ143763.1, 1.560 bp) | cox1 *A. leari* (montagem do pipeline, 1.548 bp) |
| **`genome_range`** | 15.500–18.500 bp | 15.000–18.500 bp | 15.500–18.500 bp |
| **K-mers** | 39, 33 | 33 | 39, 33 |
| **Código genético** | 2 (vertebrado) | 2 (vertebrado) | 2 (vertebrado) |
| **MITOS2** | Sim (RefSeq89m) | Não (perfil de teste) | Sim (RefSeq89m) |
| **Perfil** | `-profile a_leari` | `-profile test` | `-profile a_hyacinthinus` |

### 5.2 Aquisição e Preparo dos Dados

Para a arara-azul-de-lear, o download via `prefetch` obteve o arquivo `.sra` de 11 GB. O módulo `SRA_DOWNLOAD` utilizou o `fastq-dump -X 25000000` para extrair diretamente apenas os 25 milhões de reads recomendados pelo Pilot QC, gerando dois arquivos FASTQ de aproximadamente 7,4 GB cada (R1 e R2) — sem necessidade de converter o dataset completo de 118,5M reads. O Pilot QC, ativado automaticamente por não haver `sra_max_reads` definido no perfil, analisou 500K reads piloto e recomendou 25M reads como volume ideal para atingir a cobertura-alvo.

> **Figura 15** — Redução de volume em cada etapa do pipeline para a arara-azul-de-lear: tamanho dos arquivos por etapa (esquerda) e composição de bases após trimming (direita).

Na Figura 15, o painel esquerdo evidencia a redução progressiva do volume de dados em cada etapa do pipeline, desde os 87 GB teóricos do dataset completo até os 17 KB da montagem circularizada — uma redução de aproximadamente 5 milhões de vezes. Destaca-se que, com a extração direta via `fastq-dump -X`, o pipeline gera apenas ~14,8 GB de FASTQs (25M reads) em vez dos ~87 GB que seriam produzidos pela conversão integral. O painel direito confirma a distribuição equilibrada das quatro bases nucleotídicas (A, T, G, C) após a trimagem, indicando ausência de viés composicional que pudesse comprometer a montagem.

Para *D. bifasciatus*, o dataset compacto (12M reads, ~1,2 GB) foi processado integralmente sem necessidade de truncagem, utilizando como semente o gene cox1 da própria referência publicada (PZ143763.1, 1.560 bp) e `genome_range` de 15.000–18.500 bp.

### 5.3 Controle de Qualidade e Pré-processamento

Os relatórios do FastQC para a arara-azul-de-lear indicaram qualidade per-base consistente (medianas superiores a Phred 30 ao longo de toda a extensão das leituras), padrão esperado para dados HiSeq X Ten. O Trim Galore identificou adaptadores Illumina TruSeq em 81,4% das leituras e, após o trimming, manteve 71,9% das bases originais. Esse percentual elevado de adaptador é atribuído ao tamanho curto dos insertos da biblioteca, característico de preparações para HiSeq X Ten, em que os fragmentos frequentemente ultrapassam o comprimento das leituras.

> **Figura 16** — Relatórios FastQC para a arara-azul-de-lear. Acima: R1 — qualidade *per-base* (esq.) e conteúdo de adaptadores (dir.). Abaixo: R2 — qualidade *per-base* (esq.) e conteúdo de adaptadores (dir.).

Nos painéis superiores da Figura 16, observa-se que a qualidade per-base das leituras R1 e R2 permanece consistentemente acima do limiar Phred 30 (faixa verde), com uma leve degradação nas posições finais — comportamento esperado na tecnologia Illumina, em que a incorporação de fluoróforos torna-se menos eficiente ao longo da leitura. Os painéis inferiores revelam a presença de adaptadores Illumina TruSeq, cuja curva ascendente a partir da posição ~80 bp confirma a necessidade da etapa de trimagem executada pelo Trim Galore. Não foi executado um segundo FastQC após o Trim Galore, pois a validação da qualidade pós-trimagem é indiretamente assegurada pelo sucesso da montagem: a circularização completa do mitogenoma pelo NOVOPlasty e a anotação íntegra dos 37 genes pelo MITOS2 constituem evidências mais robustas da qualidade dos dados do que um relatório de QC adicional.

### 5.4 Montagem do Mitogenoma

A montagem pelo NOVOPlasty resultou em genomas mitocondriais circularizados para ambas as espécies:

| Métrica | *A. leari* | *D. bifasciatus* | *A. hyacinthinus* |
|---|---|---|---|
| **Tamanho** | 16.986 bp | ~16.800 bp (circularizado) | 16.999 bp |
| **Circularização** | Sim | Sim | Sim |
| **K-mer ótimo** | 39 | 33 | 39 |
| **Cobertura média** | 176× | não medida | não medida |
| **Reads alinhadas** | 19.906 / 11.842.862 (0,17%) | — | — |
| **Insert size medido** | 210 bp (param: 300) | — | — |
| **Referência** | 16.999 bp (*A. hyacinthinus*) — Δ 13 bp | 16.805 bp (PZ143763.1) — Δ ~5 bp | 16.999 bp (OR209186.1) — Δ 0 bp, 19 SNPs |

Para *D. bifasciatus*, o dataset compacto (12M reads, ~1,2 GB) foi processado integralmente sem truncagem; a circularização ocorreu na primeira tentativa com k-mer 33. A etapa de anotação funcional (MITOS2) não foi executada para esta espécie no perfil de teste, que foi configurado sem o parâmetro `mitos2_db` para demonstrar a modularidade do pipeline — a Fase 1 (montagem) opera independentemente da Fase 2 (anotação). A validação da montagem pode ser realizada por BLAST contra a referência PZ143763.1.

O resultado da arara-azul-de-lear merece destaque: o genoma montado de 16.986 bp difere em apenas 13 pares de bases do mitogenoma de referência de *A. hyacinthinus* (16.999 bp, NC_082165.1), demonstrando alta congruência interspecífica esperada para o gênero *Anodorhynchus*. A cobertura média de 176× é amplamente superior ao mínimo recomendado de 30–50× para montagens confiáveis (DIERCKXSENS; MARDULYN; SMITS, 2017), conferindo alta confiabilidade ao resultado.

A proporção de leituras mitocondriais de 0,17% é condizente com a estimativa teórica para dados de genoma total de aves, onde o número de cópias mitocondriais por célula é tipicamente na ordem de centenas a milhares, mas o genoma nuclear (~1,2 Gb em psitacídeos) é muito maior que o mitocondrial (~17 kb).

O log do NOVOPlasty reportou uma fração de subamostragem (*subsampled fraction*) de 27,1%, indicando que o algoritmo processou apenas 27,1% das leituras de entrada (~3,2M de 11,8M reads pós-trimming) antes de obter a circularização. Esse valor reflete a eficiência do *seed-and-extend*: como o NOVOPlasty constrói o genoma iterativamente a partir da semente, ele pode alcançar a montagem completa sem necessidade de processar todas as leituras disponíveis, resultando em economia computacional proporcional à fração não processada.

> **Figura 17** — Mapa circular do mitogenoma da arara-azul-de-lear (*A. leari*), mostrando a posição dos 37 genes (13 CDS, 22 tRNAs, 2 rRNAs) e a região controle (D-loop). Gerado automaticamente pelo pipeline via Biopython GenomeDiagram.

O mapa circular da Figura 17 confirma visualmente a integridade da montagem: os 37 genes estão distribuídos ao longo de toda a molécula sem lacunas, e a região controle (D-loop) ocupa a posição esperada entre os genes tRNA-Pro e tRNA-Phe. As cores distinguem os grupos funcionais — genes codificadores de proteínas, tRNAs e rRNAs — facilitando a identificação da organização gênica típica dos Psittaciformes.

### 5.5 Anotação Funcional

A anotação pelo MITOS2 identificou o conjunto completo de genes esperado para um genoma mitocondrial de ave:

| Tipo | Quantidade | Genes |
|---|---|---|
| **CDS** | 13 | nad1, nad2, nad3, nad4, nad4l, nad5, nad6, cox1, cox2, cox3, atp6, atp8, cob |
| **tRNA** | 22 | tRNA-Phe, tRNA-Val, tRNA-Leu1, tRNA-Ile, tRNA-Gln, tRNA-Met, tRNA-Trp, tRNA-Ala, tRNA-Asn, tRNA-Cys, tRNA-Tyr, tRNA-Ser1, tRNA-Asp, tRNA-Lys, tRNA-Gly, tRNA-Arg, tRNA-His, tRNA-Ser2, tRNA-Leu2, tRNA-Glu, tRNA-Thr, tRNA-Pro |
| **rRNA** | 2 | rrnS (12S), rrnL (16S) |
| **Origem de replicação** | 2 | OH_0, OH_1 |

A ordem gênica identificada é compatível com o padrão conservado de Psittaciformes, corroborando a confiabilidade da montagem. A comparação com a anotação de *A. hyacinthinus* (NC_082165.1) confirmou a conservação da organização gênica dentro do gênero. A completude e a consistência da anotação — 13 CDS, 22 tRNAs com anticódons corretos, 2 rRNAs nas posições esperadas e nenhum gene truncado ou duplicado — servem como validação biológica indireta da montagem: um genoma incorretamente circularizado ou contendo inserções e deleções artifactuais produziria genes com códons de parada prematuros, tRNAs com estruturas secundárias anômalas ou regiões intergênicas atipicamente longas — nenhum dos quais foi observado.

**Caso especial — Frameshift do ND3:** Um aspecto biologicamente relevante revelado pela anotação é o tratamento do gene ND3. O MITOS2 reportou este gene como dois ORFs distintos (nad3_0 e nad3_1) com sobreposição de 5 bp, refletindo o frameshift insertion documentado por Mindell et al. (1998) em genomas mitocondriais de aves. Esse fenômeno, em que uma inserção de ~1 nucleotídeo altera o quadro de leitura no interior do gene, é corrigido in vivo por mecanismos de programmed ribosomal frameshift ou edição do mRNA. Na referência *A. hyacinthinus* (NC_082165.1), o ND3 é anotado com `join()`:

```
CDS   join(9504..9677,9679..9854)
      /gene="ND3"
      /note="frameshift mechanism unknown (Mindell et al., 1998)"
```

Para a conversão dos resultados ao formato GenBank Flat File, o script `gff2genbank.py` do pipeline detecta automaticamente a presença de nad3_0 e nad3_1 e os unifica em um único CDS com `join()`, aplicando o qualificador `/exception=ribosomal slippage` — exatamente como exigido pelo GenBank para submissão. Quando a espécie não apresenta esse frameshift (e.g., invertebrados), o ND3 é anotado normalmente como um CDS contínuo.

> **Figura 18** — Painel com estruturas secundárias dos 22 tRNAs e 2 rRNAs preditos pelo MITOS2/RNAplot para *A. leari*. Nota-se a ausência do braço DHU no tRNA-Ser(AGY), padrão conservado em metazoários.

Na Figura 18, observa-se que todos os tRNAs apresentam a estrutura canônica em trevo (*cloverleaf*), com quatro braços distintos (aceptor, DHU, anticódon e TΨC), exceto o tRNA-Ser(AGY), que exibe a ausência do braço DHU — uma característica conservada em metazoários há mais de 600 milhões de anos. As estruturas dos rRNAs 12S e 16S exibem os padrões de dobramentos esperados, com múltiplas hélices e loops internos. A geração automática dessas estruturas pelo pipeline, integrando MITOS2 e RNAplot, constitui um diferencial em relação às ferramentas analisadas no Capítulo 3, que não incluem predição de estruturas secundárias.

### 5.6 Compilação dos Entregáveis e Formato GenBank

O módulo COMPILE_SUMMARY reuniu automaticamente todos os produtos da análise em uma pasta organizada com 14 categorias de entregáveis, conforme detalhado na Tabela da Seção 4.11. Destaca-se a geração automática do GenBank Flat File (`.gbk`), que é o formato padrão exigido pelo NCBI para submissão de genomas. O arquivo gerado para a arara-azul-de-lear contém:

```
LOCUS       A_leari                16986 bp    DNA     circular VRT 12-APR-2026
DEFINITION  Anodorhynchus leari mitochondrion, complete genome.
ACCESSION   A_leari
VERSION     A_leari
FEATURES             Location/Qualifiers
     source          1..16986
                     /organism="Anodorhynchus leari"
                     /organelle="mitochondrion"
                     /mol_type="genomic DNA"
     ...
     CDS             join(10919..11098,11094..11267)
                     /gene="ND3"
                     /exception="ribosomal slippage"
                     /note="programmed frameshift; Mindell et al. (1998)"
     ...
ORIGIN
        1 gtttacgcta taaagcgtta gatataactg ...
//
```

Esse formato pode ser diretamente submetido ao GenBank via `table2asn` ou carregado em ferramentas de visualização como OGDRAW para geração de mapas circulares de alta qualidade.

### 5.7 Análise de Desempenho e Eficiência Computacional

A avaliação de desempenho de um pipeline bioinformático é essencial para determinar sua viabilidade prática em diferentes cenários de uso — desde laptops pessoais até servidores de computação em nuvem. Em bioinformática, pipelines frequentemente processam volumes massivos de dados (dezenas a centenas de gigabytes), e compreender onde o tempo e os recursos computacionais são consumidos permite identificar gargalos, dimensionar infraestrutura adequada e estimar custos operacionais. Essa análise é particularmente relevante para laboratórios com recursos limitados, que precisam avaliar se dispõem de hardware suficiente antes de iniciar uma execução.

Nesta seção, são apresentadas métricas detalhadas de execução do pipeline para a arara-azul-de-lear, incluindo tempo, consumo de CPU, memória RAM e operações de entrada/saída (I/O) — isto é, volume de dados lidos e escritos em disco — por etapa. A partir dessas métricas, discutem-se a eficiência de cada módulo, a identificação de gargalos, o impacto do pré-processamento na qualidade da montagem, e o custo estimado de execução em infraestrutura de nuvem.

A execução foi realizada em um laptop com processador Intel Core i7-1165G7 (4 núcleos físicos, 8 threads via *hyperthreading*, frequência base de 2,80 GHz), 16 GB de memória RAM DDR4 e armazenamento em estado sólido (SSD NVMe). O tempo total decorrido — também chamado *wall-clock time*, que mede o tempo real percebido pelo usuário, incluindo esperas por I/O e rede — foi de **29 minutos e 29 segundos**, com consumo acumulado de **2,1 horas de CPU** (tempo efetivo de processamento somado de todos os núcleos). A Tabela 5 apresenta as métricas detalhadas por etapa.

**Tabela 5 — Métricas de desempenho por etapa do pipeline (execução arara-azul-de-lear, 25M reads)**

| Etapa | Tempo (min) | % do total | CPU (%) | RAM pico (GB) | I/O leitura | I/O escrita |
|-------|-------------|------------|---------|---------------|-------------|-------------|
| SRA_PILOT | 2,2 | 6,5% | 9% | 0,15 | 7,2 GB | 744 MB |
| SRA_DOWNLOAD | 15,7 | 47,5% | 19% | 0,10 | 12,9 GB | 31,6 GB |
| FASTQC | 3,6 | 10,8% | 200% | 1,10 | 20,9 GB | 3,9 MB |
| TRIM_GALORE | 5,8 | 17,6% | 217% | 0,12 | 83,1 GB | 74,8 GB |
| NOVOPLASTY | 2,0 | 6,1% | 101% | 6,80 | 8,6 GB | 43 KB |
| MITOS2 | 3,8 | 11,4% | 103% | 0,47 | 158 MB | 43 MB |
| **Total** | **33,0** | **100%** | — | — | — | — |

> **Figura 19** — Tempo de execução por etapa do pipeline para a arara-azul-de-lear, com anotação da porcentagem do tempo total e classificação por tipo de operação (I/O-bound vs CPU-bound).

A Figura 19 evidencia que o gargalo dominante é o download do SRA (47,5% do tempo total), uma etapa limitada pela largura de banda da rede e não pela capacidade computacional local. As etapas de processamento intensivo (Trim Galore e NOVOPlasty) representam, conjuntamente, apenas 24% do tempo, demonstrando que a montagem de mitogenomas não exige alto poder de processamento quando o volume de dados é adequadamente reduzido pelo Pilot QC.

#### 5.7.1 Gargalo de I/O: Download e Conversão SRA

O download do SRA é o gargalo dominante, consumindo **47,5% do tempo total** (15,7 minutos). A etapa engloba duas operações sequenciais: `prefetch` (download do arquivo `.sra` comprimido, ~11 GB) e `fastq-dump -X 25000000` (extração direta dos primeiros 25M spots para FASTQ, gerando ~14,8 GB). O volume total de I/O desta etapa (12,9 GB lidos + 31,6 GB escritos) é significativamente inferior ao que seria produzido pela conversão integral do dataset de 118M reads (~87 GB de FASTQ). A adoção do `fastq-dump -X`, que extrai diretamente apenas os spots necessários, evita a gravação desnecessária e reduz tanto o tempo quanto o consumo de armazenamento.

Esse comportamento reforça a importância do módulo Pilot QC, que utiliza `fastq-dump -X 500000` para amostrar apenas as primeiras 500K reads diretamente da origem (sem download completo), consumindo apenas 2,2 minutos e 744 MB. Sem o Pilot, o usuário precisaria decidir arbitrariamente quantas reads baixar, arriscando sub ou superamostragem.

#### 5.7.2 Eficiência do Pré-processamento

O Trim Galore consumiu 5,8 minutos (17,6% do tempo), processando 50 milhões de reads (25M pares). A taxa de remoção de adaptadores foi de 81,4% (R1), com remoção de bases por qualidade inferior a 1% (37M bp em R1). Do total de 25M pares, 1.853.399 (7,41%) foram removidos por ficarem abaixo de 50 bp após trimming — indicando dímeros de adaptador. O volume de dados foi reduzido de 7,5 Gbp brutos para 5,39 Gbp filtrados (71,9% de retenção em bases).

Uma questão pertinente é: **o trimming é necessário ou é perda de tempo?** A resposta é inequívoca: com 81,4% dos reads contendo adaptador Illumina TruSeq, executar a montagem sem trimming resultaria em k-mers quiméricos (parte sequência biológica, parte adaptador sintético) que impediriam a extensão correta das sobreposições no algoritmo *seed-and-extend*. Os 5,8 minutos de trimming (~18% do pipeline) previnem artefatos de montagem que poderiam resultar em falha de circularização ou contigs espúrios. Em contraste, a montagem pós-trimming completou em apenas 2 minutos, com circularização na primeira tentativa.

O FASTQC rodou em paralelo com o Trim Galore (ambos consomem os mesmos FASTQs), não adicionando tempo ao caminho crítico. O pico de memória do FASTQC (1,1 GB) é significativo para máquinas com pouca RAM, mas justificável pela geração de relatórios de qualidade essenciais para verificação posterior.

#### 5.7.3 Consumo de Memória do NOVOPlasty

O NOVOPlasty apresentou o maior consumo de memória entre todas as etapas: **6,8 GB de RAM** (pico). Esse consumo é diretamente proporcional ao volume de leituras de entrada: com 25M reads (23,1M após trimming), a construção do *hash* de k-mers em memória resulta em um pico compatível com computadores de 8 GB de RAM. O NOVOPlasty, escrito em Perl, não possui gerenciamento granular de memória, e o parâmetro `Max memory` (configurado automaticamente como 7 GB com base na RAM disponível) serve apenas como limite superior sugerido, não como restrição efetiva.

A fração "subsampled" de 27,1% reportada no log indica que o NOVOPlasty sub-amostrou internamente os dados de entrada, processando efetivamente ~3,2M reads dos 11,8M disponíveis. Mesmo com essa sub-amostragem, a cobertura final de 176× é 3,5× superior ao mínimo recomendado de 50× (DIERCKXSENS; MARDULYN; SMITS, 2017), conferindo alta confiabilidade à montagem.

**Recomendação prática:** para execução em máquinas com 8 GB de RAM, o pipeline é viável sem ajustes adicionais: o pico de 6,8 GB de RAM do NOVOPlasty é compatível com o limite de 7 GB auto-detectado, deixando margem para o sistema operacional e Docker.

#### 5.7.4 Cascata de Redução de Dados

Uma análise particularmente reveladora é a redução progressiva do volume de dados ao longo do pipeline:

**Tabela 6 — Cascata de redução de dados ao longo do pipeline**

| Etapa | Volume | Fator de redução | Acumulado |
|-------|--------|-------------------|-----------|
| SRA completo (118,5M reads) | ~87 GB | — | 1× |
| Após `fastq-dump -X` (25M reads) | ~14,8 GB | 5,9× | 5,9× |
| Após trimming (23,1M reads) | ~5,4 GB | 2,7× | 16,1× |
| Montagem circularizada | 17 KB | 317.000× | ~5.100.000× |
| Summary (todos entregáveis) | 1,6 MB | — | ~54.000× |

O dado original (87 GB de FASTQ bruto) é reduzido em mais de **cinco milhões de vezes** até a montagem final (17 KB). Essa cascata demonstra que a maior parte do sequenciamento WGS é irrelevante para o genoma mitocondrial — apenas 0,17% das reads alinham ao mitogenoma. A existência do Pilot QC como primeira etapa evita o download desnecessário: em vez de baixar os 87 GB completos, o pipeline determina automaticamente que 25M reads (~14,8 GB) são suficientes, extraindo-os diretamente com `fastq-dump -X` sem necessidade de conversão integral.

> **Figura 20** — Cascata de redução de dados ao longo do pipeline: de 87 GB de FASTQ bruto até 17 KB de genoma circularizado — uma redução de aproximadamente 5 milhões de vezes.

O formato de funil da Figura 20 ilustra a lógica central do pipeline: cada etapa descarta dados irrelevantes para a montagem mitocondrial, reduzindo progressivamente o volume processado. O salto mais expressivo ocorre na etapa de montagem, em que o NOVOPlasty seleciona, dentre milhões de leituras, apenas aquelas que pertencem ao genoma mitocondrial — reduzindo 5,4 GB de leituras trimadas a um genoma de 17 KB. Essa redução de mais de cinco milhões de vezes viabiliza a execução completa do pipeline em um notebook com 8 GB de RAM.

#### 5.7.5 Pico de Armazenamento e Custo Computacional

Um aspecto frequentemente negligenciado em pipelines bioinformáticos é o **consumo de armazenamento durante a execução**, que pode ser significativamente maior que o tamanho dos resultados finais. Isso ocorre porque as etapas intermediárias geram arquivos temporários volumosos (FASTQs descomprimidos, arquivos de cache) que são consumidos pelas etapas seguintes e depois descartados.

No caso da execução para a arara-azul-de-lear, o pico de uso de disco foi de **~16 GB** no diretório de trabalho, dominado pela etapa SRA_DOWNLOAD (~14,8 GB de FASTQs descomprimidos e cache `.sra`). Em contraste, ao final da execução, os resultados úteis ocupam apenas **~67 MB** em `results/`. Os arquivos intermediários podem ser removidos com o comando `nextflow clean`, liberando o espaço ocupado. Essa diferença de ~240× entre o pico de armazenamento e os resultados finais é uma característica importante a considerar ao dimensionar a infraestrutura.

**Custo computacional em hardware doméstico.** Todo o pipeline foi executado em um notebook pessoal com processador de 4 núcleos, 16 GB de RAM e disco SSD convencional. Os requisitos computacionais reais observados foram:

- **Processamento:** 33,0 minutos de tempo total (*wall-clock*), com pico de uso de 4 threads durante o trimming. Nenhuma etapa exigiu processamento de alto desempenho — o processo mais intensivo (NOVOPlasty) utilizou apenas 1 thread.
- **Armazenamento:** ~16 GB de pico no diretório de trabalho, reduzidos a ~67 MB de resultados finais após limpeza dos arquivos temporários. Um disco com pelo menos 30 GB livres é suficiente com margem de segurança.
- **Memória RAM:** pico de 6,8 GB durante a montagem com NOVOPlasty, dentro dos 16 GB disponíveis. Um computador com 8 GB de RAM executa o pipeline sem dificuldade, pois o pico de 6,8 GB é compatível com o limite de 7 GB auto-detectado pelo NOVOPlasty.
- **Rede:** o download do arquivo SRA (~11 GB) foi o principal consumidor de banda, com duração variável conforme a velocidade da conexão.

Esses requisitos modestos demonstram que **a montagem e anotação de mitogenomas não exige infraestrutura computacional especializada**. Qualquer computador pessoal de categoria intermediária, com 8 GB de RAM, disco SSD e acesso à internet, é capaz de executar o pipeline completo em pouco mais de 30 minutos. Essa acessibilidade é particularmente relevante para estudantes e grupos de pesquisa com recursos limitados, que podem reproduzir integralmente os resultados deste trabalho sem investimento adicional em hardware.

#### 5.7.6 Caminho Crítico e Otimização Implementada

O **caminho crítico** de um pipeline — ou seja, a sequência de etapas dependentes que determina o tempo mínimo de execução — é: SRA_PILOT → SRA_DOWNLOAD → TRIM_GALORE → NOVOPLASTY → MITOS2 → COMPILE_SUMMARY. Etapas que executam em paralelo, como o FASTQC (que processa os mesmos FASTQs que o Trim Galore ao mesmo tempo), **não adicionam tempo ao caminho crítico**, exceto quando excedem a duração da etapa paralela. No caso observado, o FASTQC (3,6 min) completou antes do Trim Galore (5,8 min) e, portanto, não impactou o tempo total.

A análise do caminho crítico revela que **65,0% do tempo total** (21,5 de 33,0 min) é consumido por apenas duas etapas: download SRA (15,7 min) e trimming (5,8 min). A montagem em si — o objetivo central do pipeline — leva apenas 2,0 minutos (6,1%), demonstrando que o gargalo real não é computacional, mas de transferência e transformação de dados.

A otimização implementada — substituição do `fasterq-dump` + truncamento por `fastq-dump -X`, combinada com a redução automática do volume de reads pelo Pilot QC — é uma decisão arquitetural com impacto direto no desempenho. Sem essa otimização, o pipeline precisaria converter integralmente o arquivo `.sra` de 118,5M reads para FASTQ (~87 GB) antes de truncar, gerando centenas de gigabytes de escrita em disco e aumentando significativamente o tempo de download. Com o `fastq-dump -X`, apenas os 25M reads necessários são extraídos diretamente, reduzindo a escrita em disco para 31,6 GB e o tempo de download para 15,7 min.

#### 5.7.7 Síntese da Análise de Desempenho

A análise de desempenho revelou três conclusões principais sobre o comportamento computacional do pipeline:

**O gargalo é I/O, não processamento.** Das sete etapas do pipeline, a montagem do mitogenoma — objetivo central do trabalho — consome apenas 6,1% do tempo total (2,0 min). O restante é dominado por operações de transferência e transformação de dados: o download SRA representa 47,6% (15,7 min) e o trimming responde por 17,6% (5,8 min). Em termos práticos, isso significa que melhorias na velocidade de rede ou na estratégia de amostragem têm impacto muito maior no tempo total do que otimizações no algoritmo de montagem.

**A etapa de trimming justifica seu custo.** Apesar de consumir 5,8 minutos, o trimming não é dispensável: 81,4% das reads continham sequências de adaptadores que, se mantidas, poderiam gerar extensões quiméricas no grafo do NOVOPlasty. Os 5,8 minutos investidos previnem falhas de circularização que exigiriam re-execução completa do pipeline — um custo muito superior.

**Sem as estratégias implementadas no pipeline, a montagem de um mitogenoma a partir de dados públicos seria inviável em hardware doméstico.** O dataset original de *A. leari* (SRR28399504) contém 118,5 milhões de pares de leituras, que, quando convertidos para o formato FASTQ, ocupam aproximadamente 87 GB de espaço em disco. Carregar e processar esse volume na íntegra — como seria necessário em uma abordagem manual sem truncamento — exigiria centenas de gigabytes de armazenamento temporário (leitura, escrita intermediária e saída), além de memória RAM proporcional ao volume de dados nos estágios de trimming e montagem. Um notebook convencional com 16 GB de RAM e 256 GB de SSD simplesmente não comportaria essa carga; seria necessário recorrer a um servidor dedicado ou a uma infraestrutura de computação em nuvem.

O pipeline resolve esse problema por meio de uma cadeia de reduções progressivas que viabilizam a execução em hardware modesto. O **Pilot QC** analisa uma amostra mínima (500 mil leituras, ~500 MB) e calcula exatamente quantas leituras são necessárias para atingir a cobertura-alvo, evitando o download do dataset inteiro. O **`fastq-dump -X`** extrai diretamente apenas os 25 milhões de reads necessários, reduzindo o volume de 87 GB para ~14,8 GB sem necessidade de conversão integral do arquivo SRA. O **trimming** remove adaptadores e bases de baixa qualidade, reduzindo para ~5,4 GB de dados úteis. Finalmente, o **NOVOPlasty** sub-amostra internamente e opera apenas sobre as leituras que alinham à semente, consumindo uma fração do total. O resultado: de 87 GB de dados brutos, o pipeline extrai um mitogenoma de **17 KB** — uma redução de aproximadamente 5 milhões de vezes.

Essa cadeia de decisões é o que torna possível executar o ciclo completo — do dado bruto público ao mitogenoma anotado com 14 categorias de entregáveis — em aproximadamente 33 minutos, em um notebook com 4 núcleos e 8 GB de RAM, sem necessidade de infraestrutura especializada. A viabilidade em hardware doméstico é particularmente relevante para estudantes de graduação e pesquisadores em instituições com recursos computacionais limitados, que de outra forma dependeriam de acesso a servidores ou clusters para realizar análises genômicas.

### 5.8 Discussão sobre a Robustez da Abordagem

A execução bem-sucedida do pipeline para três espécies distintas — duas aves congêneres (Psittaciformes) e um peixe (Perciformes) —, utilizando datasets de plataformas diferentes (HiSeq X Ten, NovaSeq X Plus e NovaSeq X) e tipos de dados distintos (WGS e metagenômico), demonstra a generalização da abordagem. A obtenção de genomas circularizados em todos os casos, com o conjunto completo de 37 genes identificados nos casos anotados, valida tanto o fluxo de montagem quanto a integração de ferramentas.

A cobertura de 176× obtida para *A. leari* com apenas 0,17% de leituras mitocondriais em 25 milhões de reads demonstra que o DNA mitocondrial é eficientemente recuperado pelo NOVOPlasty mesmo quando representa uma fração diminuta do total de dados, confirmando as observações da literatura (DIERCKXSENS; MARDULYN; SMITS, 2017).

A combinação metodológica adotada neste trabalho — NOVOPlasty para montagem e MITOS2 para anotação — corresponde ao padrão predominante na literatura recente de genômica mitocondrial, empregada em estudos abrangendo desde 213 mitogenomas de lagartos (ZHAN et al., 2024) até espécies individuais de peixes (KUNDU et al., 2024), arraias (GUERREIRO et al., 2025) e invertebrados (SHI et al., 2024; TAO et al., 2024). Essa convergência metodológica independente — em que grupos de pesquisa de diferentes países, trabalhando com táxons distintos, adotam as mesmas ferramentas — constitui uma validação externa da abordagem e reforça a confiabilidade das escolhas implementadas no pipeline.

A lógica iterativa de montagem com re-seeding automático e múltiplos k-mers mostrou-se eficaz, circularizando na primeira tentativa (k-mer 39) para *A. leari*, na primeira tentativa (k-mer 33) para *D. bifasciatus* e na primeira tentativa (k-mer 39) para *A. hyacinthinus*. Essa robustez é particularmente relevante para cenários em que a *seed* disponível é de uma espécie distante, ou quando regiões repetitivas dificultam a circularização.

O módulo Pilot QC foi validado por meio de execução isolada do script `pilot_qc.sh` sobre uma subamostra de 100.000 leituras de *A. leari*, extraída do dataset SRR28399504. A análise piloto detectou qualidade Q30 de 86,2%, taxa de adaptador de 66,4% e fração mitocondrial estimada de 0,17% (via correspondência de k-mers com a semente cox1). Com esses valores, a recomendação automática foi de 36 milhões de leituras para atingir 500× de cobertura — resultado coerente com os 176× efetivamente obtidos com 25 milhões de reads na execução principal. A configuração de produção do perfil `a_leari` não define `sra_max_reads`, ativando o Pilot QC automaticamente em execuções futuras. Essa funcionalidade é particularmente valiosa para pesquisadores trabalhando com espécies pouco caracterizadas, onde não há estimativa prévia da fração mitocondrial.

Em comparação com os trabalhos relacionados, o pipeline desenvolvido apresenta vantagens claras em termos de reprodutibilidade, portabilidade e acessibilidade computacional. Enquanto o MitoHiFi (ULIANO-SILVA et al., 2023) utiliza Conda para gerenciamento de dependências, sujeito a variações de ambiente, e o MITObim (HAHN et al., 2013) e o GetOrganelle (JIN et al., 2020) dependem de configuração manual, o pipeline aqui proposto encapsula todas as ferramentas em containers Docker com versões pinadas e orquestra sua execução via Nextflow, garantindo resultados idênticos em qualquer ambiente. Adicionalmente, as soluções existentes não incorporam mecanismos de redução automática do volume de dados, pressupondo que o usuário disponha de infraestrutura com armazenamento e memória suficientes para processar datasets completos — requisito que exclui a maioria dos computadores pessoais. O pipeline aqui proposto, ao integrar o Pilot QC e a extração seletiva com `fastq-dump -X`, elimina essa barreira e permite que o ciclo completo — do dado bruto ao mitogenoma anotado — seja executado em um notebook convencional. A adição da Fase 2 (anotação funcional com MITOS2 e geração automática de entregáveis em formato GenBank) vai além do que a maioria dos pipelines de montagem de organelas oferece, cobrindo o ciclo completo desde o dado bruto até a submissão ao banco de dados.

Os resultados metodológicos apresentados neste capítulo demonstram que o pipeline é funcional, eficiente e acessível.

### 5.9 Validação Cruzada com *Anodorhynchus hyacinthinus*

A terceira execução do pipeline constituiu um teste de **validação cruzada** particularmente robusto: a montagem do mitogenoma de *A. hyacinthinus* a partir de um dataset metagenômico (SRR36400750), utilizando como semente o gene cox1 extraído da montagem de *A. leari* — inversamente ao que foi feito na execução principal, em que a semente de *A. leari* veio de *A. hyacinthinus*. Essa configuração recíproca testa a capacidade do pipeline de operar com sementes de espécies congêneres em ambas as direções, e a existência de um mitogenoma de referência depositado no GenBank (OR209186.1, 16.999 bp) permite validação quantitativa nucleotídeo a nucleotídeo.

#### 5.9.1 Desafio: Dataset Metagenômico de Grande Escala

O dataset SRR36400750 apresenta dois desafios adicionais em relação às execuções anteriores. Primeiro, o **volume**: 145,8 milhões de spots (291,7M reads), totalizando 44 Gbases — o maior dataset processado pelo pipeline. Segundo, a **natureza metagenômica**: o DNA foi extraído de fluido pericárdico da ave, contendo uma mistura de DNA do hospedeiro, DNA bacteriano e potencialmente DNA viral, o que implica uma fração mitocondrial ainda menor do que em dados WGS convencionais. Esse cenário testa a capacidade do NOVOPlasty de isolar e montar o genoma mitocondrial a partir de um fundo genômico complexo e heterogêneo.

A execução revelou uma limitação importante no módulo de download original do pipeline, que utilizava `fasterq-dump` para converter o arquivo `.sra` integralmente para FASTQ antes de truncar com `head`. Para um dataset de 145,8M spots, essa abordagem gerava ~124 GB de arquivos FASTQ intermediários — excedendo o armazenamento disponível. O problema foi corrigido pela adoção de `fastq-dump -X` quando o parâmetro `sra_max_reads` está definido, que extrai diretamente apenas os primeiros *N* spots do arquivo `.sra` sem necessidade de gravação integral (Seção 5.9.4).

#### 5.9.2 Comparação com a Referência NCBI

A montagem pelo NOVOPlasty produziu um genoma circularizado de **16.999 bp** — *tamanho idêntico* à referência OR209186.1. A comparação nucleotídeo a nucleotídeo requereu o alinhamento por rotação, uma vez que o mitogenoma é circular e o NOVOPlasty inicia a linearização em uma posição arbitrária determinada pela semente. A montagem obtida inicia na posição 15.572 da referência NCBI, correspondente à região de controle (D-loop), enquanto a sequência de referência inicia no tRNA-Phe — refletindo convenções distintas de linearização.

**Tabela 7 — Comparação nucleotídeo a nucleotídeo entre a montagem do pipeline e a referência NCBI para *A. hyacinthinus*.**

| Métrica | Valor |
|---------|-------|
| Tamanho da montagem | 16.999 bp |
| Tamanho da referência (OR209186.1) | 16.999 bp |
| Diferença de tamanho | 0 bp |
| Offset de rotação | 15.571 bp (início na D-loop) |
| Identidade estrita (bases idênticas) | 16.979 / 16.999 (99,882%) |
| Mismatches reais (incompatíveis) | 19 |
| Diferenças por ambiguidade (compatíveis) | 1 |
| Identidade efetiva (excluindo ambiguidades) | 16.980 / 16.999 (99,888%) |

Das 19 diferenças nucleotídicas observadas, 11 concentram-se na região de controle (D-loop, posições 15.572–16.999 da referência) — a região mais variável do genoma mitocondrial, onde a taxa de mutação é sabidamente 2–5 vezes superior às regiões codificantes (AQUADRO; GREENBERG, 1983). As 8 diferenças restantes distribuem-se ao longo dos ~15 kb de genes codificantes e são predominantemente transições (A↔G ou C↔T), padrão consistente com **variação intraespecífica** entre indivíduos distintos da mesma espécie, e não com erros de montagem. A única diferença por ambiguidade (NCBI=T, montagem=Y na posição 15.730) indica incerteza de sequenciamento na referência original.

Esse resultado — 99,888% de identidade com uma referência independente, com diferenças concentradas na região hipervariável e consistentes com polimorfismo intraespecífico — constitui a **validação mais rigorosa do pipeline**: a montagem *de novo* a partir de dados metagenômicos reproduz, nucleotídeo a nucleotídeo, uma sequência obtida por outro grupo de pesquisa com metodologia, amostras e plataformas distintas.

> **Figura 21** — Mapa circular de polimorfismos entre a montagem do pipeline e a referência NCBI (OR209186.1) para *A. hyacinthinus*. No anel interno, marcadores vermelhos indicam SNPs na região de controle (D-loop, 11/19) e marcadores azuis indicam SNPs em regiões codificantes (8/19). O triângulo amarelo indica a única posição com código de ambiguidade. A concentração de SNPs na D-loop é consistente com a maior taxa de evolução dessa região em genomas mitocondriais.

Na Figura 21, o anel externo representa os genes do mitogenoma, coloridos por complexo funcional, enquanto o anel interno exibe a posição de cada polimorfismo detectado. A assimetria na distribuição dos SNPs é visualmente evidente: os 11 marcadores vermelhos (D-loop) agrupam-se em um arco restrito correspondente a menos de 10% do genoma, enquanto os 8 marcadores azuis (genes codificantes e rRNAs) distribuem-se ao longo dos ~90% restantes.

#### 5.9.3 Otimização do Módulo de Download para Datasets de Grande Escala

A execução com o dataset metagenômico de 145,8M spots revelou um gargalo no módulo `SRA_DOWNLOAD`: o `fasterq-dump`, embora mais rápido que o `fastq-dump` por ser multi-thread, não oferece opção nativa para limitar o número de spots extraídos. A abordagem anterior — extrair todo o dataset e truncar com `head` — era viável para datasets de 118,5M spots (~87 GB de FASTQ), mas impraticável para 145,8M spots (~124 GB), que excediam o armazenamento disponível.

A solução implementada utiliza condicionalmente o `fastq-dump` com a flag `-X` (maxSpotId) quando `sra_max_reads` está definido, extraindo diretamente apenas os primeiros *N* spots do arquivo `.sra`. Embora o `fastq-dump` seja single-thread, a redução do volume processado (25M spots vs. 145,8M) compensa amplamente a perda de paralelismo: a extração processa apenas 17% dos dados, com uso de disco proporcional (~7 GB em vez de ~124 GB). Quando `sra_max_reads` não está definido (i.e., sem limite), o pipeline continua utilizando o `fasterq-dump` multi-thread para máximo desempenho.

#### 5.9.4 Significado da Validação Cruzada

A execução com *A. hyacinthinus* valida o pipeline em três dimensões complementares:

1. **Correção biológica**: a montagem reproduz uma referência independente com 99,888% de identidade, com diferenças explicáveis por variação natural.
2. **Robustez operacional**: o pipeline processa com sucesso um dataset metagenômico de grande escala (145,8M spots), demonstrando capacidade de isolar o genoma mitocondrial de um fundo genômico complexo.
3. **Reciprocidade de sementes**: a semente de *A. leari* funciona para montar *A. hyacinthinus* (e vice-versa), validando a estratégia de uso de sementes congêneres quando a espécie-alvo não possui sequência de referência própria.

O capítulo seguinte dedica-se integralmente à espécie-alvo deste trabalho — a arara-azul-de-lear —, explorando em profundidade o mitogenoma obtido, sua organização gênica e as implicações para a conservação da espécie.

---

## Capítulo 6 — O Mitogenoma da Arara-azul-de-lear

### 6.1 A Espécie

A arara-azul-de-lear (*Anodorhynchus leari* Bonaparte, 1856) é um psitacídeo de grande porte endêmico do bioma Caatinga, no nordeste do Brasil. Pertencente à ordem Psittaciformes, família Psittacidae, gênero *Anodorhynchus*, a espécie compartilha o gênero com apenas duas outras araras: a arara-azul-grande (*A. hyacinthinus*) e a provavelmente extinta arara-azul-pequena (*A. glaucus*) (COLLAR, 1997). A descoberta formal da espécie para a ciência ocorreu em 1856, mas sua área de ocorrência em vida livre só foi localizada por Helmut Sick em 1978, no Raso da Catarina, norte do estado da Bahia (SICK, 1979).

A distribuição geográfica da arara-azul-de-lear é extremamente restrita. A espécie ocorre em uma área de aproximadamente 1.500 km² no semiárido baiano, com núcleos populacionais concentrados nas regiões de Toca Velha e Serra Branca, nos municípios de Canudos, Jeremoabo e Euclides da Cunha (ICMBio, 2022). Os indivíduos nidificam em paredões de arenito, onde escavam cavidades para reprodução, e se alimentam predominantemente de frutos de palmeira licuri (*Syagrus coronata*), estabelecendo uma relação ecológica de forte dependência com essa espécie vegetal (SICK, 1979; ICMBio, 2022).

A população estimada atualmente é de cerca de 1.700 indivíduos em vida livre (ICMBio, 2022), representando uma recuperação significativa em relação às estimativas de aproximadamente 200 indivíduos na década de 1980 (SICK, 1979; ICMBio, 2022). Apesar dessa tendência positiva, a espécie permanece classificada como **Em Perigo (EN)** pela União Internacional para a Conservação da Natureza (IUCN) e consta no Plano de Ação Nacional para a Conservação da Arara-azul-de-lear — PAN Arara-azul-de-lear (ICMBio, 2022), que coordena esforços de proteção do habitat, combate ao tráfico e monitoramento populacional.

As principais ameaças à espécie incluem: (i) o tráfico ilegal de filhotes para o comércio de animais de estimação — historicamente a pressão mais severa, com indivíduos alcançando valores extremamente elevados no mercado internacional (ICMBio, 2022); (ii) a degradação e fragmentação do habitat, especialmente pela remoção de palmeiras licuri para exploração da cera e do coco; (iii) a baixa variabilidade genética decorrente do gargalo populacional histórico; e (iv) a área de distribuição extremamente restrita, que torna a espécie vulnerável a eventos estocásticos como secas prolongadas ou epizootias (surtos de doenças infecciosas que acometem populações animais, análogos a epidemias em humanos).

Nesse contexto, a disponibilização de dados genômicos de alta qualidade torna-se uma ferramenta estratégica para a conservação. A genômica de conservação utiliza dados moleculares para subsidiar decisões de manejo, como a avaliação da diversidade genética remanescente, a identificação de unidades evolutivamente significativas e o desenvolvimento de marcadores para identificação forense de espécimes apreendidos do tráfico. O sequenciamento do genoma mitocondrial completo da arara-azul-de-lear representa uma contribuição direta para esse campo.

> **[SUGESTÃO DE FIGURA — Distribuição Geográfica]**: Mapa da distribuição geográfica da arara-azul-de-lear no nordeste da Bahia (Raso da Catarina, Boqueirão da Onça, região de Canudos). Indicar os paredões de arenito de nidificação e as áreas de ocorrência do licuri (*Syagrus coronata*). Pode-se basear nos mapas do PAN ICMBio (2022) com devida atribuição.

> **[SUGESTÃO DE FIGURA — Fotografia da Espécie]**: Fotografia da arara-azul-de-lear em ambiente natural, evidenciando a plumagem azul-cobalto diagnóstica e a mancha perioftálmica amarela que a distingue de *A. hyacinthinus*. Fonte sugerida: acervo ICMBio ou Wikimedia Commons (licença CC).

### 6.2 Caracterização Geral do Mitogenoma

O mitogenoma da arara-azul-de-lear obtido neste trabalho possui **16.986 pares de bases** e apresenta topologia circular — resultado confirmado pela convergência das extremidades da sequência durante a montagem (circularização). A molécula segue a organização típica de mitogenomas de vertebrados, com 37 genes funcionais e uma região controle (D-loop).

A composição nucleotídica apresenta as seguintes proporções:

| Base | Quantidade | Proporção |
|------|-----------|-----------|
| Adenina (A) | 5.153 | 30,3% |
| Timina (T) | 3.978 | 23,4% |
| Guanina (G) | 2.372 | 14,0% |
| Citosina (C) | 5.474 | 32,2% |
| **Total** | **16.986** | **100%** |

O conteúdo GC é de **46,2%**. Em mitogenomas, o DNA de fita dupla circular possui duas fitas com composição nucleotídica distinta: a **fita pesada** (H-strand, do inglês *heavy*) é mais rica em purinas (A e G), enquanto a **fita leve** (L-strand, *light*) é mais rica em pirimidinas (T e C) — essa nomenclatura histórica refere-se à diferença de densidade observada em ultracentrifugação. A assimetria composicional entre as fitas é quantificada por duas métricas: o **AT-skew**, calculado como $(A - T) / (A + T)$, e o **GC-skew**, calculado como $(G - C) / (G + C)$. Na arara-azul-de-lear, o AT-skew positivo (+0,129) indica excesso de adenina sobre timina na fita pesada, enquanto o GC-skew negativo (−0,395) indica forte excesso de citosina sobre guanina. Esses padrões de assimetria, chamados de *strand bias*, são consequência dos mecanismos de replicação assimétrica do DNA mitocondrial, que expõem a fita pesada como fita simples por períodos prolongados, tornando-a mais suscetível a deaminações espontâneas de citosina e adenina (BOORE, 1999).

Esses valores são condizentes com padrões observados em outros vertebrados. Estudos recentes de genômica mitocondrial comparativa demonstram que o AT-skew positivo e o GC-skew negativo na fita pesada constituem uma assinatura universal de mitogenomas de vertebrados, embora com variações de magnitude: arraias de água doce (*Potamotrygon leopoldi*, 17.504 bp) apresentam AT-skew de +0,13 e GC-skew de −0,40 (GUERREIRO et al., 2025) — valores muito próximos aos da arara-azul-de-lear (+0,129 e −0,395); peixes teléosteos como *Cephalopholis taeniops* apresentam AT-skew de +0,052 e GC-skew de −0,277 (KUNDU et al., 2024); anfibios como hilideo *Dryophytes japonicus* mostram AT-skew próximo a zero (−0,001) e GC-skew de −0,273 (HONG et al., 2024). A semelhanca notavel entre os valores de skew da arara-azul-de-lear e da arraia *P. leopoldi* — uma espécie igualmente ameaçada e endêmica brasileira — ilustra como as pressões seletivas sobre a composição nucleotídica mitocondrial são conservadas independentemente da história de vida e da posição filogenética do organismo.

O conteúdo gênico completo é composto por:

| Categoria | Quantidade | Detalhes |
|-----------|-----------|---------|
| Genes codificadores de proteínas | 13 | ND1-6, ND4L, COX1-3, ATP6, ATP8, CYTB |
| RNAs transportadores (tRNA) | 22 | Um para cada aminoácido, com exceção de leucina e serina (dois cada) |
| RNAs ribossomais (rRNA) | 2 | 12S rRNA (rrnS, 970 bp) e 16S rRNA (rrnL, 1.572 bp) |
| Região controle (D-loop) | 2 segmentos | OH_0 (262 bp) e OH_1 (47 bp) |

A **região controle**, também chamada **D-loop** (do inglês *displacement loop*), é o único segmento não codificante significativo do mitogenoma. Embora não codifique genes, essa região é essencial porque contém as sequências promotoras que controlam o início da replicação do DNA mitocondrial e da transcrição dos genes (BOORE, 1999). A D-loop é também a região com maior taxa de variação entre indivíduos da mesma espécie, o que a torna um marcador particularmente útil para estudos de diversidade genética intrapopulacional. |

> **[SUGESTÃO DE FIGURA]**: Mapa circular do mitogenoma da arara-azul-de-lear, mostrando a posição e direção de transcrição de todos os genes. Genes codificadores na trilha externa, tRNAs e rRNAs na trilha intermediária, conteúdo GC na trilha interna. Utilizar `09_circular_map.svg` ou gerar versão de publicação no OGDRAW com o arquivo `A_leari.gbk`.

### 6.3 Genes Codificadores de Proteínas

Os 13 genes codificadores de proteínas (CDS — *Coding Sequences*, sequências codificantes) compreendem subunidades de quatro dos cinco complexos da **cadeia respiratória mitocondrial** — a sequência de cinco grandes complexos proteicos (numerados I a V) localizados na membrana interna da mitocôndria, responsáveis por converter a energia dos alimentos em ATP (a "moeda energética" da célula) por meio da fosforilação oxidativa: sete subunidades do complexo I — NADH desidrogenase (ND1, ND2, ND3, ND4, ND4L, ND5 e ND6); três subunidades do complexo IV — citocromo c oxidase (COX1, COX2 e COX3); duas subunidades do complexo V — ATP sintase (ATP6 e ATP8); e uma subunidade do complexo III — citocromo b (CYTB). A ausência de subunidades do complexo II (succinato desidrogenase) é característica conservada em todos os mitogenomas de metazoários.

Doze dos 13 CDS são transcritos a partir da fita pesada (H-strand); apenas o ND6 é codificado na fita leve (L-strand) — padrão universal em mitogenomas aviários. A tabela a seguir apresenta os códons de início e término de cada gene:

| Gene | Códon de Início | Códon de Término | Tamanho (bp) | Tamanho (aa) |
|------|----------------|-----------------|-------------|-------------|
| ND1 | ATG | TA(+A) | 974 | 324 |
| ND2 | ATA | TAA | 1.041 | 347 |
| COX1 | GTG | AGG | 1.548 | 516 |
| COX2 | ATG | TAA | 684 | 228 |
| ATP8 | ATG | TAA | 168 | 56 |
| ATP6 | ATG | TAA | 684 | 228 |
| COX3 | ATG | T(+AA) | 784 | 261 |
| ND4L | ATG | TAA | 297 | 99 |
| ND4 | ATG | T(+AA) | 1.378 | 459 |
| ND5 | ATA | TAA | 1.809 | 603 |
| CYTB | ATG | TAA | 1.140 | 380 |
| ND6 | ATG | TAG | 513 | 171 |
| ND3 | ATA | GAA | 349* | 118* |

*Tabela 7 — Códons de início e término dos 13 genes codificadores de proteínas.*

Observam-se padrões nos códons de início e término que são típicos de mitogenomas de vertebrados:

(i) O códon de início mais frequente é ATG (9 de 13 genes). Três genes utilizam ATA como códon de início alternativo (ND2, ND5, ND3), e o COX1 utiliza GTG — ambas são variantes aceitas no **código genético mitocondrial de vertebrados**. As mitocôndrias utilizam um código genético ligeiramente diferente do código universal empregado pelo núcleo celular: por exemplo, ATA codifica metionina (em vez de isoleucina), TGA codifica triptofano (em vez de ser um códon de parada) e AGA/AGG funcionam como códons de parada (em vez de codificarem arginina) (ANDERSON et al., 1981; BOORE, 1999). Essa variação é catalogada como Tabela 2 no sistema de codificação do NCBI.

(ii) Três genes (ND1, COX3 e ND4) terminam com códons de parada incompletos — consistindo apenas em "TA" ou "T" em vez do códon completo "TAA". Isso ocorre porque, no DNA mitocondrial, os genes são transcritos como um único RNA longo (**mRNA policistrônico** — isto é, um mRNA que contém múltiplos genes em sequência). Quando esse RNA é clivado em mRNAs individuais, a extremidade de cada mRNA recebe uma cauda de adeninas (A) por um processo chamado **poliadenilação** — a adição pós-transcricional de múltiplos nucleotídeos de adenina. Essa cauda completa o códon de parada: "T" + "AA" (da poliadenilação) = "TAA" funcional (ANDERSON et al., 1981).

(iii) O ND6 utiliza TAG como códon de término, menos frequente porém funcional.

**Seleção purificante e pressão evolutiva nos genes codificadores.** Uma análise recorrente em estudos de genômica mitocondrial é a razão Ka/Ks (também denominada dN/dS) — a proporção entre **substituições não-sinônimas** (Ka, mutações que alteram o aminoácido codificado) e **substituições sinônimas** (Ks, mutações silenciosas que não alteram o aminoácido). Quando Ka/Ks < 1, indica que mutações que alteram a proteína estão sendo eliminadas pela seleção natural (**seleção purificante**), evidenciando que a proteína desempenha função crítica e não tolera alterações; quando Ka/Ks > 1, indica **seleção positiva** (adaptação ativa a novas condições); quando Ka/Ks = 1, indica evolução neutra.

A literatura recente revela um padrão universal: **todos os 13 genes codificadores mitocondriais de todos os organismos analisados apresentam Ka/Ks < 1**, independentemente do táxon — aves, répteis, peixes, moluscos, insetos, crustáceos, cnidários (ZHAN et al., 2024; GUERREIRO et al., 2025; SHI et al., 2024; YE et al., 2024; KUNDU et al., 2024; ALBOASUD; JEONG; LEE, 2024; ZHOU et al., 2024; TAO et al., 2024). Isso confirma que os genes mitocondriais codificadores de proteínas estão sob forte constração funcional: suas proteínas são componentes essenciais da cadeia respiratória, e qualquer mutação deleteria é eliminada por seleção natural.

Consistentemente entre todos os estudos, o **COX1 é o gene mais conservado** (menor Ka/Ks, tipicamente 0,02–0,07) e o **ATP8 é o mais variável** (maior Ka/Ks, tipicamente 0,5–0,8). O COX1 codifica a subunidade principal do complexo IV (citocromo c oxidase), o último complexo da cadeia respiratória, cuja função de transferência de elétrons para o oxigênio é absolutamente crítica — razão pela qual é utilizado como marcador para código de barras de DNA. O ATP8, por outro lado, codifica a menor subunidade do complexo V (ATP sintase), com apenas 56 aminoácidos, e tolera mais variações sem perda de função. Notavelmente, em ramshorn snails (Planorbidae), o ATP8 chega a apresentar Ka/Ks > 1 em alguns gêneros, sugerindo seleção positiva (TAO et al., 2024) — fenômeno raro entre metazoários.

No mitogenoma da arara-azul-de-lear, a conservação do COX1 é evidenciada pela sua identidade com *A. hyacinthinus* e pela adoção universal do gene como marcador forense (cf. seção 6.6). A análise formal de Ka/Ks entre *A. leari* e *A. hyacinthinus*, embora não realizada neste trabalho, constitui uma perspectiva futura promissora: considerando a divergência de apenas 13 bp entre os dois mitogenomas, espera-se que virtualmente todas as substituições sejam sinônimas (Ks >> Ka), confirmando seleção purificante extrema no gênero *Anodorhynchus*.

**Uso de códons (RSCU).** Outra métrica informativa é o **RSCU** (*Relative Synonymous Codon Usage* — uso relativo de códons sinônmos), que mede a preferência por determinados códons entre os que codificam o mesmo aminoácido. Em mitogenomas, há um viés sistemático em favor de códons com A ou T na terceira posição, reflexo da composição nucleotídica rica em AT. Shi et al. (2024) demonstraram em afídeos que esse viés é dirigido predominantemente por **seleção natural** (e não apenas por pressão mutacional), conforme evidenciado por análises de ENC-plot, PR2-bias e neutralidade. Nos vertebrados, o padrão é menos extremo, mas a preferência por códons terminados em A ou T na terceira posição permanece detectável (KUNDU et al., 2024; GUERREIRO et al., 2025). A geração de tabelas de RSCU para o mitogenoma da arara-azul-de-lear constitui uma perspectiva futura integrável ao pipeline.

#### O Frameshift do Gene ND3: Uma Assinatura Evolutiva Conservada

O achado mais notável da anotação do mitogenoma da arara-azul-de-lear é a presença de um nucleotídeo extra no gene ND3, que interrompe o quadro de leitura (*reading frame*) em sua região central. Esse fenômeno, denominado *frameshift* do ND3, constitui uma das mais intrigantes particularidades da genômica mitocondrial de vertebrados e merece análise detalhada.

**Descoberta e descrição.** O frameshift do ND3 foi identificado pela primeira vez por Mindell, Sorenson e Dimcheff (1998), que sequenciaram o gene ND3 de 52 espécies de aves e 10 de tartarugas. Os autores observaram que, em uma posição conservada — após o nucleotídeo 174 contado a partir do códon de início —, há uma inserção de um único nucleotídeo que desloca o quadro de leitura em +1. Essa inserção divide o gene em dois quadros de leitura abertos (ORFs) consecutivos: o primeiro ORF codifica os 58 aminoácidos iniciais da proteína ND3 (posições 1–174 do mRNA), e o segundo ORF, no quadro deslocado, codifica os ~59 aminoácidos restantes. No mitogenoma da arara-azul-de-lear, o MITOS2 anotou automaticamente essas duas regiões como `nad3_0` (174 bp) e `nad3_1` (180 bp), com 5 bp de sobreposição na região do sítio de deslizamento — resultado plenamente consistente com o padrão descrito na literatura.

**Mecanismo molecular: duas hipóteses.** A produção de uma proteína ND3 funcional a partir de um gene com quadro de leitura interrompido exige um mecanismo de correção, e duas hipóteses principais foram propostas:

(i) *Deslizamento ribossomal programado (+1 PRF).* Neste modelo, o ribossomo mitocondrial, ao traduzir o mRNA do ND3, encontra o sítio de frameshift e "desliza" um nucleotídeo para frente (+1), continuando a tradução no novo quadro de leitura. Esse mecanismo, conhecido como *programmed +1 ribosomal frameshifting*, difere fundamentalmente do deslizamento −1 (−1 PRF) amplamente estudado em retrovírus como HIV e RSV (JACKS et al., 1988). No −1 PRF viral, o ribossomo retrocede um nucleotídeo em uma sequência escorregadia (*slippery sequence*) com motivo X_XXY_YYH, auxiliado por estruturas secundárias de RNA como pseudonós (*pseudoknots*) a jusante (BRIERLEY, 1995; ATKINS et al., 2016). No +1 PRF do ND3 mitocondrial, o mecanismo proposto envolve uma pausa ribossomal causada pela baixa disponibilidade de um tRNA raro correspondente ao códon no sítio A, favorecendo cineticamente o deslizamento para frente, onde um tRNA mais abundante pode parear com o códon no novo quadro (HARGER; MESKAUSKAS; DINMAN, 2002).

(ii) *Edição pós-transcricional do mRNA.* Neste modelo alternativo, o nucleotídeo extra é removido do mRNA por um mecanismo de edição de RNA antes da tradução, restaurando o quadro de leitura contínuo. Evidências de edição de RNA mitocondrial em plantas e protistas são abundantes, mas sua ocorrência em mitocôndrias de vertebrados é mais limitada. Estudos com cDNA de ND3 em tartarugas sugeriram que ao menos parte dos transcritos pode ser editada (RUSSELL; BECKENBACH, 2008), embora a questão permaneça em aberto. É possível, inclusive, que ambos os mecanismos coexistam: parte dos transcritos sendo editados e parte sendo traduzidos via frameshifting, em um sistema de redundância funcional.

> **[SUGESTÃO DE FIGURA]**: Diagrama esquemático do frameshift do ND3, mostrando: (a) a posição da inserção do nucleotídeo extra após a posição 174; (b) os dois ORFs resultantes (nad3_0 e nad3_1); (c) o mecanismo de +1 PRF (deslizamento ribossomal) versus edição de RNA. Comparar com o modelo de −1 PRF viral para evidenciar as diferenças.

**Distribuição taxonômica e significado evolutivo.** O frameshift do ND3 não é um artefato nem uma peculiaridade restrita a psitacídeos: trata-se de uma característica compartilhada por toda a classe Aves (todas as ordens investigadas), pela ordem Testudines (tartarugas e cágados) e pela ordem Crocodylia (crocodilos e jacarés) (MINDELL; SORENSON; DIMCHEFF, 1998). Russell e Beckenbach (2008) ampliaram significativamente esse panorama ao sequenciar o mitogenoma da tartaruga *Trachemys scripta*, revelando que as tartarugas exibem uma diversidade de sítios de frameshifting sem paralelo entre os vertebrados: foram documentadas inserções de frameshift em **seis sítios separados, em três genes distintos** (nad3, nad4l e cob), a maioria únicos de linhagens particulares. Além disso, os autores identificaram em *Chelodina parkeri* (tartaruga-australiana) a **redefinição de um códon de parada** — o primeiro exemplo documentado desse fenômeno em mitocôndrias de vertebrados —, demonstrando que os sistemas de tradução mitocondrial de quelônios toleram não apenas frameshifting, mas também recodificação de códons. O fato de que a maioria dos sítios de frameshift em tartarugas é específica de linhagens particulares e provavelmente recente sugere que o frameshifting não tem papel regulatório, mas é tolerado sob condições específicas que combinam uma sequência nucleotídica particular com um sistema de tradução amenável ao deslizamento ribossomal (RUSSELL; BECKENBACH, 2008).

Notavelmente, o frameshift está **ausente** em mamíferos, anfíbios e na maioria dos peixes. Essa distribuição filogenética mapeia-se precisamente ao clado **Archelosauria** (Aves + Crocodylia + Testudines), sugerindo que a inserção do nucleotídeo extra ocorreu uma única vez no ancestral comum desses grupos, há aproximadamente 250 milhões de anos, no início do Mesozoico — e foi mantida por seleção purificadora desde então. A conservação de um frameshift por um período tão extenso de tempo evolutivo, atravessando centenas de milhões de gerações sem ser eliminado por seleção natural, constitui uma forte evidência de que o mecanismo de correção (seja por frameshifting ribossomal, edição de RNA ou ambos) é altamente eficiente e funcionalmente indispensável.

**Resultado no mitogenoma da arara-azul-de-lear.** A detecção automática do frameshift pelo MITOS2, que anotou o ND3 como duas regiões separadas (`nad3_0` e `nad3_1`), confirma que o pipeline desenvolvido neste trabalho é capaz de identificar corretamente essa assinatura evolutiva sem intervenção manual. O resultado é plenamente consistente com o padrão verificado no congênere *A. hyacinthinus* (GenBank NC_082165.1), cuja anotação oficial registra `CDS join(9504..9677,9679..9854)` para o gene ND3. Na conversão para formato GenBank realizada pelo script `generate_genbank.py`, o ND3 foi representado com as coordenadas unidas pela anotação `join()` e o qualificador `/exception="ribosomal slippage"`, conforme as diretrizes do NCBI para genes com frameshift programado. Essa padronização é essencial para que a submissão seja aceita pelos bancos de dados internacionais e para que a anotação reflita corretamente a biologia subjacente.

### 6.4 RNAs Transportadores e Ribossomais

Os 22 genes de tRNA identificados codificam transportadores para os 20 aminoácidos padrão. A nomenclatura dos tRNAs segue o padrão "tRNA-Xxx", onde "Xxx" é a abreviatura de três letras do aminoácido que o tRNA carrega: por exemplo, tRNA-Phe transporta fenilalanina, tRNA-Val transporta valina, tRNA-Ile transporta isoleucina, e assim por diante. Dois aminoácidos — leucina e serina — possuem dois tRNAs cada, chamados **isoaceptores** (tRNAs diferentes que carregam o mesmo aminoácido, mas reconhecem códons distintos): tRNA-Leu(UUR) e tRNA-Leu(CUN) para leucina, e tRNA-Ser(UCN) e tRNA-Ser(AGY) para serina. As siglas entre parênteses indicam os códons reconhecidos por cada isoaceptor. Essa configuração com 22 tRNAs é universal em mitogenomas de metazoários. Quatorze tRNAs são transcritos a partir da fita pesada e oito a partir da fita leve.

Cada tRNA foi caracterizado quanto ao seu anticódon e aminoácido correspondente:

A tabela a seguir apresenta as características de cada tRNA, incluindo seu tamanho em pares de bases, o **anticódon** — a trinca de nucleotídeos no tRNA que se pareia com o códon correspondente no mRNA durante a tradução — e o aminoácido que transporta. Os anticódons são apresentados na direção 5'→3' da molécula de tRNA; por serem complementares, leem-se na direção oposta ao códon do mRNA.

| tRNA | Tamanho (bp) | Anticódon | Aminoácido |
|------|-------------|-----------|------------|
| tRNA-Phe | 67 | GAA | Fenilalanina |
| tRNA-Val | 69 | TAC | Valina |
| tRNA-Leu(UUR) | 74 | TAA | Leucina |
| tRNA-Ile | 71 | GAT | Isoleucina |
| tRNA-Gln | 70 | TTG | Glutamina |
| tRNA-Met | 68 | CAT | Metionina |
| tRNA-Trp | 70 | TCA | Triptofano |
| tRNA-Ala | 68 | TGC | Alanina |
| tRNA-Asn | 73 | GTT | Asparagina |
| tRNA-Cys | 66 | GCA | Cisteína |
| tRNA-Tyr | 70 | GTA | Tirosina |
| tRNA-Ser(UCN) | 75 | TGA | Serina |
| tRNA-Asp | 68 | GTC | Aspartato |
| tRNA-Lys | 68 | TTT | Lisina |
| tRNA-Gly | 68 | TCC | Glicina |
| tRNA-Arg | 69 | TCG | Arginina |
| tRNA-His | 68 | GTG | Histidina |
| tRNA-Ser(AGY) | 65 | GCT | Serina |
| tRNA-Leu(CUN) | 70 | TAG | Leucina |
| tRNA-Thr | 67 | TGT | Treonina |
| tRNA-Pro | 68 | TGG | Prolina |
| tRNA-Glu | 68 | TTC | Glutamato |

*Tabela 8 — Anticódons dos 22 genes de tRNA do mitogenoma da arara-azul-de-lear.*

Os tamanhos dos tRNAs variam de 65 bp (tRNA-Ser(AGY)) a 75 bp (tRNA-Ser(UCN)), com tamanho médio de 69 bp — dentro do intervalo típico de 60–80 bp observado em mitogenomas aviários. As estruturas secundárias dos 22 tRNAs foram preditas computacionalmente utilizando o programa RNAplot do pacote ViennaRNA (LORENZ et al., 2011), a partir dos modelos de covariância do MITOS2. Todos os tRNAs apresentaram a **estrutura canônica em trevo** (*cloverleaf*) — a conformação tridimensional característica dos tRNAs, em que a molécula se dobra sobre si mesma formando quatro "braços" que lembram uma folha de trevo: o braço aceptor (onde o aminoácido se liga), o braço do anticódon (que reconhece o códon do mRNA), o braço TΨC (envolvido na ligação ao ribossomo) e o **braço DHU** (nomeado pela presença de diidrouridina, um nucleotídeo modificado, e envolvido na estabilidade estrutural). O tRNA-Ser(AGY) constitui a exceção: carece do braço DHU, adotando uma estrutura simplificada em forma de "D". Essa anomalia estrutural não é uma peculiaridade da arara-azul-de-lear nem das aves em geral: trata-se de uma característica conservada em **virtualmente todos os metazoários** e considerada um relícto evolutivo ancestral mantido desde a origem dos animais multicelulares (BOORE, 1999). Nos 15 artigos analisados neste trabalho, abrangendo táxons tão diversos quanto corais hexacorálios (WEI et al., 2024), moluscos bivalves (LI et al., 2023; KARTAVTSEV; MASALKOVA, 2024), estrelas-do-mar (ALBOASUD; JEONG; LEE, 2024), afídeos (SHI et al., 2024), ramshorn snails (TAO et al., 2024), tubarões (YE et al., 2024), arraias (GUERREIRO et al., 2025), lagartos (ZHAN et al., 2024), anfíbios (HONG et al., 2024) e peixes teléosteos (KUNDU et al., 2024), **todos** reportam a ausência do braço DHU especificamente no tRNA-Ser(AGY). Essa conservação por mais de 600 milhões de anos, atravessando todas as grandes divisões do reino animal, indica que o tRNA-Ser(AGY) sem braço DHU não é disfuncional — pelo contrário, a estrutura simplificada em "D" é reconhecida pelo ribossomo mitocondrial de modo alternativo, possivelmente envolvendo interações compensatórias com outros componentes da maquinaria de tradução. A presença desse padrão no mitogenoma da arara-azul-de-lear, portanto, corrobora a integridade da anotação realizada pelo pipeline.

> **[SUGESTÃO DE FIGURA]**: Painel 4×6 com as estruturas secundárias em trevo dos 22 tRNAs + 2 rRNAs, geradas em formato SVG. Utilizar os arquivos da pasta `structure_svgs/tRNA/` e `structure_svgs/rRNA/`. Destacar o tRNA-Ser(AGY) como o único sem braço DHU.

Os dois genes de rRNA — 12S rRNA (rrnS, 970 bp) e 16S rRNA (rrnL, 1.572 bp) — estão localizados entre o tRNA-Phe e o tRNA-Leu(UUR), separados pelo tRNA-Val. Essa organização é inalterada em todos os Psittaciformes sequenciados até o momento. Os rRNAs mitocondriais compõem os ribossomos mitocondriais (55S, formados por subunidade 28S com 12S rRNA e subunidade 39S com 16S rRNA), responsáveis pela tradução dos 13 mRNAs mitocondriais em proteínas dentro da organela.

### 6.5 Organização Gênica e Comparação com *Anodorhynchus hyacinthinus*

A ordem gênica do mitogenoma da arara-azul-de-lear segue o arranjo padrão de aves neognatas (DESJARDINS; MORAIS, 1990):

```
D-loop — trnF — rrnS — trnV — rrnL — trnL2 — nad1 — trnI — trnQ(-) — trnM — nad2 —
trnW — trnA(-) — trnN(-) — trnC(-) — trnY(-) — cox1 — trnS2(-) — trnD — cox2 — trnK —
atp8 — atp6 — cox3 — trnG — nad3 — trnR — nad4l — nad4 — trnH — trnS1 — trnL1 —
nad5 — cob — trnT — trnP(-) — nad6(-) — trnE(-) — D-loop
```

(Os genes marcados com "(-)" são transcritos a partir da fita leve.)

A comparação direta com o mitogenoma de *A. hyacinthinus* (NC_082165.1, 16.999 bp) revela:

| Característica | Arara-azul-de-lear | Arara-azul-grande |
|---------------|-------------------|-------------------|
| Tamanho total | 16.986 bp | 16.999 bp |
| Diferença | - | Δ 13 bp |
| Genes codificadores | 13 | 13 |
| tRNAs | 22 | 22 |
| rRNAs | 2 | 2 |
| Ordem gênica | Padrão aviário | Padrão aviário |
| Frameshift ND3 | Presente | Presente |
| CDS na fita leve | ND6 | ND6 |

*Tabela 9 — Comparação entre os mitogenomas de A. leari e A. hyacinthinus.*

A diferença de apenas 13 bp entre os dois mitogenomas — aproximadamente 0,08% do comprimento total — sugere uma **divergência recente** entre as duas espécies, consistente com estimativas filogenéticas que posicionam a separação do gênero *Anodorhynchus* no Plioceno tardio (TAVARES et al., 2006). A conservação do frameshift do ND3, da ordem gênica e da composição geral confirmam que os dois mitogenomas são sintatênicos — isto é, compartilham a mesma organização estrutural, sem rearranjos cromossômicos detectáveis.

Essa elevada similaridade tem implicações práticas para a conservação. A literatura recente alerta para limitações do COX1 como marcador único de identificação: Ye et al. (2024), ao analisar tubarões do gênero *Scoliodon*, demonstraram que o COX1 sozinho é insuficiente para distinguir espécies crípticas muito próximas, sendo necessário recorrer a múltiplos genes (CYTB, ND1, ND2, ND5) ou ao mitogenoma completo para resolução taxonômica confiável. Esse achado é diretamente relevante para o gênero *Anodorhynchus*: com apenas 13 bp de divergência total, marcadores mitocondriais baseados em regiões conservadas (como COX1 e CYTB) provavelmente não terão resolução suficiente para distinguir espécimes das duas espécies em análises forenses.

**Região controle (D-loop) como marcador diferencial.** A região controle é o segmento do mitogenoma com maior taxa de variação interespecífica e intraespecífica. Estudos recentes confirmam esse padrão de forma consistente: em tubarões *Scoliodon*, a D-loop apresenta a maior variabilidade entre todas as regiões do mitogenoma (YE et al., 2024); em arraias *Potamotrygon*, a região controle contém repetições em tandem específicas da espécie (GUERREIRO et al., 2025); em afídeos, o comprimento da D-loop varia geograficamente dentro da mesma espécie, de 650 a 1.109 bp, refletindo diferenças no número de cópias de repetições em tandem (SHI et al., 2024). A comparação detalhada da região controle entre *A. leari* e *A. hyacinthinus* — incluindo a identificação de domínios conservados centrais (CCDs), blocos de sequência conservados (CSBs) e repetições em tandem específicas — constitui uma perspectiva futura estratégica para o desenvolvimento de marcadores diagnósticos de alta resolução para identificação forense das duas espécies. A disponibilidade de mitogenomas completos de ambas permite identificar regiões com maior variabilidade (como a região controle/D-loop) que podem servir como marcadores diagnósticos mais informativos.

> **[SUGESTÃO DE FIGURA — Diagrama de Sintenia]**: Diagrama de colinearidade comparando a ordem gênica dos mitogenomas de *A. leari* e *A. hyacinthinus*, com linhas conectoras gene-a-gene mostrando a correspondência total. Evidenciar a conservação completa da ordem gênica (sintatênia) entre as duas espécies. Formato similar ao gráfico de sintenia do PDF de referência da disciplina.

### 6.6 Implicações para a Conservação

A reconstrução do mitogenoma completo da arara-azul-de-lear contribui diretamente para três frentes da genômica aplicada à conservação. Para contextualizar essas aplicações, é necessário definir alguns conceitos da genética da conservação:

- **Haplótipo**: variante genética definida por um conjunto específico de mutações em uma região do DNA, compartilhada por um grupo de indivíduos que herdaram essa sequência de um ancestral comum. Em genética mitocondrial, cada haplótipo representa uma linhagem materna distinta. O número de haplótipos em uma população é uma medida direta de sua diversidade genética.
- **Diversidade nucleotídica** ($\pi$): a proporção média de nucleotídeos que diferem entre dois indivíduos da mesma população. Valores baixos indicam pouca variação genética, típica de populações que passaram por gargalos.
- **Depressão endogâmica** (*inbreeding depression*): redução da aptidão biológica (sobrevivência, fertilidade, resistência a doenças) causada pelo acúmulo de **alelos** (variantes alternativas de um mesmo gene) deletérios recessivos quando indivíduos aparentados se reproduzem entre si — situação inevitável em populações muito pequenas.
- **Filogeografia**: disciplina que estuda a distribuição geográfica de linhagens genéticas (haplótipos) dentro de uma espécie ou grupo de espécies, buscando reconstruir os processos históricos (glaciações, fragmentação de habitat, dispersão) que moldaram essa distribuição.
- **Relógio molecular**: método que estima tempos de divergência evolutiva entre espécies com base na taxa de acumulação de mutações. Partindo do pressuposto de que mutações se acumulam a uma taxa aproximadamente constante ao longo do tempo, a diferença genética entre duas sequências pode ser convertida em uma estimativa de tempo desde o ancestral comum.
- **Especiação**: processo pelo qual uma população ancestral se divide em duas ou mais espécies distintas, tipicamente por isolamento geográfico ou ecológico prolongado.

**Identificação forense de espécimes.** O tráfico ilegal continua sendo uma das maiores ameaças à arara-azul-de-lear. Um paralelo relevante é o trabalho de Guerreiro et al. (2025) com a arraia *Potamotrygon leopoldi*, espécie endêmica do rio Xingu (Amazônia brasileira) igualmente ameaçada. Os autores sequenciaram o mitogenoma completo (17.504 bp) e identificaram sítios sob seleção positiva em genes como ND1, ND4, ND5, COX1 e COX2, vinculados à adaptação ao ambiente de água doce. Esse estudo demonstra como o mitogenoma completo permite investigar não apenas identidade taxonômica, mas também processos adaptativos — uma possibilidade que se abre para a arara-azul-de-lear em relação às pressões ambientais do semiárido baiano.

A disponibilidade de uma referência genômica mitocondrial completa permite o desenvolvimento de protocolos de **código de barras de DNA** (*DNA barcoding*) — técnica que utiliza uma sequência curta e padronizada de DNA (tipicamente o gene COX1) como identificador taxonômico, de forma análoga a um código de barras comercial (HEBERT et al., 2003). O conceito de DNA barcoding surgiu como resposta ao progressivo declinio da expertise taxonômica tradicional: conforme argumentam Hebert et al. (2003), a capacidade de identificação biológica baseada em morfologia está em colapso, e a única perspectiva sustentável para uma capacidade de identificação em larga escala reside na construção de sistemas que empreguem sequências de DNA como "códigos de barras" de táxons, sendo o gene COI (citocromo c oxidase I) o núcleo de um sistema global de bioidentificação para animais. Os marcadores COX1 (1.548 bp) e CYTB (1.140 bp) podem ser utilizados como primeira abordagem por laboratórios forenses para confirmar a identidade taxonômica de espécimes apreendidos. Entretanto, conforme demonstrado por Ye et al. (2024) em tubarões *Scoliodon*, o COX1 isoladamente pode ser insuficiente para distinguir espécies muito próximas — limitacao que reforça a importância de dispor do mitogenoma completo como referência. Os mesmos autores demonstraram que modelos de aprendizado de máquina (Random Forest, SVM, MLP) treinados com múltiplos genes mitocondriais atingem AUC > 0,90 para identificação de espécies crípticas — estratégia que poderia ser adaptada para distinguir espécimes de *A. leari* e *A. hyacinthinus* em contexto forense. Com o mitogenoma completo de ambas as espécies agora disponível, esses marcadores multigenênicos podem ser utilizados — ovos, penas, tecidos ou indivíduos em cativeiro —, subsidiando processos judiciais e ações de fiscalização. A diferenciação em relação a *A. hyacinthinus*, cuja morfologia é semelhante, é particularmente relevante dado que as duas espécies possuem status de conservação e regulamentações distintas.

**Avaliação da diversidade genética.** A população da arara-azul-de-lear sofreu um **gargalo populacional** severo nas décadas de 1970-1990 — isto é, uma redução drástica e súbita no número de indivíduos, que elimina grande parte da variação genética e cujos efeitos persistem por muitas gerações —, quando chegou a cerca de 200 indivíduos. Gargalos dessa magnitude tipicamente reduzem a variabilidade genética, aumentando a vulnerabilidade a doenças e a depressão endogâmica (cf. definições acima). A sequência mitocondrial completa fornece a base para estudos de diversidade genética mitocondrial: ao sequenciar a região controle (D-loop) de múltiplos indivíduos, é possível estimar o número de haplótipos circulantes, a diversidade nucleotídica e o grau de estruturação genética (diferenciação entre subpopulações que trocam poucos migrantes) entre os núcleos populacionais de Toca Velha e Serra Branca. Esses dados são insumos diretos para o PAN Arara-azul-de-lear (ICMBio, 2022) e para decisões de manejo como translocações e programas de reprodução em cativeiro.

**Filogeografia e evolução do gênero.** O mitogenoma completo permite **análises filogenéticas** robustas — reconstruções das relações de parentesco evolutivo entre espécies, representadas em diagramas ramificados chamados árvores filogenéticas —, baseadas em múltiplos genes simultaneamente, para reconstruir as relações evolutivas dentro de Psittaciformes. A divergência de apenas 13 bp em relação a *A. hyacinthinus* pode ser utilizada para calibrar relógios moleculares (cf. definições acima) e estimar tempos de divergência, contribuindo para a compreensão da biogeografia histórica — o estudo de como eventos geológicos e climáticos do passado moldaram a distribuição atual das espécies — do gênero *Anodorhynchus* e dos processos de especiação na Caatinga.

### 6.7 Disponibilidade dos Dados

Os dados produzidos neste trabalho estão disponíveis nos seguintes formatos e repositórios:

- **GenBank Flat File** (`A_leari.gbk`): arquivo no formato padrão do NCBI, pronto para submissão, contendo a sequência completa anotada com todos os 37 genes, códons de início e término, qualificadores de produto e o tratamento correto do frameshift do ND3 (anotação `join()` com `/exception="ribosomal slippage"`).
- **Feature Table** (`A_leari.tbl`) e **FASTA formatado** (`A_leari.fsa`): arquivos complementares para processamento via `table2asn` do NCBI.
- **Repositório GitHub**: todo o código-fonte do pipeline, incluindo Dockerfiles, módulos Nextflow, scripts de análise e configurações de perfis, está disponível publicamente, permitindo a reprodução integral das análises e a aplicação a outras espécies (WILKINSON et al., 2016).

A disponibilidade pública contribui para a democratização do acesso à bioinformática, permitindo que outros grupos de pesquisa reproduzam as análises ou apliquem o pipeline a outras espécies — aves, peixes, invertebrados ou qualquer metazoário — necessitando apenas de um acesso SRA e uma sequência-semente de gene mitocondrial. A generalização do pipeline foi verificada pela execução bem-sucedida em dois táxons distintos (Psittaciformes e Perciformes), reforçando a aderência aos princípios FAIR (*Findable, Accessible, Interoperable, and Reusable*).

---

## Capítulo 7 — Considerações Finais

A Figura 24 sintetiza o percurso completo do trabalho desenvolvido, desde os dados brutos até os entregáveis finais.

> **Figura 24** — Infográfico-resumo do trabalho: do problema inicial (87 GB de dados brutos), passando pelas 7 etapas do pipeline, até o resultado final (16.986 bp circularizado, 37 genes, GenBank-ready) e as 6 contribuições originais. A validação cruzada com *A. hyacinthinus* (99,888% de identidade com a referência NCBI) reforça a confiabilidade da abordagem.

Conforme representado na Figura 24, o pipeline transforma 87 GB de dados brutos de sequenciamento em um mitogenoma anotado de 16.986 bp contendo 37 genes, por meio de sete etapas automatizadas. O infográfico destaca as seis contribuições originais deste trabalho: a integração completa do fluxo analítico, a containerização com Docker, a orquestração via Nextflow, o módulo Pilot QC para redução automática de dados, a geração de entregáveis em formato GenBank e a viabilidade de execução em hardware doméstico.

O presente trabalho teve como objetivo propor e implementar um pipeline de bioinformática voltado à montagem e anotação funcional de genomas mitocondriais, concebido segundo princípios de automação, reprodutibilidade, ciência aberta e acessibilidade computacional. A integração de sete ferramentas consolidadas — SRA-Toolkit, FastQC, Trim Galore com Cutadapt, NOVOPlasty, MITOS2, ViennaRNA e Biopython (COCK et al., 2009) — em um fluxo containerizado (Docker) e orquestrado por meio do Nextflow constitui a principal contribuição metodológica desta proposta, ao oferecer uma solução transparente, portável e replicável para análises genômicas.

Um resultado particularmente significativo deste trabalho é a demonstração de que a montagem e anotação completa de mitogenomas a partir de dados públicos é viável em hardware doméstico. O dataset original utilizado totaliza 87 GB de leituras brutas — volume que, sem as estratégias de redução implementadas, exigiria centenas de gigabytes de armazenamento temporário e servidores dedicados. A cadeia de reduções progressivas incorporada ao pipeline (Pilot QC, extração seletiva com `fastq-dump -X`, *trimming*) reduz esse volume em aproximadamente 5 milhões de vezes, viabilizando a execução em notebooks com 8 GB de RAM em pouco mais de 30 minutos. Essa acessibilidade é especialmente relevante para estudantes, pesquisadores iniciantes e instituições com recursos computacionais limitados, que anteriormente dependeriam de acesso a servidores ou *clusters* para realizar análises genômicas equivalentes.

Do ponto de vista biológico, o pipeline foi aplicado com sucesso a três espécies: a arara-azul-de-lear (*A. leari*, 16.986 bp circularizado, 176× de cobertura, 37 genes), a espécie de validação *Diploprion bifasciatus* (circularizado) e a arara-azul-grande (*A. hyacinthinus*, 16.999 bp circularizado a partir de dados metagenômicos). Esta última execução constitui a validação mais rigorosa do pipeline: a montagem *de novo* reproduziu a referência NCBI (OR209186.1) com 99,888% de identidade nucleotídica, com as 19 diferenças concentradas na região hipervariável (D-loop) e consistentes com polimorfismo intraespecífico. A reciprocidade de sementes — *A. leari* montada com semente de *A. hyacinthinus* e vice-versa — valida a estratégia de uso de sequências congêneres quando a espécie-alvo não possui referência própria. Conforme detalhado no Capítulo 6, a caracterização completa do mitogenoma da arara-azul-de-lear revelou composição nucleotídica, organização gênica e assinaturas moleculares (como o *frameshift* do ND3) consistentes com a espécie congênere *A. hyacinthinus*, confirmando a qualidade da montagem e a aplicabilidade dos dados para estudos de conservação, identificação forense e filogeografia.

Dentre as contribuições originais do pipeline, destacam-se: (i) o módulo Pilot QC, que determina automaticamente o volume de dados necessário para atingir a cobertura desejada, evitando processamento desnecessário e viabilizando a execução em computadores com armazenamento limitado; (ii) a cadeia de reduções progressivas (Pilot QC → `fastq-dump -X` → *trimming*) que transforma 87 GB de dados brutos em um mitogenoma de 17 KB, eliminando a necessidade de servidores dedicados; (iii) a lógica iterativa de montagem com re-seeding automático e múltiplos k-mers; (iv) a geração automática de estruturas secundárias de tRNAs e rRNAs via RNAplot/ViennaRNA, funcionalidade anteriormente restrita à interface web do MITOS2; (v) a conversão automática para GenBank Flat File com tratamento do frameshift do ND3; e (vi) a compilação automática de 14 categorias de entregáveis em pasta organizada.

Apesar de seu potencial, algumas limitações precisam ser reconhecidas. O sucesso da montagem depende da qualidade da seed inicial e da profundidade de cobertura dos dados de sequenciamento. A anotação do MITOS2, embora automatizada, pode requerer curadoria manual em casos de estruturas gênicas atípicas (como o frameshift do ND3). Adicionalmente, o tamanho dos containers Docker (~2–4 GB por imagem) pode representar uma barreira em ambientes com espaço de armazenamento limitado.

Como perspectivas futuras, propõem-se: (i) a submissão do mitogenoma da arara-azul-de-lear ao GenBank; (ii) a integração de análises filogenéticas automáticas ao pipeline (MAFFT + RAxML/IQ-TREE); (iii) a geração automática de tabelas de RSCU (*Relative Synonymous Codon Usage*) e de razão Ka/Ks para análise de pressão seletiva sobre os genes codificadores — análises que se tornaram padrão em estudos recentes de genômica mitocondrial (ZHAN et al., 2024; KUNDU et al., 2024; SHI et al., 2024); (iv) a integração de tRNAscan-SE como módulo opcional de verificação dos tRNAs preditos pelo MITOS2, seguindo a prática metodológica adotada em estudos recentes (ZHAN et al., 2024; KUNDU et al., 2024; SHI et al., 2024); (v) a geração de gráficos de sintenia comparativa entre espécies congêneres; e (vi) a publicação do pipeline no nf-core como recurso comunitário.

Em síntese, o pipeline desenvolvido representa uma contribuição metodológica significativa para a bioinformática de organelas, ao unir rigor técnico, aplicabilidade biológica e boas práticas de reprodutibilidade computacional. Mais do que atender a um caso de estudo específico, constitui um recurso com potencial de ser reutilizado, adaptado e expandido por diferentes grupos de pesquisa, colaborando para a consolidação de uma ciência mais transparente, padronizada e sustentável.

---

## 8. Referências

ALBERTS, B. et al. Molecular Biology of the Cell. 7. ed. New York: W.W. Norton, 2022.

ANDERSON, S. et al. Sequence and organization of the human mitochondrial genome. Nature, v. 290, p. 457-465, 1981.

ANDREWS, S. FastQC: a quality control tool for high throughput sequence data. 2010. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Acesso em: 31 ago. 2025.

ASAKAWA, S.; HIMENO, H.; MIURA, K.; WATANABE, K. Nucleotide sequence and gene organization of the starfish Asterina pectinifera mitochondrial genome. Genetics, v. 140, n. 3, p. 1047-1060, 1995.

ATKINS, J. F. et al. Ribosomal frameshifting and transcriptional slippage: from genetic steganography and cryptography to adventitious use. Nucleic Acids Research, v. 44, n. 15, p. 7007-7078, 2016. DOI: 10.1093/nar/gkw530.

BABRAHAM BIOINFORMATICS. Trim Galore! 2019. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. Acesso em: 31 ago. 2025.

BAKER, M. 1,500 scientists lift the lid on reproducibility. Nature, v. 533, n. 7604, p. 452-454, 2016.

BOORE, J. L. Animal mitochondrial genomes. Nucleic Acids Research, v. 27, n. 8, p. 1767-1780, 1999.

BRIERLEY, I. Ribosomal frameshifting on viral RNAs. Journal of General Virology, v. 76, pt. 8, p. 1885-1892, 1995. DOI: 10.1099/0022-1317-76-8-1885.

BRITISH ECOLOGICAL SOCIETY. A guide to reproducible code in ecology and evolution. London: British Ecological Society, 2017.

CALABRESE, C. et al. MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing. Bioinformatics, v. 30, n. 21, p. 3115-3117, 2014. DOI: 10.1093/bioinformatics/btu483.

COCK, P. J. A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, v. 25, n. 11, p. 1422-1423, 2009. DOI: 10.1093/bioinformatics/btp163.

COLLAR, N. J. Family Psittacidae (parrots). In: DEL HOYO, J.; ELLIOTT, A.; SARGATAL, J. (eds.). Handbook of the Birds of the World. v. 4. Barcelona: Lynx Edicions, 1997. p. 280-477.

DESJARDINS, P.; MORAIS, R. Sequence and gene organization of the chicken mitochondrial genome: a novel gene order in higher vertebrates. Journal of Molecular Biology, v. 212, n. 4, p. 599-634, 1990. DOI: 10.1016/0022-2836(90)90225-B.

DIERCKXSENS, N.; MARDULYN, P.; SMITS, G. NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, v. 45, n. 4, p. e18, 2017.

DI TOMMASO, P. et al. Nextflow enables reproducible computational workflows. Nature Biotechnology, v. 35, n. 4, p. 316-319, 2017.

DONATH, A. et al. Improved annotation of protein-coding gene boundaries in metazoan mitochondrial genomes. Molecular Ecology Resources, v. 19, n. 4, p. 609-615, 2019. DOI: 10.1111/1755-0998.12985.

EKBLOM, R.; WOLF, J. B. W. A field guide to whole-genome sequencing, assembly and annotation. Evolutionary Applications, v. 7, n. 9, p. 1026-1042, 2014.

EWING, B. et al. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research, v. 8, n. 3, p. 175-185, 1998. DOI: 10.1101/gr.8.3.175.

EWELS, P. A. et al. The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, v. 38, n. 3, p. 276-278, 2020.

GRÜNING, B. et al. Practical computational reproducibility in the life sciences. Cell Systems, v. 6, n. 6, p. 631-635, 2018.

HAHN, C.; BACHMANN, L.; CHEVREUX, B. Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads — a baiting and iterative mapping approach. Nucleic Acids Research, v. 41, n. 13, p. e129, 2013. DOI: 10.1093/nar/gkt371.

HARGER, J. W.; MESKAUSKAS, A.; DINMAN, J. D. An "integrated model" of programmed ribosomal frameshifting. Trends in Biochemical Sciences, v. 27, n. 9, p. 448-454, 2002. DOI: 10.1016/S0968-0004(02)02149-7.

HEBERT, P. D. N. et al. Biological identifications through DNA barcodes. Proceedings of the Royal Society of London B, v. 270, n. 1512, p. 313-321, 2003. DOI: 10.1098/rspb.2002.2218.

ICMBio — INSTITUTO CHICO MENDES DE CONSERVAÇÃO DA BIODIVERSIDADE. Plano de Ação Nacional para a Conservação da Arara-azul-de-lear. Brasília: ICMBio, 2022.

IRIDIAN GENOMES. Anodorhynchus leari whole genome sequencing. NCBI Sequence Read Archive, acesso SRR28399504, 2024. Disponível em: https://www.ncbi.nlm.nih.gov/sra/SRR28399504. Acesso em: 10 mar. 2026.

JACKS, T. et al. Characterization of ribosomal frameshifting in HIV-1 gag-pol expression. Nature, v. 331, n. 6153, p. 280-283, 1988. DOI: 10.1038/331280a0.

JIN, J.-J. et al. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, v. 21, n. 1, p. 241, 2020. DOI: 10.1186/s13059-020-02154-5.

KÖSTER, J.; RAHMANN, S. Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, v. 28, n. 19, p. 2520-2522, 2012.

LARMAN, C.; BASILI, V. R. Iterative and incremental development: a brief history. IEEE Computer, v. 36, n. 6, p. 47-56, 2003. DOI: 10.1109/MC.2003.1204375.

LEINONEN, R.; SUGAWARA, H.; SHUMWAY, M. The Sequence Read Archive. Nucleic Acids Research, v. 39, supl. 1, p. D19-D21, 2011. DOI: 10.1093/nar/gkq1019.

LORENZ, R. et al. ViennaRNA Package 2.0. Algorithms for Molecular Biology, v. 6, art. 26, 2011. DOI: 10.1186/1748-7188-6-26.

MARTIN, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, v. 17, n. 1, p. 10-12, 2011. DOI: 10.14806/ej.17.1.200.

MERKEL, D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal, n. 239, p. 2, 2014.

MILLER, J. R.; KOREN, S.; SUTTON, G. Assembly algorithms for next-generation sequencing data. Genomics, v. 95, n. 6, p. 315-327, 2010.

MINDELL, D. P.; SORENSON, M. D.; DIMCHEFF, D. E. An extra nucleotide is not translated in mitochondrial ND3 of some birds and turtles. Molecular Biology and Evolution, v. 15, n. 11, p. 1568-1571, 1998.

NCBI — NATIONAL CENTER FOR BIOTECHNOLOGY INFORMATION. SRA Toolkit Documentation. Bethesda: National Library of Medicine, 2023. Disponível em: https://github.com/ncbi/sra-tools. Acesso em: 10 mar. 2026.

PENG, R. D. Reproducible research in computational science. Science, v. 334, n. 6060, p. 1226-1227, 2011. DOI: 10.1126/science.1213847.

RUSSELL, R. D.; BECKENBACH, A. T. Recoding of translation in turtle mitochondrial genomes: programmed frameshift mutations and evidence of a modified genetic code. Journal of Molecular Evolution, v. 67, n. 6, p. 682-695, 2008. DOI: 10.1007/s00239-008-9179-0.

SANDVE, G. K. et al. Ten simple rules for reproducible computational research. PLoS Computational Biology, v. 9, n. 10, p. e1003285, 2013.

SICK, H. Discovery of the home of the Indigo Macaw in Brazil. American Birds, v. 33, n. 3, p. 235-236, 1979.

SOMMERVILLE, I. Software Engineering. 10. ed. Harlow: Pearson Education, 2016.

TAVARES, E. S. et al. Phylogenetic relationships and historical biogeography of Neotropical parrots (Psittaciformes: Psittacidae: Arini) inferred from mitochondrial and nuclear DNA sequences. Systematic Biology, v. 55, n. 3, p. 454-470, 2006. DOI: 10.1080/10635150600697390.

ULIANO-SILVA, M. et al. MitoHiFi: a pipeline for mitochondrial genome assembly from PacBio HiFi reads. BMC Bioinformatics, v. 24, art. 288, 2023. DOI: 10.1186/s12859-023-05385-y.

WILKINSON, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data, v. 3, p. 160018, 2016.

ALBOASUD, M.; JEONG, H.; LEE, T. Complete Mitochondrial Genomes and Phylogenetic Analysis of Genus *Henricia* (Asteroidea: Spinulosida: Echinasteridae). *International Journal of Molecular Sciences*, v. 25, art. 5575, 2024. DOI: 10.3390/ijms25115575.

GUERREIRO, S. L. M. et al. Analysis of the Entire Mitogenome of the Threatened Freshwater Stingray *Potamotrygon leopoldi* (Myliobatiformes: Potamotrygonidae) and Comprehensive Phylogenetic Assessment in the Xingu River, Brazilian Amazon. *International Journal of Molecular Sciences*, v. 26, art. 8252, 2025. DOI: 10.3390/ijms26178252.

HONG, Y.-H. et al. Differential Mitochondrial Genome Expression of Four Hylid Frog Species under Low-Temperature Stress and Its Relationship with Amphibian Temperature Adaptation. *International Journal of Molecular Sciences*, v. 25, art. 5967, 2024. DOI: 10.3390/ijms25115967.

KARTAVTSEV, Y. P.; MASALKOVA, N. A. Structure, Evolution, and Mitochondrial Genome Analysis of Mussel Species (Bivalvia, Mytilidae). *International Journal of Molecular Sciences*, v. 25, art. 6902, 2024. DOI: 10.3390/ijms25136902.

KUNDU, S. et al. Mitogenomic Characterization and Phylogenetic Placement of African Hind, *Cephalopholis taeniops*: Shedding Light on the Evolution of Groupers (Serranidae: Epinephelinae). *International Journal of Molecular Sciences*, v. 25, art. 1822, 2024. DOI: 10.3390/ijms25031822.

LI, F. et al. The Complete Mitochondrial Genomes of Two Rock Scallops (Bivalvia: Spondylidae) Indicate Extensive Gene Rearrangements and Adaptive Evolution Compared with Pectinidae. *International Journal of Molecular Sciences*, v. 24, art. 13844, 2023. DOI: 10.3390/ijms241813844.

SHI, A. et al. Characterization, Codon Usage Pattern and Phylogenetic Implications of the Waterlily Aphid *Rhopalosiphum nymphaeae* (Hemiptera: Aphididae) Mitochondrial Genome. *International Journal of Molecular Sciences*, v. 25, art. 11336, 2024. DOI: 10.3390/ijms252111336.

TAO, K. et al. Comparative Mitogenome Analyses of Fifteen Ramshorn Snails and Insights into the Phylogeny of Planorbidae (Gastropoda: Hygrophila). *International Journal of Molecular Sciences*, v. 25, art. 2279, 2024. DOI: 10.3390/ijms25042279.

WEI, Z. et al. The Mitogenomic Landscape of Hexacorallia Corals: Insight into Their Slow Evolution. *International Journal of Molecular Sciences*, v. 25, art. 8218, 2024. DOI: 10.3390/ijms25158218.

YE, P. et al. Potential Cryptic Diversity in the Genus *Scoliodon* (Carcharhiniformes: Carcharhinidae): Insights from Mitochondrial Genome Sequencing. *International Journal of Molecular Sciences*, v. 25, art. 11851, 2024. DOI: 10.3390/ijms252111851.

ZHAN, L. et al. The Phylogenetic Relationships of Major Lizard Families Using Mitochondrial Genomes and Selection Pressure Analyses in Anguimorpha. *International Journal of Molecular Sciences*, v. 25, art. 8464, 2024. DOI: 10.3390/ijms25158464.

ZHOU, Z. et al. The Complete Mitochondrial Genome and Phylogenetic Analysis of the Freshwater Shellfish *Novaculina chinensis* (Bivalvia: Pharidae). *International Journal of Molecular Sciences*, v. 25, art. 67, 2024. DOI: 10.3390/ijms25010067.

WILSON, G. et al. Best practices for scientific computing. PLoS Biology, v. 12, n. 1, e1001745, 2014. DOI: 10.1371/journal.pbio.1001745.
