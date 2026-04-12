# TCC — Desenvolvimento e Implementação de um Pipeline de Bioinformática para Montagem e Análise de Mitogenomas

**Autor:** Matheus Sobreira Benevides
**Curso:** Bacharelado em Ciência da Computação — UERN
**Orientador:** Prof. Wilfredo
**Prof. TCC:** Prof. Carlos

---

## Capítulo 1 — Introdução

### 1.1 Introdução

O estudo de genomas mitocondriais ocupa uma posição central na biologia evolutiva, genética de populações e sistemática, uma vez que essas moléculas oferecem informações altamente informativas sobre parentesco, variação intra e interespecífica e processos de conservação. O DNA mitocondrial, por suas características estruturais e funcionais — como a herança predominantemente materna, a ausência de recombinação e a elevada taxa evolutiva em comparação ao DNA nuclear —, tornou-se um marcador molecular amplamente utilizado em diferentes contextos biológicos (BOORE, 1999; EKBLOM; WOLF, 2014).

O avanço das tecnologias de sequenciamento de nova geração (*Next-Generation Sequencing* — NGS) ampliou o acesso a grandes volumes de dados, favorecendo a obtenção de mitogenomas completos e de alta qualidade. No entanto, a reconstrução dessas moléculas a partir de dados brutos permanece um desafio técnico, especialmente no caso de leituras curtas, que exigem algoritmos sofisticados para superar problemas como regiões repetitivas e possíveis falhas de montagem. Nesse cenário, ferramentas especializadas, como o NOVOPlasty, têm desempenhado papel fundamental ao implementar estratégias de *seed-and-extend* capazes de recuperar organelas de forma eficiente a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017).

Embora a diversidade de softwares tenha contribuído para avanços significativos na montagem de organelas, a ausência de padronização metodológica ainda compromete a reprodutibilidade das análises. A multiplicidade de versões, dependências e configurações de ambiente gera obstáculos práticos para replicar experimentos, fenômeno amplamente conhecido como *dependency hell* (GRÜNING et al., 2018). Esse problema se insere no contexto mais amplo da crise de reprodutibilidade científica (BAKER, 2016), que tem motivado a adoção, na bioinformática, de práticas da engenharia de software, como a conteinerização de ambientes computacionais com Docker (MERKEL, 2014) e a orquestração de fluxos com sistemas de gerenciamento de workflows, a exemplo do Nextflow (DI TOMMASO et al., 2017).

Além da montagem propriamente dita, a anotação funcional dos genes mitocondriais constitui uma etapa igualmente relevante, que permite identificar os genes codificadores de proteínas, RNAs transportadores e RNAs ribossômicos presentes na molécula reconstruída. Ferramentas como o MITOS2 (DONATH et al., 2019) automatizam essa tarefa, mas sua execução — normalmente restrita a interfaces web — permanece desvinculada dos pipelines de montagem, criando uma descontinuidade no fluxo analítico.

Dessa forma, a presente pesquisa propõe-se a contribuir para esse campo ao desenvolver um pipeline automatizado, reprodutível e portável para a montagem e anotação de mitogenomas. A integração de ferramentas bioinformáticas consolidadas em um ecossistema containerizado e orquestrado garante não apenas consistência metodológica, mas também acessibilidade e transparência, em consonância com os princípios de ciência aberta e os critérios FAIR (*Findable, Accessible, Interoperable, Reusable*) (WILKINSON et al., 2016).

### 1.2 Justificativa

A análise de genomas mitocondriais (mitogenomas) representa um componente essencial para a biologia evolutiva, a genética de populações e a sistemática, dado que essas moléculas fornecem informações altamente informativas sobre relações de parentesco, variação intra e interespecífica e processos de conservação. Sua relevância científica foi consolidada em trabalhos seminais que destacaram o papel estrutural e evolutivo dos mitogenomas como marcadores moleculares (BOORE, 1999). Em seguida, estudos enfatizaram sua aplicabilidade em investigações de genética de populações e evolução molecular (EKBLOM; WOLF, 2014), enquanto contribuições mais recentes demonstraram seu potencial em análises de biodiversidade e conservação por meio de avanços metodológicos em montagem e anotação de organelas (ULIANO-SILVA et al., 2023).

Embora a obtenção de mitogenomas completos tenha sido amplamente favorecida pelo advento do Sequenciamento de Nova Geração (NGS), a execução de pipelines de montagem ainda enfrenta desafios técnicos significativos. A diversidade de ferramentas, versões e dependências envolvidas gera um cenário de difícil replicação, conhecido como *dependency hell*, que compromete a reprodutibilidade e a padronização das análises (SANDVE et al., 2013; BAKER, 2016; GRÜNING et al., 2018).

Esse problema se torna especialmente relevante quando aplicado ao estudo de espécies ameaçadas de extinção, nas quais a disponibilidade de amostras biológicas é limitada e a confiabilidade dos resultados genômicos é crucial para fundamentar decisões de conservação. É o caso da arara-azul-de-lear (*Anodorhynchus leari*), ave endêmica da Caatinga nordestina do Brasil, classificada como "Em Perigo" pela IUCN Red List. Com população estimada em pouco mais de 1.700 indivíduos na natureza (ICMBio, 2022), a espécie depende de dados genômicos acessíveis e validáveis para subsidiar programas de manejo, reprodução em cativeiro e monitoramento populacional. Embora seu mitogenoma já tenha sido previamente reconstruído em atividade acadêmica, essa montagem foi realizada de forma manual, sem padronização metodológica e sem disponibilização pública do fluxo de trabalho — o que impede sua verificação e replicação por outros pesquisadores.

Um desafio adicional reside na completude dos pipelines disponíveis. A maioria das ferramentas e workflows publicados concentra-se exclusivamente na etapa de montagem, sem integrar a anotação funcional nem a preparação dos resultados para submissão a bancos de dados públicos como o GenBank. Essa lacuna obriga o pesquisador a recorrer a interfaces web externas (e.g., MITOS2 via UseGalaxy) e a conversões manuais de formato, introduzindo etapas não rastreáveis que comprometem a reprodutibilidade do ciclo completo de análise.

Nesse contexto, torna-se necessário o desenvolvimento de soluções que conciliem robustez biológica e consistência computacional. A adoção de práticas modernas da engenharia de software, como a containerização e os sistemas de gerenciamento de workflows, tem se mostrado fundamental para superar limitações de reprodutibilidade, além de permitir que fluxos analíticos sejam portáveis, escaláveis e transparentes (DI TOMMASO et al., 2017; EWELS et al., 2020).

Assim, este trabalho justifica-se pela proposta de construir um pipeline automatizado e reprodutível que cubra o ciclo completo — desde a aquisição dos dados brutos até a geração de arquivos prontos para submissão ao GenBank —, unindo boas práticas de ciência aberta e engenharia computacional à relevância biológica dos mitogenomas. A aplicação ao caso concreto de *A. leari* demonstra não apenas a viabilidade técnica da abordagem, mas também sua contribuição direta para a conservação de uma espécie ameaçada da biodiversidade brasileira.

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

g) Validar o pipeline automatizado por meio de estudos de caso com dados públicos de sequenciamento de pelo menos duas espécies de grupos taxonômicos distintos, comparando os genomas montados e anotados com resultados de referência.

h) Disponibilizar o pipeline completo e sua documentação em repositório público no GitHub, assegurando sua acessibilidade e reprodutibilidade pela comunidade científica.

Dessa forma, este capítulo estabeleceu o contexto biológico e computacional do trabalho, destacando a relevância dos genomas mitocondriais, os desafios de reprodutibilidade na bioinformática e a proposta de um pipeline automatizado e portável. O Capítulo 2 apresenta, a seguir, os fundamentos teóricos necessários para compreender em maior profundidade os conceitos que sustentam a execução do pipeline.

---

## Capítulo 2 — Referencial Teórico

O presente capítulo reúne os fundamentos conceituais necessários para compreender a proposta metodológica deste trabalho. São abordados, inicialmente, os princípios da montagem *de novo* de genomas e as especificidades do DNA mitocondrial, que justificam sua relevância como objeto de estudo. Em seguida, discute-se a questão da reprodutibilidade em ciência computacional e as soluções contemporâneas oferecidas pela engenharia de software, com destaque para a tecnologia de conteinerização. Por fim, apresentam-se os sistemas de gerenciamento de workflows científicos, cuja integração com containers constitui a base para a construção de pipelines modernos, portáveis e reprodutíveis.

### 2.1 Montagem *De Novo* de Genomas e a Análise Mitocondrial

Em bioinformática, a montagem de genomas consiste no processo de reconstrução da sequência original de DNA a partir de fragmentos obtidos por tecnologias de Sequenciamento de Nova Geração (NGS). O método mais amplamente utilizado para gerar esses fragmentos é o *shotgun sequencing*, no qual o DNA genômico é clivado aleatoriamente em milhões de pequenos pedaços que, posteriormente, são sequenciados. Essa estratégia garante alta cobertura do genoma, mas resulta em dados fragmentados que precisam ser cuidadosamente alinhados e sobrepostos para permitir a reconstrução da sequência completa (EKBLOM; WOLF, 2014). A complexidade computacional dessa tarefa motivou o desenvolvimento de algoritmos especializados, baseados em técnicas como grafos de De Bruijn e estratégias de sobreposição-layout-consenso (MILLER; KOREN; SUTTON, 2010).

**Figura 1** - Esquema do processo de shotgun sequencing.
Fonte: EKBLOM; WOLF (2014).

Na montagem genômica, fragmentos sobrepostos são inicialmente agrupados em contigs, sequências contínuas que representam trechos do genoma. Quando existe informação adicional que permite inferir a ordem relativa entre diferentes contigs - como dados de leituras pareadas (paired-end ou mate-pair) -, esses podem ser organizados em estruturas mais amplas chamadas scaffolds. Assim, enquanto os contigs representam sequências consolidadas diretamente a partir das leituras, os scaffolds constituem hipóteses mais abrangentes sobre a organização genômica, servindo como passo intermediário até a obtenção de montagens completas e anotadas.

As tecnologias de NGS diferem substancialmente no comprimento das leituras produzidas. Sequenciadores Illumina, por exemplo, geram leituras curtas de alta acurácia (50 a 300 pares de bases), que exigem algoritmos sofisticados para lidar com regiões repetitivas. Em contrapartida, tecnologias como PacBio e Oxford Nanopore produzem leituras longas, que podem alcançar milhares de pares de bases em uma única sequência, favorecendo a continuidade da montagem, ainda que a custos mais elevados e com menor rendimento global (EKBLOM; WOLF, 2014).

**Figura 2** - Comparação entre montagens SR e LR.
Fonte: ORELLANA et al. (2023).

A Figura 2 - Comparação de metagenome-assembled genomes (MAGs) recuperados a partir de tecnologias de short reads (Illumina) e long reads (PacBio HiFi). Genomas derivados de long reads apresentam maior continuidade (N50 mais elevado), menor fragmentação (menos contigs) e maior número de genes preditos em comparação aos derivados de short reads (ORELLANA et al., 2023) - evidencia como a escolha da tecnologia influencia diretamente a qualidade da montagem. Enquanto os short reads permitem recuperar um número maior de genomas devido à profundidade de sequenciamento, eles produzem montagens mais fragmentadas. Já os long reads oferecem menor diversidade detectada, mas geram montagens mais contínuas e completas, preservando genes essenciais, como o 16S rRNA, e representando melhor regiões repetitivas.

O DNA mitocondrial (mtDNA) reúne um conjunto de particularidades que o tornam um alvo privilegiado para a montagem *de novo*. Em metazoários, ele corresponde a uma molécula circular de fita dupla, com tamanho geralmente entre 12 e 22 mil pares de bases, codificando 37 genes distribuídos em 13 proteínas associadas à fosforilação oxidativa, 22 RNAs transportadores (tRNAs) e 2 RNAs ribossômicos (rRNAs) (ANDERSON et al., 1981; ASAKAWA et al., 1995; BOORE, 1999). O sequenciamento completo do genoma mitocondrial humano por Anderson et al. (1981) representou um marco histórico, estabelecendo o modelo de referência para a organização gênica mitocondrial dos metazoários. A Figura 3, a seguir, apresenta a representação circular típica de um genoma mitocondrial, evidenciando essa organização gênica compacta.

**Figura 3** - Representação circular do genoma mitocondrial.
Fonte: Do próprio autor (2024).

A ausência de íntrons, a estrutura compacta e a elevada taxa de cópias por célula tornam o mtDNA especialmente acessível, facilitando sua recuperação mesmo em experimentos voltados ao genoma nuclear, nos quais ele aparece como subproduto. Além disso, características funcionais como a herança predominantemente materna, a ausência geral de recombinação e a taxa evolutiva superior à do DNA nuclear ampliam sua utilidade em estudos biológicos. Essas propriedades explicam sua ampla aplicação em análises de genética de populações, filogeografia, sistemática e conservação da biodiversidade.

**Figura 4** - Exemplo de arquivo FASTA.
Fonte: Do próprio autor (2024).

A figura 4 representa um trecho de arquivo FASTA do genoma mitocondrial da Anodorhynchus leari, representando a forma textual dos dados brutos utilizados como entrada em pipelines de montagem genômica. Esse exemplo evidencia como os dados genômicos são inicialmente armazenados como cadeias de nucleotídeos em arquivos de texto simples, que posteriormente são processados e reconstruídos em estruturas circulares biologicamente coerentes. A contraposição entre o dado cru (FASTA) e a representação organizada (Figura 3) demonstra a importância dos algoritmos de montagem para transformar informação bruta em conhecimento biológico.

Para explorar essas vantagens, foram desenvolvidos algoritmos específicos de montagem. Uma das abordagens mais empregadas é a *seed-and-extend*, em que uma sequência inicial conhecida (*seed*) — que pode ser um gene ou até mesmo um fragmento de organela de espécie próxima — serve como ponto de partida para estender iterativamente a montagem até que a molécula circular seja reconstruída. O NOVOPlasty é um dos principais programas baseados nessa estratégia, destacando-se pela eficiência na recuperação de genomas mitocondriais e cloroplastidiais a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017). Avanços recentes, como o pipeline MitoHiFi, exploram leituras longas de alta fidelidade, evidenciando a constante evolução das ferramentas voltadas à montagem de organelas (ULIANO-SILVA et al., 2023).

Complementarmente à montagem, a anotação funcional dos genes mitocondriais permite identificar e classificar os elementos codificantes da molécula reconstruída. Para metazoários, ferramentas como o MITOS2 combinam alinhamentos de homologia com modelos de covariância (Infernal) para predizer genes codificadores de proteínas, tRNAs e rRNAs, produzindo anotações em formatos padronizados como GFF3 (DONATH et al., 2019). A integração da etapa de anotação em pipelines automatizados de montagem permanece, contudo, uma lacuna na maioria das soluções disponíveis na literatura.

### 2.2 Reprodutibilidade em Ciência Computacional

A reprodutibilidade é um pilar do método científico, consistindo na capacidade de um pesquisador independente replicar um experimento, com os mesmos materiais e métodos, e obter resultados consistentes. Em pesquisas computacionais, isso se traduz na possibilidade de executar o mesmo código e software sobre os mesmos dados e chegar ao mesmo resultado (PENG, 2011; SANDVE et al., 2013).

Contudo, a ciência moderna enfrenta o que tem sido amplamente denominado como uma "crise de reprodutibilidade", na qual um número significativo de estudos publicados se mostra difícil ou impossível de ser replicado (BAKER, 2016). No domínio da bioinformática, essa crise é particularmente acentuada devido à alta complexidade dos pipelines de análise. Um único pipeline pode envolver dezenas de ferramentas de software distintas, cada uma desenvolvida por diferentes grupos, em diferentes linguagens de programação e com um conjunto único e, por vezes, conflitante de dependências de bibliotecas e do próprio sistema operacional. A tarefa de instalar e configurar manualmente esse ecossistema de software para replicar uma análise é extremamente suscetível a erros, uma condição frequentemente descrita como dependency hell (inferno de dependências) (GRÜNING et al., 2018).

Pequenas variações na versão de uma ferramenta ou de uma biblioteca subjacente podem levar a resultados drasticamente diferentes, comprometendo a validade e a confiabilidade da pesquisa. Essa barreira técnica não apenas dificulta a verificação dos resultados, mas também impede a reutilização e a adaptação de métodos computacionais, contrariando os princípios FAIR (*Findable, Accessible, Interoperable, and Reusable*), que visam maximizar o valor dos dados e das análises científicas (SANDVE et al., 2013; WILKINSON et al., 2016; GRÜNING et al., 2018; BRITISH ECOLOGICAL SOCIETY, 2017).

### 2.3 Tecnologia de Conteinerização com Docker

Como solução para o problema do dependency hell e para garantir a reprodutibilidade computacional, a comunidade científica tem adotado práticas e tecnologias da engenharia de software, com destaque para a virtualização em nível de sistema operacional, também conhecida como containerização. A plataforma Docker se estabeleceu como a principal ferramenta desse paradigma, permitindo o empacotamento de uma aplicação e de todo o seu ambiente de execução - incluindo bibliotecas, códigos-fonte e dependências - em uma unidade padronizada e isolada chamada container (MERKEL, 2014; GRÜNING et al., 2018).

O funcionamento do Docker baseia-se em dois conceitos centrais: a imagem e o container. A imagem é um template estático e imutável, construído a partir de um arquivo de texto com instruções sequenciais chamado Dockerfile. Esse arquivo serve como uma "receita" que define o sistema operacional base, as dependências a serem instaladas e os comandos a serem executados, criando um ambiente de software completo e autocontido. O container, por sua vez, é uma instância executável e isolada de uma imagem. Ao contrário de máquinas virtuais tradicionais, os containers compartilham o kernel do sistema operacional hospedeiro, o que os torna leves e eficientes para iniciar (BRITISH ECOLOGICAL SOCIETY, 2017).

Essa arquitetura garante que um software containerizado se comporte de maneira idêntica em qualquer ambiente que possua o Docker instalado - seja em computadores pessoais, servidores de produção ou ambientes de nuvem - resolvendo de forma eficaz o problema da inconsistência de ambientes e sendo um passo fundamental para a criação de pipelines científicos portáveis e reprodutíveis (DI TOMMASO et al., 2017).

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

Um trabalho de destaque na área é o MitoHiFi, um pipeline desenvolvido para a montagem de genomas mitocondriais a partir de dados de leituras longas de alta fidelidade (PacBio HiFi). Proposto por Uliano-Silva et al. (2023), o objetivo principal era criar um método que aproveitasse a extensão e a precisão deste tipo de dado para gerar mitogenomas completos e corretamente circularizados com maior eficiência. A metodologia do MitoHiFi integra ferramentas como o hifiasm-meta para montagem e o GetOrganelle para a etapa de finalização, sendo orquestrada por meio de scripts em Shell e utilizando Conda para o gerenciamento das dependências de software.

**Figura 4** - Página oficial do repositório MitoHiFi no GitHub.
Fonte: Uliano-Silva et al. (2023)

A abordagem do MitoHiFi é extremamente atual por focar em uma tecnologia de sequenciamento de ponta. No entanto, sua implementação via scripts e Conda, embora eficaz para seu propósito, pode apresentar limitações de portabilidade e reprodutibilidade quando comparada a soluções baseadas em containers e sistemas de workflow. O diferencial do presente TCC reside precisamente neste ponto: enquanto o MitoHiFi inova na aplicação de um novo tipo de dado biológico, nossa proposta foca em uma arquitetura de software mais robusta, utilizando Docker e Nextflow, para garantir que o pipeline de montagem (baseado em leituras curtas) seja universalmente reprodutível e portável.

### 3.2 MITObim: Estratégia de Baiting and Iterative Mapping

O MITObim, desenvolvido por Hahn et al. (2013), foi um dos primeiros métodos amplamente utilizados para montagem de organelas a partir de dados de NGS. Sua estratégia de baiting and iterative mapping consiste em iniciar o processo com uma sequência guia (um fragmento do genoma alvo ou de uma espécie próxima) e, a partir dela, recuperar leituras homólogas, refinando o alinhamento de forma iterativa até a montagem completa do mitogenoma.

**Figura 5** - Página oficial do repositório MITObim no GitHub.
Fonte: Hahn et al. (2013).

Embora tenha sido pioneiro na aplicação desta abordagem, o MITObim apresenta limitações em termos de eficiência computacional e escalabilidade, especialmente quando aplicado a conjuntos de dados maiores ou de maior complexidade. Apesar de sua relevância histórica, vem sendo gradualmente substituído por ferramentas mais modernas, como o NOVOPlasty e o GetOrganelle.

### 3.3 GetOrganelle: Montagem Baseada em Grafos

O GetOrganelle, proposto por Jin et al. (2020), é uma das ferramentas mais recentes e populares para montagem de organelas. Ele utiliza grafos de montagem derivados de dados de NGS para isolar e reconstruir genomas de organelas, oferecendo resultados robustos tanto para mitocôndrias quanto para cloroplastos. A ferramenta tem se destacado por sua precisão na recuperação de sequências circulares e por oferecer resultados consistentes mesmo em organismos com maior complexidade genômica.

**Figura 6** - Página oficial do repositório GetOrganelle no GitHub.
Fonte: Jin et al. (2020).

No entanto, assim como outras ferramentas tradicionais, o GetOrganelle depende de um ambiente de software adequado para ser executado, geralmente configurado manualmente pelo usuário ou via gerenciadores como Conda, o que ainda pode limitar sua reprodutibilidade plena em diferentes contextos computacionais.

### 3.4 MToolBox: Foco em Genomas Mitocondriais Humanos

O MToolBox, apresentado por Calabrese et al. (2014), foi uma das primeiras pipelines integradas voltadas especificamente para o estudo de genomas mitocondriais humanos. Ele combina montagem, anotação e análise funcional em um único fluxo, permitindo identificar variantes mitocondriais relevantes para estudos biomédicos e populacionais.

Apesar de sua contribuição para o campo, o MToolBox é fortemente orientado ao contexto humano e não é tão amplamente aplicável a outros organismos. Além disso, por ter sido desenvolvido antes da adoção massiva de containers e workflows, carece de recursos de portabilidade e reprodutibilidade presentes em soluções mais modernas.

**Figura 7** - Página oficial do repositório MToolBox no GitHub.
Fonte: Calabrese et al. (2014).

### 3.5 NOVOPlasty: Estratégia Seed-and-Extend

O NOVOPlasty, descrito por Dierckxsens et al. (2017), é atualmente uma das ferramentas mais utilizadas para montagem de genomas mitocondriais e cloroplastidiais a partir de dados de NGS de leituras curtas. Baseado na estratégia seed-and-extend, o programa inicia a montagem a partir de uma sequência semente fornecida pelo usuário e estende iterativamente o genoma até que a molécula circular esteja completa. Sua popularidade decorre da facilidade de uso e da capacidade de gerar montagens confiáveis mesmo a partir de dados de genoma total.

No entanto, o NOVOPlasty, como software isolado, ainda demanda configuração manual do ambiente e execução direta pelo usuário, o que pode dificultar sua integração em fluxos mais complexos e comprometer a reprodutibilidade quando diferentes versões ou configurações são utilizadas. É justamente neste ponto que se insere o presente TCC: ao invés de propor um novo montador, o objetivo é integrar o NOVOPlasty em um ecossistema de containers e workflow (Docker + Nextflow), assegurando que seu uso seja portável, automatizado e reprodutível em diferentes ambientes de execução.

**Figura 8** - Página oficial do repositório NOVOPlasty no GitHub.
Fonte: Dierckxsens et al. (2017).

### 3.6 Análise Comparativa

A Tabela a seguir sintetiza as características das ferramentas e pipelines analisados, confrontando-as com a proposta deste trabalho. A comparação considera critérios relevantes para a reprodutibilidade, portabilidade e completude do fluxo analítico.

| Critério | MitoHiFi | MITObim | GetOrganelle | MToolBox | NOVOPlasty | **Este trabalho** |
|---|---|---|---|---|---|---|
| **Tipo de dado** | Long reads (HiFi) | Short reads | Short/Long reads | Short reads | Short reads | Short reads |
| **Montagem** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Anotação funcional** | — | — | — | ✓ (humano) | — | ✓ (metazoários) |
| **GenBank Flat File** | — | — | — | — | — | ✓ |
| **Containerização** | Parcial (Conda) | — | Conda | — | — | Docker (5 imagens) |
| **Workflow manager** | — | — | — | — | — | Nextflow DSL2 |
| **Pilot QC automático** | — | — | — | — | — | ✓ |
| **Estrutura secundária RNA** | — | — | — | — | — | ✓ (ViennaRNA) |
| **Código aberto** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

A análise comparativa evidencia que, embora cada ferramenta contribua com inovações relevantes para domínios específicos, nenhuma oferece um fluxo integrado que cubra todas as etapas — desde a aquisição dos dados brutos até a geração de arquivos prontos para submissão ao GenBank — em um ambiente completamente containerizado e orquestrado. A proposta deste trabalho preenche essa lacuna ao combinar montagem, anotação funcional, geração de entregáveis e práticas de reprodutibilidade em um único pipeline modular.

A análise dos trabalhos relacionados permitiu identificar diferentes abordagens para a montagem de organelas, incluindo pipelines baseados em leituras longas, estratégias iterativas e montadores especializados como o NOVOPlasty. Embora cada solução apresente contribuições relevantes, permanecem lacunas quanto à padronização, portabilidade, completude do fluxo analítico e reprodutibilidade. Nesse contexto, o Capítulo 4 descreve a metodologia adotada neste trabalho, detalhando a concepção do pipeline proposto e as ferramentas que o compõem.

---

## Capítulo 4 — Metodologia

Nesta seção, é apresentado o plano metodológico para o desenvolvimento do pipeline de montagem de mitogenomas. A abordagem combina etapas de análise bioinformática com práticas consolidadas de engenharia de software, seguindo um processo de desenvolvimento sistemático descrito na seção 4.1, antes de detalhar a arquitetura do pipeline e suas etapas constituintes.

### 4.1 Processo de Desenvolvimento de Software

O desenvolvimento do pipeline seguiu uma abordagem de **desenvolvimento iterativo e incremental**, metodologia consolidada na engenharia de software em que o sistema é construído em ciclos sucessivos, cada um produzindo uma versão funcional e testável do produto (LARMAN; BASILI, 2003). Essa abordagem se contrapõe ao modelo cascata (*waterfall*), no qual todas as fases — levantamento de requisitos, projeto, implementação e teste — são executadas sequencialmente e o produto só é validado ao final. No contexto de software científico, em que os requisitos frequentemente evoluem à medida que os resultados são analisados, o modelo iterativo é amplamente recomendado (WILSON et al., 2014).

O ciclo de desenvolvimento foi organizado em duas fases principais, cada uma constituindo uma iteração completa com entrega funcional:

- **Fase 1 — Montagem**: implementação dos módulos de download (SRA-Toolkit), controle de qualidade (FastQC, Trim Galore) e montagem (NOVOPlasty), com containerização em Docker e orquestração via Nextflow. Ao final desta fase, o pipeline já era capaz de produzir mitogenomas completos e circularizados a partir de um acesso SRA.
- **Fase 2 — Anotação e entregáveis**: adição dos módulos de análise piloto de qualidade (Pilot QC), anotação funcional (MITOS2), geração de estruturas secundárias de RNA (ViennaRNA/RNAplot), produção de arquivos para submissão ao GenBank (GenBank Flat File) e compilação automática de entregáveis.

Cada fase produziu um pipeline funcional e independente. A Fase 1 foi preservada como branch separada no repositório Git (`fase1-montagem`), permitindo uso e validação independente por terceiros — uma prática de versionamento semântico que garante rastreabilidade e facilita a manutenção (SOMMERVILLE, 2016).

#### 4.1.1 Validação por Espécie-Controle

Um elemento central da estratégia de verificação e validação (V&V) foi o emprego de uma **espécie-controle** com mitogenoma previamente depositado em banco de dados público. Antes de aplicar o pipeline à espécie-alvo (*Anodorhynchus leari*), cuja sequência mitocondrial era inédita, o sistema foi integralmente validado com *Diploprion bifasciatus*, cujo mitogenoma de referência (PZ143763.1) está depositado no GenBank. Essa estratégia permitiu comparar a montagem obtida pelo pipeline com a sequência de referência conhecida, verificando a correção funcional do sistema em condições controladas antes de aplicá-lo a dados sem referência prévia.

No desenvolvimento de software científico, a validação com dados de referência conhecida é considerada uma das práticas mais importantes para garantir a confiabilidade dos resultados computacionais (WILSON et al., 2014). Essa abordagem é análoga ao conceito de *oracle testing* na engenharia de software, em que uma fonte externa de verdade é utilizada para validar as saídas do sistema sob teste (SOMMERVILLE, 2016).

#### 4.1.2 Tratamento de Erros e Correções Iterativas

O modelo iterativo também se manifestou na resolução de problemas técnicos encontrados durante a execução do pipeline. Exemplos concretos incluem:

1. **Incompatibilidade do SRA-Toolkit**: a descoberta de que o `fasterq-dump` do SRA-Toolkit 3.x não implementa o flag `-X` para limitação de leituras na origem, corrigida pela adoção de truncamento pós-conversão utilizando o utilitário `head`;
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

De forma geral, o pipeline é composto por sete macroetapas principais. A primeira é a **análise piloto de qualidade (Pilot QC)**, um estágio opcional que analisa automaticamente uma subamostra dos dados para determinar o volume ideal de leituras a serem processadas, evitando downloads desnecessários em datasets grandes. A segunda etapa é a **aquisição e preparo dos dados**, responsável por obter datasets públicos em repositórios como o Sequence Read Archive (SRA) do NCBI (LEINONEN; SUGAWARA; SHUMWAY, 2011). Em seguida, realiza-se o **controle de qualidade e pré-processamento**, que envolve a avaliação das leituras por meio do FastQC (ANDREWS, 2010) e a remoção de adaptadores e bases de baixa qualidade utilizando o Trim Galore (BABRAHAM BIOINFORMATICS, 2019), que por sua vez utiliza o Cutadapt (MARTIN, 2011) internamente.

A quarta etapa corresponde à **montagem do mitogenoma**, realizada com o NOVOPlasty, ferramenta especializada que utiliza a estratégia seed-and-extend para reconstruir a molécula circular completa a partir de leituras curtas (DIERCKXSENS; MARDULYN; SMITS, 2017). A quinta etapa, de **anotação funcional**, identifica os genes do genoma montado, utilizando o MITOS2 (DONATH et al., 2019) para a predição de genes codificadores de proteínas, rRNAs e tRNAs, com geração de estruturas secundárias previstas via ViennaRNA (LORENZ et al., 2011). A sexta etapa consiste na **compilação dos entregáveis**, em que um script automatizado reúne todos os produtos da análise em uma pasta organizada, incluindo arquivos no formato GenBank Flat File para submissão ao NCBI. Por fim, a sétima etapa concentra-se na **automação e reprodutibilidade**, em que cada componente do pipeline é containerizado com Docker (MERKEL, 2014) e orquestrado via Nextflow (DI TOMMASO et al., 2017).

> **[SUGESTÃO DE FIGURA]**: Atualizar a Figura 9 (visão geral do pipeline) para incluir as 7 etapas: Pilot QC → SRA Download → FastQC → Trim Galore → NOVOPlasty → MITOS2 → Compile Summary. Usar o formato de fluxograma já existente, adicionando as etapas novas com cores distintas para as fases 1 e 2.

**Figura 9** - Visão geral do Pipeline implementado.
Fonte: Do próprio autor (2026).

A Figura 9 apresenta a visão geral do pipeline, representada em formato de fluxograma, ilustrando o encadeamento das etapas e destacando a integração entre análise bioinformática e práticas de engenharia de software. Deve-se observar que a lógica condicional do pipeline permite que a etapa de anotação (MITOS2) seja executada apenas quando o banco de dados de referência estiver disponível, tornando a Fase 1 (montagem) independente da Fase 2 (anotação). Adicionalmente, o workflow aplica um **filtro de circularização**: apenas montagens cujo nome de arquivo inicia com `Circularized_assembly` são encaminhadas ao MITOS2 e à compilação de entregáveis. Quando o NOVOPlasty não obtém circularização (produzindo apenas contigs parciais), a anotação é automaticamente omitida, evitando análises sobre montagens incompletas que poderiam gerar resultados biologicamente incorretos. Essa decisão de design garante que os entregáveis finais representem exclusivamente genomas mitocondriais completos e circulares.

### 4.3 Ferramentas e Tecnologias Utilizadas

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

#### 4.3.1 SRA-Toolkit

O SRA-Toolkit é o pacote oficial de utilitários do NCBI para manipulação de dados depositados no Sequence Read Archive (SRA), o principal repositório público de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). No pipeline desenvolvido, dois utilitários do SRA-Toolkit são empregados em sequência: o `prefetch`, que realiza o download eficiente do arquivo `.sra` compactado, e o `fasterq-dump`, que converte esse arquivo para o formato FASTQ paired-end, compatível com as etapas subsequentes do pipeline.

Uma decisão de implementação relevante diz respeito à estratégia de truncamento de dados. Em datasets de grande volume — como o utilizado para *Anodorhynchus leari*, com 118,5 milhões de pares de leituras (SRR28399504) —, o processamento integral seria desnecessário e computacionalmente oneroso, dado que o genoma mitocondrial de aves tem cerca de 17 kb e apenas uma fração diminuta das leituras (~0,17%) corresponde ao mtDNA. O pipeline implementa, portanto, uma truncagem pós-conversão utilizando o utilitário `head`, limitando os arquivos FASTQ ao número de leituras definido pelo parâmetro `sra_max_reads`. Essa abordagem foi adotada porque o `fasterq-dump` do SRA-Toolkit 3.x não implementa o flag `-X` para limitação na origem, diferentemente do `fastq-dump` legado, cuja execução é significativamente mais lenta.

A escolha de processar a conversão completa seguida de truncamento, em vez de utilizar o `fastq-dump` com `-X`, é justificada pela diferença de desempenho: o `fasterq-dump` é multi-threaded e até 10 vezes mais rápido que o `fastq-dump` em datasets grandes (NCBI, 2023). O custo adicional de disco temporário é compensado pela economia de tempo e pela confiabilidade da operação, que é automatizada integralmente pelo workflow.

#### 4.3.2 FastQC

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

#### 4.3.3 Cutadapt e Trim Galore

Após a avaliação inicial da qualidade das leituras, pode ser necessário realizar etapas de filtragem para remover adaptadores residuais e bases de baixa qualidade. Nesse contexto, duas ferramentas complementares são utilizadas no pipeline: Cutadapt e Trim Galore.

O Cutadapt é um programa desenvolvido por Martin (2011) especificamente para identificar e remover sequências de adaptadores em dados de alto rendimento. As leituras geradas por tecnologias como Illumina frequentemente contêm fragmentos de adaptadores ligados durante a preparação da biblioteca, que, se não forem removidos, podem comprometer a montagem ou o mapeamento subsequente. O Cutadapt utiliza algoritmos de busca eficiente para detectar esses fragmentos nas extremidades das leituras e removê-los, além de permitir o corte de bases de baixa qualidade de acordo com pontuações Phred. Como entrada, o programa recebe arquivos FASTQ, e como saída gera novos FASTQ já filtrados, prontos para análises subsequentes.

O Trim Galore é um wrapper que automatiza a execução do Cutadapt em conjunto com o FastQC, simplificando o processo de trimming e tornando-o mais acessível ao usuário (BABRAHAM BIOINFORMATICS, 2019). Com ele, é possível remover adaptadores mesmo quando a sequência exata não é previamente conhecida, o que é particularmente útil em experimentos com bibliotecas de origem diversa. Além disso, o Trim Galore aplica rotinas de filtragem por qualidade e integra relatórios do FastQC para avaliar o impacto do corte.

Entre as principais vantagens do Cutadapt e do Trim Galore estão: (i) a remoção precisa de adaptadores, prevenindo erros de montagem, (ii) a flexibilidade de configuração de parâmetros, e (iii) a integração automatizada que reduz a intervenção manual. Suas limitações incluem a necessidade de escolher cuidadosamente os limiares de qualidade, já que cortes excessivos podem resultar em perda de informação biológica relevante, enquanto cortes insuficientes podem deixar contaminantes residuais.

A escolha por essas ferramentas neste trabalho se justifica pelo seu amplo uso na comunidade científica e pela confiabilidade de seus resultados. A integração de Cutadapt e Trim Galore garante que as leituras que avançam para a montagem apresentem qualidade adequada, reduzindo o risco de artefatos e aumentando a robustez dos resultados obtidos pelo pipeline.

#### 4.3.4 NOVOPlasty

O NOVOPlasty é um montador especializado desenvolvido para reconstrução de genomas de organelas, como mitocôndrias e cloroplastos, a partir de dados de sequenciamento de genoma total. Publicado por Dierckxsens, Mardulyn e Smits (2017), tornou-se uma das ferramentas mais amplamente utilizadas para esse fim devido à sua eficiência e facilidade de uso.

Do ponto de vista técnico, o NOVOPlasty baseia-se na estratégia seed-and-extend. O processo de montagem é iniciado a partir de uma sequência fornecida pelo usuário — denominada seed — que pode ser um fragmento do genoma da organela em estudo, uma sequência de uma espécie próxima ou até mesmo uma sequência completa de organela de outro organismo. A partir dessa seed, o software realiza uma extensão bidirecional, identificando leituras sobrepostas e construindo gradualmente a sequência circular completa do genoma. Para acelerar o processo, o programa utiliza tabelas de hash, que permitem buscas rápidas por sobreposições entre as leituras.

O NOVOPlasty recebe como entrada arquivos FASTQ contendo leituras pareadas, bem como um arquivo FASTA com a sequência seed. A saída consiste em um arquivo FASTA representando o genoma circularizado, acompanhado de arquivos de log e, em alguns casos, múltiplas opções de montagem (quando diferentes caminhos são possíveis no grafo de montagem).

Entre as principais vantagens do NOVOPlasty destacam-se: (i) sua alta eficiência na recuperação de genomas completos de organelas mesmo a partir de dados de genoma total, (ii) a baixa exigência de parâmetros complexos, o que o torna acessível a usuários com diferentes níveis de experiência, e (iii) a capacidade de lidar com diferentes tipos de sementes, permitindo flexibilidade em estudos de espécies com genomas pouco conhecidos. No entanto, como limitações, pode-se destacar: (i) a dependência da qualidade e da escolha adequada da seed, que influencia diretamente o sucesso da montagem, e (ii) o risco de circularizações incorretas em casos de regiões altamente repetitivas ou coberturas desbalanceadas.

No pipeline implementado, o módulo NOVOPlasty incorpora uma lógica iterativa com re-seeding automático: caso a primeira execução não resulte em circularização, o maior contig obtido é automaticamente utilizado como nova semente em uma segunda tentativa, e múltiplos valores de k-mer são experimentados em sequência. Essa estratégia aumenta a taxa de sucesso da montagem sem requerer intervenção manual do usuário.

#### 4.3.5 MITOS2

O MITOS2 é uma ferramenta desenvolvida para a anotação funcional de genomas mitocondriais em metazoários. Apresentada por Donath et al. (2019), constitui a evolução da primeira versão do MITOS, oferecendo maior precisão na predição de genes, integração com bases de dados atualizadas e maior flexibilidade de execução.

Tecnicamente, o MITOS2 combina métodos de homologia e algoritmos de predição ab initio para identificar as principais classes de genes presentes em genomas mitocondriais. O sistema utiliza alinhamentos com modelos de substituição específicos para proteínas codificadas por mitocôndrias, além de perfis de RNA para identificar rRNAs e tRNAs por meio de modelos de covariância (Infernal). As entradas consistem em um arquivo FASTA contendo a sequência do genoma mitocondrial a ser anotado, e as saídas incluem: (i) arquivos tabulares com as coordenadas dos genes anotados, (ii) representações gráficas do genoma mostrando a organização gênica, e (iii) arquivos de anotação em formatos padrão, como GFF3 e BED.

A versão empregada neste trabalho (v2.1.9) foi instalada via Conda em um container Docker baseado em Miniconda 24.1.2-0, utilizando o banco de dados RefSeq89m obtido do Zenodo (DONATH et al., 2019). Essa abordagem garante que a ferramenta, normalmente executada via interface web, funcione de forma autônoma em linha de comando e completamente integrada ao workflow Nextflow.

Como extensão funcional, o pipeline integra automaticamente a geração de representações visuais da estrutura secundária dos tRNAs e rRNAs preditos. Para isso, o módulo MITOS2 extrai as sequências e estruturas dot-bracket dos arquivos de saída e as processa com o RNAplot, componente do pacote ViennaRNA (LORENZ et al., 2011), gerando arquivos SVG individuais para cada RNA não codificante. Essa funcionalidade, originalmente disponível apenas na versão web do MITOS2 (UseGalaxy), foi reproduzida localmente no pipeline, garantindo que os diagramas de estrutura secundária sejam gerados sem dependência de serviços externos.

Entre as principais vantagens do MITOS2 destacam-se: (i) a automação completa da etapa de anotação, reduzindo a necessidade de curadoria manual extensiva; (ii) a compatibilidade com diferentes grupos taxonômicos de metazoários; e (iii) a produção de relatórios gráficos de fácil interpretação. Por outro lado, como limitações, pode apresentar divergências em regiões com estruturas gênicas incomuns — como o caso do gene ND3 em aves, que apresenta um frameshift insertion documentado (MINDELL et al., 1998), levando o MITOS2 a reportar o gene como dois ORFs distintos (nad3_0 e nad3_1) em vez de um único CDS com `join()`. Esse tipo de situação requer validação complementar e tratamento específico na conversão para formatos de submissão ao GenBank.

#### 4.3.6 Docker

No pipeline proposto, cada ferramenta foi encapsulada em um container Docker, assegurando que todas as dependências fossem executadas em ambientes controlados e consistentes (MERKEL, 2014). O pipeline utiliza cinco imagens Docker distintas, cada uma correspondendo a um estágio do processamento:

| Imagem | Base | Ferramenta | Dockerfile |
|---|---|---|---|
| `mitogenome-pipeline/sra-tools:1.0` | ubuntu:22.04 | SRA-Toolkit 3.0.10 | `docker/sra-tools/Dockerfile` |
| `mitogenome-pipeline/fastqc:1.0` | ubuntu:22.04 | FastQC 0.12.1 | `docker/fastqc/Dockerfile` |
| `mitogenome-pipeline/trim-galore:1.0` | ubuntu:22.04 | Trim Galore 0.6.10 + Cutadapt 4.6 | `docker/trim-galore/Dockerfile` |
| `mitogenome-pipeline/novoplasty:1.0` | ubuntu:22.04 | NOVOPlasty 4.3.1 | `docker/novoplasty/Dockerfile` |
| `mitogenome-pipeline/mitos2:1.0` | miniconda3:24.1.2-0 | MITOS2 2.1.9 + ViennaRNA | `docker/mitos2/Dockerfile` |

Essa abordagem eliminou a necessidade de configurações manuais complexas e garantiu que softwares como FastQC, Cutadapt, NOVOPlasty e MITOS2 operassem de forma uniforme em diferentes sistemas operacionais. Com isso, a montagem e a anotação de mitogenomas tornam-se replicáveis em qualquer ambiente computacional que disponha do Docker instalado, ampliando a portabilidade do pipeline.

#### 4.3.7 Nextflow

A orquestração do pipeline foi implementada no Nextflow DSL2, que estruturou os processos em módulos independentes (um arquivo `.nf` por etapa), conectando automaticamente suas entradas e saídas (DI TOMMASO et al., 2017). Essa organização permitiu automatizar a lógica iterativa da montagem, como a possibilidade de reaproveitar contigs parciais do NOVOPlasty como novas seeds em execuções subsequentes. Além disso, o Nextflow possibilitou a paralelização de tarefas, a retomada de execuções interrompidas (via flag `-resume`) e a integração direta com containers Docker. Dessa forma, o pipeline desenvolvido combina simplicidade de uso com escalabilidade e reprodutibilidade, adequando-se tanto a computadores locais quanto a servidores de alto desempenho.

O sistema de perfis do Nextflow foi utilizado para parametrizar execuções distintas sem alterar o código-fonte. Cada espécie estudada possui um arquivo de configuração próprio (e.g., `conf/a_leari.config`) contendo os parâmetros específicos — acesso SRA, semente, tamanho esperado do genoma, código genético — que são carregados automaticamente ao invocar o perfil correspondente (`nextflow run main.nf -profile a_leari`).

#### 4.3.8 GitHub

O GitHub é uma plataforma de hospedagem de código baseada no sistema de controle de versão Git, amplamente utilizada para o desenvolvimento colaborativo de software e projetos científicos. Para pipelines de bioinformática, a disponibilização em um repositório GitHub garante não apenas transparência, mas também acessibilidade e manutenção a longo prazo.

Neste trabalho, o GitHub foi utilizado como repositório público para disponibilizar o pipeline desenvolvido, assegurando sua acessibilidade, versionamento e documentação adequada. Essa prática se alinha às diretrizes de ciência aberta e reforça o compromisso com a reprodutibilidade e a transparência metodológica. O repositório está disponível em https://github.com/matheus-sobreira/mitogenome-pipeline e inclui: o código completo do workflow Nextflow, todos os Dockerfiles, scripts auxiliares, configurações de perfil por espécie, dados de semente (sequências cox1), documentação e guia de execução.

Adicionalmente, o sistema de branches do Git foi utilizado para preservar estados estáveis do pipeline: a branch `fase1-montagem` contém o estado validado da Fase 1 (apenas montagem, sem anotação), enquanto a branch `main` contém o pipeline completo com todas as etapas integradas.

### 4.4 Análise Piloto de Qualidade (Pilot QC)

Uma inovação incorporada ao pipeline é o módulo de análise piloto de qualidade, que soluciona um problema prático recorrente em montagens de organelas: a determinação do volume adequado de dados a serem processados. Em datasets de sequenciamento de genoma total, a proporção de leituras mitocondriais é tipicamente inferior a 1%, o que significa que processar integralmente um dataset de 100 GB resultaria em centenas de gigabytes de dados nucleares irrelevantes para a montagem do mitogenoma.

O módulo Pilot QC opera da seguinte forma: quando o parâmetro `sra_max_reads` não é definido pelo usuário, o pipeline baixa automaticamente uma subamostra de 500.000 leituras (aproximadamente 500 MB) utilizando o `fastq-dump -X`, que limita o download na origem sem necessidade de transferir o dataset completo. A escolha do `fastq-dump` (em vez do `fasterq-dump` utilizado no módulo de download principal) é deliberada: o `fastq-dump`, embora mais lento, suporta o flag `-X` que restringe o download na origem, evitando transferir o dataset completo. Já o `fasterq-dump`, utilizado no download principal por ser multi-threaded e até 10 vezes mais rápido (NCBI, 2023), não implementa esse flag, exigindo a estratégia de truncamento pós-conversão via `head`. Para amostras piloto de 500 K reads, a menor velocidade do `fastq-dump` é irrelevante (execução em segundos), tornando a escolha ótima para cada contexto. Essa subamostra é então submetida ao script `scripts/pilot_qc.sh`, um analisador desenvolvido inteiramente em Bash que opera diretamente sobre os arquivos FASTQ sem dependências externas além de `awk` e `grep`. O script extrai três métricas-chave:

1. **Proporção de bases com qualidade ≥ Q30**, indicando a qualidade geral das leituras;
2. **Proporção de leituras com adaptador Illumina TruSeq** residual, estimando a perda efetiva no trimming;
3. **Fração mitocondrial estimada**, calculada por correspondência de k-mers da sequência semente contra as leituras piloto.

A partir dessas métricas, o script calcula o número total de leituras necessário para atingir uma cobertura-alvo de 500× no mitogenoma, aplicando fatores de correção por qualidade e por presença de adaptador, com uma margem de segurança de 1,5×. O resultado é limitado ao intervalo de 5 a 50 milhões de leituras. Esse valor é então passado automaticamente ao módulo de download principal, que trunca os FASTQs no número calculado.

A fórmula de cálculo é:

$$R_{mito} = \frac{C_{alvo} \times G}{L \times 2}$$

$$R_{total} = \frac{R_{mito} \times F_{Q30} \times F_{adapter}}{f_{mito}} \times 1.5$$

Onde $C_{alvo}$ é a cobertura desejada, $G$ é o tamanho médio esperado do genoma, $L$ é o comprimento médio das leituras, $f_{mito}$ é a fração mitocondrial estimada, e $F_{Q30}$ e $F_{adapter}$ são fatores de correção baseados na qualidade e presença de adaptadores.

Caso o usuário opte por definir manualmente o `sra_max_reads` via parâmetro de linha de comando, a etapa Pilot QC é automaticamente pulada, preservando a flexibilidade para cenários em que o pesquisador já possui conhecimento prévio sobre seus dados.

### 4.4.1 Arquitetura Modular e Scripts do Pipeline

A organização do código segue a convenção de projetos Nextflow DSL2, em que cada etapa do pipeline corresponde a um **módulo** independente (arquivo `.nf` no diretório `modules/`), e os scripts auxiliares residem no diretório `scripts/`. Essa separação permite que cada componente seja desenvolvido, testado e mantido de forma isolada. A Tabela 2 apresenta a relação completa dos módulos, scripts e suas responsabilidades.

> **[SUGESTÃO DE TABELA — Tabela 2]**: Inserir tabela com a arquitetura modular:

| Arquivo | Tipo | Função |
|---|---|---|
| `main.nf` | Workflow | Orquestração geral: define a ordem de execução, lógica condicional (Pilot QC, MITOS2) e passagem de dados entre módulos |
| `modules/sra_download.nf` | Módulo | Download do SRA via `prefetch` + `fasterq-dump`, com truncamento por `head` |
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

### 4.5 Aquisição e Preparo dos Dados

A primeira etapa prática do pipeline consiste na aquisição e organização dos dados de sequenciamento que servirão de insumo para a montagem dos genomas mitocondriais. Para este trabalho, foram utilizados datasets públicos provenientes do Sequence Read Archive (SRA), repositório mantido pelo NCBI e reconhecido como um dos principais bancos internacionais de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). A escolha por dados públicos justifica-se pela sua ampla disponibilidade, pela diversidade de organismos representados e pela possibilidade de reproduzir e validar resultados obtidos em trabalhos prévios.

O acesso aos dados é realizado por meio do SRA-Toolkit v3.0.10 (encapsulado em container Docker), por meio de dois comandos sequenciais: `prefetch`, responsável pelo download do arquivo `.sra` compactado, e `fasterq-dump`, que converte o arquivo para o formato FASTQ paired-end. Implementou-se também a estratégia de retries automáticos (até 2 tentativas) para lidar com instabilidades de rede, comuns em downloads de datasets grandes do NCBI.

Como estudo de caso principal, foi utilizado o dataset SRR28399504 referente ao genoma de *Anodorhynchus leari* (arara-azul-de-lear), gerado por sequenciamento paired-end na plataforma Illumina HiSeq X Ten (IRIDIAN GENOMES, 2024), com leituras de 150 pares de bases e um total de 118,5 milhões de pares. Para validação do pipeline, foi utilizado o dataset SRR36182901 referente a *Diploprion bifasciatus* (peixe-sabão barrado), gerado em plataforma NovaSeq X Plus com leituras de 151 bp.

### 4.6 Controle de Qualidade e Pré-processamento

Após a aquisição e organização dos dados de sequenciamento, a etapa seguinte do pipeline consiste na avaliação da qualidade das leituras e na aplicação de procedimentos de filtragem. Essa fase é essencial para assegurar que apenas dados confiáveis avancem para a montagem do genoma mitocondrial, evitando que artefatos técnicos comprometam a acurácia do resultado final.

O primeiro passo é a execução do FastQC v0.12.1, que gera relatórios interativos em HTML contendo métricas fundamentais para a avaliação da qualidade das leituras. Em seguida, é aplicado o Trim Galore v0.6.10 (que integra o Cutadapt v4.6), com parâmetros de qualidade Phred mínimo de 20 e comprimento mínimo de leitura de 50 bp para remoção de adaptadores e bases de baixa qualidade. No módulo implementado, o número de threads de CPU foi deliberadamente limitado a 2 para evitar ultrapassar a alocação de CPU no container, dado que o Trim Galore gera internamente threads adicionais para compressão (pigz) e para o próprio Cutadapt.

Além de aumentar a confiabilidade biológica das leituras, o pré-processamento contribui significativamente para a redução do volume de dados. Na execução para *A. leari*, por exemplo, o Trim Galore detectou adaptadores Illumina TruSeq em 81,4% das leituras, e após o tricming, 71,9% das bases foram mantidas — uma redução de volume de ~28%.

### 4.7 Montagem do Mitogenoma

A montagem do genoma mitocondrial constitui a etapa central do pipeline, pois é nela que se obtém a sequência circular completa que serve de base para a anotação funcional e análises comparativas. Para essa finalidade é utilizado o NOVOPlasty v4.3.1, um dos montadores mais consolidados para organelas (DIERCKXSENS; MARDULYN; SMITS, 2017).

No pipeline implementado, o módulo NOVOPlasty incorpora diversas melhorias em relação ao uso padrão da ferramenta:

1. **Múltiplos k-mers**: o parâmetro `novoplasty_kmers` permite definir uma lista de valores de k-mer a serem tentados em sequência (e.g., `'39,33'`). O pipeline experimenta cada um até obter circularização;
2. **Re-seeding automático**: para cada k-mer, caso a primeira execução não circularize, o maior contig obtido é automaticamente extraído e utilizado como nova semente em até 5 iterações (configurável);
3. **Liberação de espaço**: após a montagem, os arquivos FASTQ trimados são automaticamente deletados para economizar espaço em disco, comportamento crucial para execuções com dados volumosos.

A seleção dos valores de k-mer segue as recomendações do próprio NOVOPlasty (DIERCKXSENS; MARDULYN; SMITS, 2017), que aceita valores ímpares no intervalo de 21 a 39. Valores maiores de k-mer (como 39) favorecem a **especificidade** da montagem — cada k-mer é mais longo e, portanto, menos ambíguo, reduzindo o risco de quimeras e extensões incorretas. Por outro lado, valores menores (como 33) aumentam a **sensibilidade**, permitindo extensões em regiões de menor cobertura ou maior divergência entre a semente e o genoma-alvo. Assim, a estratégia adotada neste pipeline inicia pelo valor máximo (39), que produz montagens mais confiáveis quando a cobertura e a qualidade dos dados são adequadas, e recorre ao valor menor (33) apenas se a circularização não for alcançada na primeira tentativa. Para *A. leari*, a circularização ocorreu na primeira tentativa com k-mer 39, indicando que a cobertura de 306× e a qualidade do HiSeq X Ten foram suficientes para o valor mais restritivo. Já para *D. bifasciatus*, a circularização ocorreu com k-mer 33, possivelmente refletindo diferenças na cobertura ou na composição nucleotídica do dataset.

**Obtenção da semente (seed).** O NOVOPlasty requer uma sequência inicial conhecida para ancorar a extensão. A estratégia adotada neste trabalho consiste em: (i) identificar no NCBI uma espécie congênere ou filogeneticamente próxima que possua mitogenoma completo depositado; (ii) extrair o gene cox1 dessa referência, por ser o marcador mitocondrial mais conservado e amplamente disponível; e (iii) salvar a região em formato FASTA. Para *A. leari*, foi utilizado o gene cox1 de *A. hyacinthinus* (NC_082165.1, posições 5359–6906, 1.548 bp), espécie congênere cuja proximidade filogenética garante homologia suficiente para o seed-and-extend funcionar eficazmente. Para *D. bifasciatus*, foi utilizado o próprio cox1 da referência publicada (PZ143763.1, 1.560 bp). O repositório do pipeline inclui um guia detalhado para obtenção de sementes (`data/seeds/COMO_OBTER_SEMENTE.md`) e exemplos para sete espécies de diferentes grupos taxonômicos.

**Seleção do intervalo de tamanho (`genome_range`).** O NOVOPlasty utiliza um intervalo esperado de tamanho do genoma para descartar extensões espúrias. A determinação desse intervalo segue a mesma lógica da semente: busca-se no NCBI o tamanho do mitogenoma de referência da espécie congênere mais próxima e aplica-se uma margem de ±1.500–2.000 bp para acomodar variação interespecífica. Para *A. leari*, a referência *A. hyacinthinus* (16.999 bp) resultou no intervalo 15.500–18.500 bp; para *D. bifasciatus*, a referência PZ143763.1 (16.805 bp) resultou em 15.000–18.500 bp.

**Parâmetro `insert_size`.** O tamanho do inserto da biblioteca é parametrizado como 300 bp (valor padrão conservador para preparações Illumina). O NOVOPlasty utiliza esse valor como estimativa inicial, mas refina-o automaticamente durante a montagem — no caso de *A. leari*, o inserto real medido foi de 213 bp, sem impacto na qualidade da circularização.

### 4.8 Anotação Funcional

Concluída a etapa de montagem, procede-se à anotação funcional do genoma mitocondrial reconstruído, utilizando o MITOS2 v2.1.9 com o banco de dados RefSeq89m. A anotação é executada automaticamente pelo pipeline quando o parâmetro `mitos2_db` está configurado, sendo condicional para permitir que a Fase 1 (montagem) opere independentemente.

O módulo MITOS2 recebe como entrada o arquivo FASTA da montagem circularizada e produz como saída: arquivo GFF3 com coordenadas dos genes, tabelas BED, sequências proteicas preditas (FAA), mapa linear do genoma (PNG) e plots de qualidade (PDF). Adicionalmente, o pipeline gera automaticamente diagramas SVG da estrutura secundária prevista para cada tRNA e rRNA identificado, utilizando o programa RNAplot do pacote ViennaRNA.

O código genético é parametrizado (`genetic_code = 2` para vertebrados mitocondriais, `genetic_code = 5` para invertebrados), permitindo que o mesmo pipeline seja utilizado para diferentes grupos taxonômicos sem modificação de código.

### 4.9 Compilação dos Entregáveis

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

### 4.10 Automação e Reprodutibilidade

A automação do pipeline e a padronização de sua execução representam os principais diferenciais desta proposta em relação a abordagens convencionais de montagem de genomas mitocondriais. Mais do que reunir ferramentas em sequência, busca-se disponibilizar um fluxo de trabalho que possa ser aplicado por diferentes grupos de pesquisa em distintos ambientes computacionais, assegurando consistência analítica, transparência metodológica e facilidade de reutilização.

Para atingir esse objetivo, cada etapa do pipeline foi encapsulada em containers Docker, garantindo que as ferramentas operem em ambientes controlados e reproduzíveis, independentemente do sistema em que forem executadas (MERKEL, 2014). A orquestração pelo Nextflow automatiza a execução dos processos, gerencia dependências e possibilita paralelização quando aplicável. O sistema ainda permite a retomada da análise a partir do ponto de falha (via `-resume`), evitando reprocessamento desnecessário e aumentando a eficiência computacional (DI TOMMASO et al., 2017).

Todo o código-fonte do workflow está disponibilizado em repositório público no GitHub (https://github.com/matheus-sobreira/mitogenome-pipeline), acompanhado de documentação detalhada, guia de execução e exemplos de uso. Essa prática de ciência aberta não apenas reforça a transparência metodológica, mas também incentiva a adaptação e o aprimoramento do pipeline por outros pesquisadores, em consonância com iniciativas comunitárias como o nf-core (EWELS et al., 2020).

---

## Capítulo 5 — Resultados e Discussão

Este capítulo apresenta os resultados obtidos com a aplicação do pipeline desenvolvido, organizados segundo as etapas do fluxo de trabalho. Diferentemente da versão anterior deste trabalho, que apresentava resultados esperados de natureza projetiva, os dados aqui descritos são produtos reais de execuções completas do pipeline, validando sua funcionalidade e robustez.

### 5.1 Execuções Realizadas

O pipeline foi executado para duas espécies, com objetivos complementares:

| | *Anodorhynchus leari* | *Diploprion bifasciatus* |
|---|---|---|
| **Objetivo** | Objeto principal do TCC | Validação do pipeline |
| **Dataset SRA** | SRR28399504 | SRR36182901 |
| **Plataforma** | Illumina HiSeq X Ten | NovaSeq X Plus |
| **Leituras** | 118,5M × 150 bp | 12,0M × 151 bp |
| **Semente** | cox1 *A. hyacinthinus* (NC_082165.1, 1.548 bp) | cox1 *D. bifasciatus* (PZ143763.1, 1.560 bp) |
| **`genome_range`** | 15.500–18.500 bp | 15.000–18.500 bp |
| **K-mers** | 39, 33 | 33 |
| **Código genético** | 2 (vertebrado) | 2 (vertebrado) |
| **MITOS2** | Sim (RefSeq89m) | Não (perfil de teste sem banco) |
| **Perfil** | `-profile a_leari` | `-profile test` |

### 5.2 Aquisição e Preparo dos Dados

Para *A. leari*, o download via `prefetch` obteve o arquivo `.sra` de 11 GB, que foi convertido pelo `fasterq-dump` em dois arquivos FASTQ de aproximadamente 43,5 GB cada (R1 e R2). Após a truncagem para 20 milhões de leituras, os arquivos foram reduzidos para ~7,4 GB cada, representando uma economia de 83% em volume de dados processados sem qualquer prejuízo à montagem, dada a redundância intrínseca dos dados de genoma total.

> **[SUGESTÃO DE FIGURA]**: Inserir aqui o gráfico `pipeline_data_reduction.png`, já gerado, mostrando a redução de volume em cada etapa do pipeline. Legendar como "Figura X — Redução do volume de dados ao longo das etapas do pipeline para *A. leari*. Painel esquerdo: tamanho dos arquivos por etapa. Painel direito: composição de bases após trimming."

Para *D. bifasciatus*, o dataset compacto (12M reads, ~1,2 GB) foi processado integralmente sem necessidade de truncagem, utilizando como semente o gene cox1 da própria referência publicada (PZ143763.1, 1.560 bp) e `genome_range` de 15.000–18.500 bp.

### 5.3 Controle de Qualidade e Pré-processamento

Os relatórios do FastQC para *A. leari* indicaram qualidade per-base consistente (medianas superiores a Phred 30 ao longo de toda a extensão das leituras), padrão esperado para dados HiSeq X Ten. O Trim Galore identificou adaptadores Illumina TruSeq em 81,4% das leituras e, após o trimming, manteve 71,9% das bases originais. Esse percentual elevado de adaptador é atribuído ao tamanho curto dos insertos da biblioteca, característico de preparações para HiSeq X Ten, em que os fragmentos frequentemente ultrapassam o comprimento das leituras.

> **[SUGESTÃO DE FIGURA]**: Capturas de tela dos relatórios FastQC (antes e depois do trimming), mostrando a qualidade per-base e o conteúdo de adaptadores. Sugestão: montar como figura composta (2×2) com R1 antes, R1 depois, R2 antes, R2 depois.

### 5.4 Montagem do Mitogenoma

A montagem pelo NOVOPlasty resultou em genomas mitocondriais circularizados para ambas as espécies:

| Métrica | *A. leari* | *D. bifasciatus* |
|---|---|---|
| **Tamanho** | 16.986 bp | ~16.800 bp (circularizado) |
| **Circularização** | Sim | Sim |
| **K-mer ótimo** | 39 | 33 |
| **Cobertura média** | 306× | não medida (dataset completo) |
| **Reads alinhadas** | 34.644 / 20.285.220 (0,17%) | — |
| **Insert size medido** | 213 bp (parametrizado: 300) | — |
| **Referência** | 16.999 bp (*A. hyacinthinus*) — Δ 13 bp | 16.805 bp (PZ143763.1) — Δ ~5 bp |

Para *D. bifasciatus*, o dataset compacto (12M reads, ~1,2 GB) foi processado integralmente sem truncagem; a circularização ocorreu na primeira tentativa com k-mer 33. A etapa de anotação funcional (MITOS2) não foi executada para esta espécie no perfil de teste, que foi configurado sem o parâmetro `mitos2_db` para demonstrar a modularidade do pipeline — a Fase 1 (montagem) opera independentemente da Fase 2 (anotação). A validação da montagem pode ser realizada por BLAST contra a referência PZ143763.1.

O resultado de *A. leari* merece destaque: o genoma montado de 16.986 bp difere em apenas 13 pares de bases do mitogenoma de referência de *A. hyacinthinus* (16.999 bp, NC_082165.1), demonstrando alta congruência interspecífica esperada para o gênero *Anodorhynchus*. A cobertura média de 306× é amplamente superior ao mínimo recomendado de 30-50× para montagens confiáveis, conferindo alta confiabilidade ao resultado.

A proporção de leituras mitocondriais de 0,17% é condizente com a estimativa teórica para dados de genoma total de aves, onde o número de cópias mitocondriais por célula é tipicamente na ordem de centenas a milhares, mas o genoma nuclear (~1,2 Gb em psitacídeos) é muito maior que o mitocondrial (~17 kb).

O log do NOVOPlasty reportou uma fração de subamostragem (*subsampled fraction*) de 61,9%, indicando que o algoritmo processou 61,9% das leituras de entrada antes de obter a circularização. Esse valor reflete a eficiência do seed-and-extend: como o NOVOPlasty constrói o genoma iterativamente a partir da semente, ele pode alcançar a montagem completa sem necessidade de processar todas as leituras disponíveis, resultando em economia computacional proporcional à fração não processada.

> **[SUGESTÃO DE FIGURA]**: Inserir aqui uma representação do genoma circularizado montado. Pode ser o mapa circular gerado pelo Biopython (09_circular_map.svg/pdf) ou um mapa gerado no OGDRAW a partir do arquivo .gbk produzido pelo pipeline.

### 5.5 Anotação Funcional

A anotação pelo MITOS2 identificou o conjunto completo de genes esperado para um genoma mitocondrial de ave:

| Tipo | Quantidade | Genes |
|---|---|---|
| **CDS** | 13 | nad1, nad2, nad3, nad4, nad4l, nad5, nad6, cox1, cox2, cox3, atp6, atp8, cob |
| **tRNA** | 22 | tRNA-Phe, tRNA-Val, tRNA-Leu1, tRNA-Ile, tRNA-Gln, tRNA-Met, tRNA-Trp, tRNA-Ala, tRNA-Asn, tRNA-Cys, tRNA-Tyr, tRNA-Ser1, tRNA-Asp, tRNA-Lys, tRNA-Gly, tRNA-Arg, tRNA-His, tRNA-Ser2, tRNA-Leu2, tRNA-Glu, tRNA-Thr, tRNA-Pro |
| **rRNA** | 2 | rrnS (12S), rrnL (16S) |
| **Origem de replicação** | 2 | OH_0, OH_1 |

A ordem gênica identificada é compatível com o padrão conservado de Psittaciformes, corroborando a confiabilidade da montagem. A comparação com a anotação de *A. hyacinthinus* (NC_082165.1) confirmou a conservação da organização gênica dentro do gênero.

**Caso especial — Frameshift do ND3:** Um aspecto biologicamente relevante revelado pela anotação é o tratamento do gene ND3. O MITOS2 reportou este gene como dois ORFs distintos (nad3_0 e nad3_1) com sobreposição de 5 bp, refletindo o frameshift insertion documentado por Mindell et al. (1998) em genomas mitocondriais de aves. Esse fenômeno, em que uma inserção de ~1 nucleotídeo altera o quadro de leitura no interior do gene, é corrigido in vivo por mecanismos de programmed ribosomal frameshift ou edição do mRNA. Na referência *A. hyacinthinus* (NC_082165.1), o ND3 é anotado com `join()`:

```
CDS   join(9504..9677,9679..9854)
      /gene="ND3"
      /note="frameshift mechanism unknown (Mindell et al., 1998)"
```

Para a conversão dos resultados ao formato GenBank Flat File, o script `gff2genbank.py` do pipeline detecta automaticamente a presença de nad3_0 e nad3_1 e os unifica em um único CDS com `join()`, aplicando o qualificador `/exception=ribosomal slippage` — exatamente como exigido pelo GenBank para submissão. Quando a espécie não apresenta esse frameshift (e.g., invertebrados), o ND3 é anotado normalmente como um CDS contínuo.

> **[SUGESTÃO DE FIGURA]**: Inserir a tabela de posição física dos genes (06_gene_positions.tsv), similar à tabela do PDF de referência da disciplina, com colunas Region, Start, Stop, Strand, Length.

> **[SUGESTÃO DE FIGURA]**: Inserir uma ou mais imagens dos SVGs de estrutura secundária dos tRNAs gerados pelo pipeline (e.g., tRNA-Phe, tRNA-Ile, tRNA-Met, tRNA-Trp), agrupando 4 diagramas por figura.

> **[SUGESTÃO DE TABELA]**: Tabela com Start/Stop codons dos CDS, similar à do PDF de referência.

> **[SUGESTÃO DE TABELA]**: Tabela com anticódons dos tRNAs.

### 5.6 Compilação dos Entregáveis e Formato GenBank

O módulo COMPILE_SUMMARY reuniu automaticamente todos os produtos da análise em uma pasta organizada com 14 categorias de entregáveis, conforme detalhado na Tabela da Seção 4.8. Destaca-se a geração automática do GenBank Flat File (`.gbk`), que é o formato padrão exigido pelo NCBI para submissão de genomas. O arquivo gerado para *A. leari* contém:

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

### 5.7 Discussão sobre a Robustez da Abordagem

A execução bem-sucedida do pipeline para duas espécies distintas — uma ave (Psittaciformes) e um peixe (Perciformes) —, utilizando datasets de plataformas diferentes (HiSeq X Ten e NovaSeq X Plus), demonstra a generalização da abordagem. A obtenção de genomas circularizados em ambos os casos, com o conjunto completo de 37 genes identificados no caso anotado, valida tanto o fluxo de montagem quanto a integração de ferramentas.

A cobertura de 306× obtida para *A. leari* com apenas 0,17% de leituras mitocondriais em 20 milhões de pares demonstra que o DNA mitocondrial é eficientemente recuperado pelo NOVOPlasty mesmo quando representa uma fração diminuta do total de dados, confirmando as observações da literatura (DIERCKXSENS; MARDULYN; SMITS, 2017).

A lógica iterativa de montagem com re-seeding automático e múltiplos k-mers mostrou-se eficaz, circularizando na primeira tentativa (k-mer 39) para *A. leari* e na primeira tentativa (k-mer 33) para *D. bifasciatus*. Essa robustez é particularmente relevante para cenários em que a seed disponível é de uma espécie distante, ou quando regiões repetitivas dificultam a circularização.

O módulo Pilot QC foi validado por meio de execução isolada do script `pilot_qc.sh` sobre uma subamostra de 100.000 leituras de *A. leari*, extraída do dataset SRR28399504. A análise piloto detectou qualidade Q30 de 86,2%, taxa de adaptador de 66,4% e fração mitocondrial estimada de 0,17% (via correspondência de k-mers com a semente cox1). Com esses valores, a recomendação automática foi de 36 milhões de leituras para atingir 500× de cobertura — resultado coerente com os 306× efetivamente obtidos com 20 milhões de leituras na execução principal. A configuração de produção do perfil `a_leari` não define `sra_max_reads`, ativando o Pilot QC automaticamente em execuções futuras. Essa funcionalidade é particularmente valiosa para pesquisadores trabalhando com espécies pouco caracterizadas, onde não há estimativa prévia da fração mitocondrial.

Em comparação com os trabalhos relacionados, o pipeline desenvolvido apresenta vantagens claras em termos de reprodutibilidade e portabilidade. Enquanto o MitoHiFi (ULIANO-SILVA et al., 2023) utiliza Conda para gerenciamento de dependências, sujeito a variações de ambiente, e o MITObim (HAHN et al., 2013) e o GetOrganelle (JIN et al., 2020) dependem de configuração manual, o pipeline aqui proposto encapsula todas as ferramentas em containers Docker com versões pinadas e orquestra sua execução via Nextflow, garantindo resultados idênticos em qualquer ambiente. A adição da Fase 2 (anotação funcional com MITOS2 e geração automática de entregáveis em formato GenBank) vão além do que a maioria dos pipelines de montagem de organelas oferece, cobrindo o ciclo completo desde o dado bruto até a submissão ao banco de dados.

### 5.8 Relevância Biológica e Científica

Do ponto de vista aplicado, a execução do pipeline possibilitou a reconstrução reprodutível do mitogenoma completo de *Anodorhynchus leari* — espécie endêmica da Caatinga brasileira, restrita ao norte da Bahia, com população estimada em cerca de 1.700 indivíduos em vida livre (ICMBio, 2022), classificada como Em Perigo (EN) pela IUCN. A disponibilização de um genoma mitocondrial completo e anotado fornece subsídios para estudos de genética de populações, filogeografia e conservação, ampliando as ferramentas disponíveis para gestão da biodiversidade.

O mitogenoma de 16.986 bp aqui obtido pode ser comparado diretamente com o de *A. hyacinthinus* (16.999 bp, NC_082165.1), a outra espécie do gênero, permitindo estudos de divergência molecular e análises filogenéticas. A diferença de apenas 13 bp entre os dois genomas sugere uma divergência recente entre as espécies, consistente com estimativas filogeográficas baseadas em marcadores mitocondriais parciais.

A disponibilização pública do pipeline e de seus resultados em repositório GitHub contribui para a democratização do acesso a métodos de bioinformática, permitindo que outros grupos de pesquisa reproduzam as análises ou apliquem o mesmo fluxo a outras espécies de interesse — sejam aves, peixes, invertebrados ou qualquer outro metazoário —, bastando fornecer o acesso SRA e uma sequência semente de gene mitocondrial. Essa característica de generalização, verificada empiricamente pela execução bem-sucedida em dois txons distintos (Psittaciformes e Perciformes), reforça o compromisso com os princípios FAIR (WILKINSON et al., 2016).

---

## Capítulo 6 — Considerações Finais

O presente trabalho teve como objetivo propor e implementar um pipeline de bioinformática voltado à montagem e anotação funcional de genomas mitocondriais, concebido segundo princípios de automação, reprodutibilidade e ciência aberta. A integração de sete ferramentas consolidadas — SRA-Toolkit, FastQC, Trim Galore com Cutadapt, NOVOPlasty, MITOS2, ViennaRNA e Biopython (COCK et al., 2009) — em um fluxo containerizado (Docker) e orquestrado por meio do Nextflow constitui a principal contribuição metodológica desta proposta, ao oferecer uma solução transparente, portável e replicável para análises genômicas.

Do ponto de vista biológico, o pipeline foi aplicado com sucesso à montagem do mitogenoma de *Anodorhynchus leari* (16.986 bp circularizado, 306× de cobertura, 37 genes identificados) e à espécie de validação *Diploprion bifasciatus* (circularizado). Os resultados demonstraram que o pipeline é capaz não apenas de reconstruir genomas mitocondriais completos a partir de dados de leituras curtas, mas também de anotá-los de forma automatizada, gerando todos os entregáveis necessários para apresentação acadêmica e submissão ao GenBank.

Dentre as contribuições originais do pipeline, destacam-se: (i) o módulo Pilot QC, que determina automaticamente o volume de dados necessário para atingir a cobertura desejada, evitando processamento desnecessário; (ii) a lógica iterativa de montagem com re-seeding automático e múltiplos k-mers; (iii) a geração automática de estruturas secundárias de tRNAs e rRNAs via RNAplot/ViennaRNA, funcionalidade anteriormente restrita à interface web do MITOS2; (iv) a conversão automática para GenBank Flat File com tratamento do frameshift do ND3; e (v) a compilação automática de 14 categorias de entregáveis em pasta organizada.

Apesar de seu potencial, algumas limitações precisam ser reconhecidas. O sucesso da montagem depende da qualidade da seed inicial e da profundidade de cobertura dos dados de sequenciamento. A anotação do MITOS2, embora automatizada, pode requerer curadoria manual em casos de estruturas gênicas atípicas (como o frameshift do ND3). Adicionalmente, o tamanho dos containers Docker (~2–4 GB por imagem) pode representar uma barreira em ambientes com espaço de armazenamento limitado.

Como perspectivas futuras, propõem-se: (i) a validação da montagem de *D. bifasciatus* por BLAST contra a referência PZ143763.1; (ii) a submissão do mitogenoma de *A. leari* ao GenBank; (iii) a integração de análises filogenéticas automáticas ao pipeline (MAFFT + RAxML/IQ-TREE); (iv) a geração de gráficos de sintenia e RSCU (Relative Synonymous Codon Usage); e (v) a publicação do pipeline no nf-core como recurso comunitário.

Em síntese, o pipeline desenvolvido representa uma contribuição metodológica significativa para a bioinformática de organelas, ao unir rigor técnico, aplicabilidade biológica e boas práticas de reprodutibilidade computacional. Mais do que atender a um caso de estudo específico, constitui um recurso com potencial de ser reutilizado, adaptado e expandido por diferentes grupos de pesquisa, colaborando para a consolidação de uma ciência mais transparente, padronizada e sustentável.

---

## 8. Referências

ANDERSON, S.; et al. Sequence and organization of the human mitochondrial genome. Nature, v. 290, p. 457-465, 1981.

ANDREWS, S. FastQC: a quality control tool for high throughput sequence data. 2010. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Acesso em: 31 ago. 2025.

ASAKAWA, S.; HIMENO, H.; MIURA, K.; WATANABE, K. Nucleotide sequence and gene organization of the starfish Asterina pectinifera mitochondrial genome. Genetics, v. 140, n. 3, p. 1047-1060, 1995.

BABRAHAM BIOINFORMATICS. Trim Galore! 2019. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. Acesso em: 31 ago. 2025.

BAKER, M. 1,500 scientists lift the lid on reproducibility. Nature, v. 533, n. 7604, p. 452-454, 2016.

BOORE, J. L. Animal mitochondrial genomes. Nucleic Acids Research, v. 27, n. 8, p. 1767-1780, 1999.

BRITISH ECOLOGICAL SOCIETY. A guide to reproducible code in ecology and evolution. London: British Ecological Society, 2017.

CALABRESE, C.; et al. MToolBox: a bioinformatics pipeline for the characterization of human mitochondrial genome variants. Bioinformatics, v. 30, n. 21, p. 3115-3117, 2014. doi: 10.1093/bioinformatics/btu483.

COCK, P. J. A.; et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, v. 25, n. 11, p. 1422-1423, 2009. doi: 10.1093/bioinformatics/btp163.

DIERCKXSENS, N.; MARDULYN, P.; SMITS, G. NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, v. 45, n. 4, p. e18, 2017.

DI TOMMASO, P.; et al. Nextflow enables reproducible computational workflows. Nature Biotechnology, v. 35, n. 4, p. 316-319, 2017.

DONATH, A.; et al. Improved annotation of protein-coding gene boundaries in metazoan mitochondrial genomes. Molecular Ecology Resources, v. 19, n. 4, p. 609-615, 2019. doi: 10.1111/1755-0998.12985.

EKBLOM, R.; WOLF, J. B. W. A field guide to whole-genome sequencing, assembly and annotation. Evolutionary Applications, v. 7, n. 9, p. 1026-1042, 2014.

EWELS, P. A.; et al. The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology, v. 38, n. 3, p. 276-278, 2020.

GRÜNING, B.; et al. Practical computational reproducibility in the life sciences. Cell Systems, v. 6, n. 6, p. 631-635, 2018.

HAHN, C.; BACHMANN, L.; CHEVREUX, B. Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads — a baiting and iterative mapping approach. Nucleic Acids Research, v. 41, n. 13, p. e129, 2013. doi: 10.1093/nar/gkt371.

ICMBio — INSTITUTO CHICO MENDES DE CONSERVAÇÃO DA BIODIVERSIDADE. Plano de Ação Nacional para a Conservação da Arara-azul-de-lear. Brasília: ICMBio, 2022.

IRIDIAN GENOMES. Anodorhynchus leari whole genome sequencing. NCBI Sequence Read Archive, acesso SRR28399504, 2024. Disponível em: https://www.ncbi.nlm.nih.gov/sra/SRR28399504. Acesso em: 10 mar. 2026.

JIN, J.-J.; et al. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology, v. 21, n. 1, p. 241, 2020. doi: 10.1186/s13059-020-02154-5.

KÖSTER, J.; RAHMANN, S. Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, v. 28, n. 19, p. 2520-2522, 2012.

LARMAN, C.; BASILI, V. R. Iterative and incremental development: a brief history. IEEE Computer, v. 36, n. 6, p. 47-56, 2003. doi: 10.1109/MC.2003.1204375.

LEINONEN, R.; SUGAWARA, H.; SHUMWAY, M. The Sequence Read Archive. Nucleic Acids Research, v. 39, supl. 1, p. D19-D21, 2011. doi: 10.1093/nar/gkq1019.

LORENZ, R.; et al. ViennaRNA Package 2.0. Algorithms for Molecular Biology, v. 6, n. 26, 2011. doi: 10.1186/1748-7188-6-26.

MARTIN, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, v. 17, n. 1, p. 10-12, 2011. doi: 10.14806/ej.17.1.200.

MERKEL, D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal, n. 239, p. 2, 2014.

MILLER, J. R.; KOREN, S.; SUTTON, G. Assembly algorithms for next-generation sequencing data. Genomics, v. 95, n. 6, p. 315-327, 2010.

MINDELL, D. P.; SORENSON, M. D.; DIMCHEFF, D. E. An extra nucleotide is not translated in mitochondrial ND3 of some birds and turtles. Molecular Biology and Evolution, v. 15, n. 11, p. 1568-1571, 1998.

NCBI — NATIONAL CENTER FOR BIOTECHNOLOGY INFORMATION. SRA Toolkit Documentation. Bethesda: National Library of Medicine, 2023. Disponível em: https://github.com/ncbi/sra-tools. Acesso em: 10 mar. 2026.

ORELLANA, L. H.; et al. Comparing genomes recovered from time-series metagenomes using long- and short-read sequencing technologies. Microbiome, v. 11, n. 1, p. 105, 2023.

PENG, R. D. Reproducible research in computational science. Science, v. 334, n. 6060, p. 1226-1227, 2011. DOI: https://doi.org/10.1126/science.1213847.

SANDVE, G. K.; et al. Ten simple rules for reproducible computational research. PLoS Computational Biology, v. 9, n. 10, p. e1003285, 2013.

SOMMERVILLE, I. Software Engineering. 10. ed. Harlow: Pearson Education, 2016.

ULIANO-SILVA, M.; et al. MitoHiFi: a pipeline for mitochondrial genome assembly from PacBio HiFi reads. BMC Bioinformatics, v. 24, n. 288, 2023. DOI: 10.1186/s12859-023-05385-y.

WILKINSON, M. D.; et al. The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data, v. 3, p. 160018, 2016.

WILSON, G.; et al. Best practices for scientific computing. PLoS Biology, v. 12, n. 1, e1001745, 2014. doi: 10.1371/journal.pbio.1001745.
