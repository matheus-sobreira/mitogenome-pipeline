# DESENVOLVIMENTO E IMPLEMENTAÇÃO DE UM PIPELINE DE BIOINFORMÁTICA PARA MONTAGEM E ANÁLISE DE MITOGENOMAS

**UNIVERSIDADE DO ESTADO DO RIO GRANDE DO NORTE**
**CAMPUS AVANÇADO DE NATAL**
**DEPARTAMENTO DE CIÊNCIA DA COMPUTAÇÃO**
**BACHARELADO EM CIÊNCIA DA COMPUTAÇÃO**

**MATHEUS SOBREIRA BENEVIDES**

Natal, 2025

---

## Observações

1. O acesso ao projeto foi compartilhado com o professor Wilfredo (Orientador) e o professor Carlos (Professor de TCC). As perguntas aqui serão direcionadas ao professor Wilfredo, porém se o professor Carlos quiser fazer algum comentário, sinta-se à vontade de responder/comentar qualquer coisa.
2. Está no docs temporariamente, e, esse projeto não está padronizado completamente em um modelo acadêmico específico (ABNT ou UERN). Estou deixando desta forma apenas para acelerar a etapa de projeto. Assim que ultrapassada o piloto e após a defesa do projeto, irei começar a transferir para o modelo UERN no overleaf (LaTeX).
3. Faltam imagens das ferramentas utilizadas na metodologia. Essa parte ficará para após a apresentação do projeto.

---

## Sumário

1. Introdução
   - 1.1 Introdução
   - 1.2 Justificativa
   - 1.3 Objetivos
     - 1.3.1 Objetivo Geral
     - 1.3.2 Objetivos Específicos
2. Referencial Teórico
   - 2.1 Montagem de Novo de Genomas e a Análise Mitocondrial
   - 2.2 Reprodutibilidade em Ciência Computacional
   - 2.3 Tecnologia de Conteinerização com Docker
   - 2.4 Sistemas de Gerenciamento de Workflows Científicos
3. Trabalhos Relacionados
   - 3.1 MitoHiFi: Montagem com Leituras Longas
   - 3.2 MITObim: Estratégia de Baiting and Iterative Mapping
   - 3.3 GetOrganelle: Montagem Baseada em Grafos
   - 3.4 MToolBox: Foco em Genomas Mitocondriais Humanos
   - 3.5 NOVOPlasty: Estratégia Seed-and-Extend
4. Metodologia
   - 4.1 Visão Geral do Pipeline
   - 4.2 Ferramentas e Tecnologias Utilizadas
     - 4.2.1 SRA-Toolkit
     - 4.2.2 FastQC
     - 4.2.3 Cutadapt e Trim Galore
     - 4.2.4 NOVOPlasty
     - 4.2.5 Mitos2
     - 4.2.6 tRNAscan-SE
     - 4.2.7 Docker
     - 4.2.8 Nextflow
     - 4.2.9 GitHub
   - 4.3 Aquisição e Preparo dos Dados
   - 4.4 Controle de Qualidade e Pré-processamento
   - 4.5 Montagem do Mitogenoma
   - 4.6 Pós-processamento e Anotação
   - 4.7 Automação e Reprodutibilidade
5. Resultados e Discussão
   - 5.1 Resultados por Etapa do Pipeline
   - 5.2 Discussão sobre a Robustez da Abordagem
   - 5.3 Relevância Biológica e Científica
6. Considerações Finais
7. Referências

---

## 1. Introdução

### 1.1 Introdução

O estudo de genomas mitocondriais ocupa uma posição central na biologia evolutiva, genética de populações e sistemática, uma vez que essas moléculas oferecem informações altamente informativas sobre parentesco, variação intra e interespecífica e processos de conservação. O DNA mitocondrial, por suas características estruturais e funcionais, como a herança predominantemente materna, a ausência de recombinação e a elevada taxa evolutiva em comparação ao DNA nuclear, tornou-se um marcador molecular amplamente utilizado em diferentes contextos biológicos (BOORE, 1999; EKBLOM; WOLF, 2014).

O avanço das tecnologias de sequenciamento de nova geração (Next-Generation Sequencing — NGS) ampliou o acesso a grandes volumes de dados, favorecendo a obtenção de mitogenomas completos e de alta qualidade. No entanto, a reconstrução dessas moléculas a partir de dados brutos permanece um desafio técnico, especialmente no caso de leituras curtas, que exigem algoritmos sofisticados para superar problemas como regiões repetitivas e possíveis falhas de montagem. Nesse cenário, ferramentas especializadas, como o NOVOPlasty, têm desempenhado papel fundamental ao implementar estratégias de *seed-and-extend* capazes de recuperar organelas de forma eficiente a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017).

Embora a diversidade de softwares tenha contribuído para avanços significativos na montagem de organelas, a ausência de padronização metodológica ainda compromete a reprodutibilidade das análises. A multiplicidade de versões, dependências e configurações de ambiente gera obstáculos práticos para replicar experimentos, fenômeno amplamente conhecido como *dependency hell* (GRÜNING et al., 2018). Esse problema se insere no contexto mais amplo da crise de reprodutibilidade científica (BAKER, 2016), que tem motivado a adoção, na bioinformática, de práticas da engenharia de software, como a conteinerização de ambientes computacionais com Docker (MERKEL, 2014) e a orquestração de fluxos com sistemas de gerenciamento de workflows, a exemplo do Nextflow (DI TOMMASO et al., 2017).

Dessa forma, a presente pesquisa propõe-se a contribuir para esse campo ao desenvolver um pipeline automatizado, reprodutível e portável para a montagem e anotação de mitogenomas. A integração de ferramentas bioinformáticas consolidadas em um ecossistema containerizado e orquestrado garante não apenas consistência metodológica, mas também acessibilidade e transparência, em consonância com os princípios de ciência aberta e os critérios FAIR (*Findable, Accessible, Interoperable, Reusable*) (WILKINSON et al., 2016). Essa proposta estabelece a base para a justificativa do trabalho, apresentada na seção seguinte, onde são discutidos em maior detalhe os desafios técnicos e científicos que motivam sua realização.

### 1.2 Justificativa

A análise de genomas mitocondriais (mitogenomas) representa um componente essencial para a biologia evolutiva, a genética de populações e a sistemática, dado que essas moléculas fornecem informações altamente informativas sobre relações de parentesco, variação intra e interespecífica e processos de conservação. Sua relevância científica foi consolidada em trabalhos seminais que destacaram o papel estrutural e evolutivo dos mitogenomas como marcadores moleculares (BOORE, 1999). Em seguida, estudos enfatizaram sua aplicabilidade em investigações de genética de populações e evolução molecular (EKBLOM; WOLF, 2014), enquanto contribuições mais recentes demonstraram seu potencial em análises de biodiversidade e conservação por meio de avanços metodológicos em montagem e anotação de organelas (ULIANO-SILVA et al., 2023).

Embora a obtenção de mitogenomas completos tenha sido amplamente favorecida pelo advento do Sequenciamento de Nova Geração (NGS), a execução de pipelines de montagem ainda enfrenta desafios técnicos significativos. A diversidade de ferramentas, versões e dependências envolvidas gera um cenário de difícil replicação, conhecido como *dependency hell*, que compromete a reprodutibilidade e a padronização das análises (SANDVE et al., 2013; BAKER, 2016; GRÜNING et al., 2018).

Esse problema se torna especialmente relevante quando aplicado ao estudo de espécies ameaçadas de extinção, nas quais a disponibilidade de amostras biológicas é limitada e a confiabilidade dos resultados genômicos é crucial para fundamentar decisões de conservação. É o caso da arara-azul-de-lear (*Anodorhynchus leari*), ave endêmica da Caatinga nordestina do Brasil, classificada como "Em Perigo" pela IUCN Red List. Com população estimada em pouco mais de 1.700 indivíduos na natureza, a espécie depende de dados genômicos acessíveis e validáveis para subsidiar programas de manejo, reprodução em cativeiro e monitoramento populacional. Embora seu mitogenoma já tenha sido previamente reconstruído em atividade acadêmica, essa montagem foi realizada de forma manual, sem padronização metodológica e sem disponibilização pública do fluxo de trabalho — o que impede sua verificação e replicação por outros pesquisadores.

Nesse contexto, torna-se necessário o desenvolvimento de soluções que conciliem robustez biológica e consistência computacional. A adoção de práticas modernas da engenharia de software, como a containerização e os sistemas de gerenciamento de workflows, tem se mostrado fundamental para superar limitações de reprodutibilidade, além de permitir que fluxos analíticos sejam portáveis, escaláveis e transparentes (DI TOMMASO et al., 2017; EWELS et al., 2020).

Assim, este trabalho justifica-se pela proposta de construir um pipeline automatizado e reprodutível para a montagem de mitogenomas, unindo boas práticas de ciência aberta e engenharia computacional à relevância biológica desses genomas. A aplicação ao caso concreto de *A. leari* demonstra não apenas a viabilidade técnica da abordagem, mas também sua contribuição direta para a conservação de uma espécie ameaçada da biodiversidade brasileira.

### 1.3 Objetivos

#### 1.3.1 Objetivo Geral

Desenvolver e validar um pipeline de bioinformática automatizado, reprodutível e portável para a montagem *de novo* de genomas mitocondriais, utilizando tecnologias de conteinerização (Docker) e sistemas de gerenciamento de workflows.

#### 1.3.2 Objetivos Específicos

a) Containerizar, por meio de Dockerfiles, cada uma das ferramentas bioinformáticas necessárias ao pipeline de montagem de mitogenomas (SRA-Toolkit, FastQC, ferramentas de *trimming*, NOVOPlasty e ferramentas de anotação).
b) Desenvolver o script de orquestração do pipeline utilizando um gerenciador de workflow (ex.: Nextflow), definindo os processos, suas entradas, saídas e interdependências.
c) Integrar os containers Docker ao workflow, garantindo que cada etapa da análise seja executada em seu próprio ambiente isolado e pré-configurado.
d) Validar o pipeline automatizado por meio de um estudo de caso com dados públicos de sequenciamento, comparando o genoma montado e anotado com resultados previamente publicados.
e) Disponibilizar o pipeline completo e sua documentação em repositório público de controle de versão (ex.: GitHub), assegurando sua acessibilidade e a reprodutibilidade na comunidade científica.

Dessa forma, este capítulo estabeleceu o contexto biológico e computacional do trabalho, destacando a relevância dos genomas mitocondriais, os desafios de reprodutibilidade na bioinformática e a proposta de um pipeline automatizado e portável. O Capítulo 2 apresenta, a seguir, os fundamentos teóricos necessários para compreender em maior profundidade os conceitos que sustentam a execução do pipeline.

---

## 2. Referencial Teórico

O presente capítulo reúne os fundamentos conceituais necessários para compreender a proposta metodológica deste trabalho. São abordados, inicialmente, os princípios da montagem *de novo* de genomas e as especificidades do DNA mitocondrial, que justificam sua relevância como objeto de estudo. Em seguida, discute-se a questão da reprodutibilidade em ciência computacional e as soluções contemporâneas oferecidas pela engenharia de software, com destaque para a tecnologia de conteinerização. Por fim, apresentam-se os sistemas de gerenciamento de workflows científicos, cuja integração com containers constitui a base para a construção de pipelines modernos, portáveis e reprodutíveis.

### 2.1 Montagem *De Novo* de Genomas e a Análise Mitocondrial

Em bioinformática, a montagem de genomas consiste no processo de reconstrução da sequência original de DNA a partir de fragmentos obtidos por tecnologias de Sequenciamento de Nova Geração (NGS). O método mais amplamente utilizado para gerar esses fragmentos é o *shotgun sequencing*, no qual o DNA genômico é clivado aleatoriamente em milhões de pequenos pedaços que, posteriormente, são sequenciados. Essa estratégia garante alta cobertura do genoma, mas resulta em dados fragmentados que precisam ser cuidadosamente alinhados e sobrepostos para permitir a reconstrução da sequência completa (EKBLOM; WOLF, 2014).

**Figura 1** — Esquema do processo de *shotgun sequencing*.
*Fonte: EKBLOM; WOLF (2014).*

Na montagem genômica, fragmentos sobrepostos são inicialmente agrupados em *contigs*, sequências contínuas que representam trechos do genoma. Quando existe informação adicional que permite inferir a ordem relativa entre diferentes contigs — como dados de leituras pareadas (*paired-end* ou *mate-pair*) —, esses podem ser organizados em estruturas mais amplas chamadas *scaffolds*. Assim, enquanto os contigs representam sequências consolidadas diretamente a partir das leituras, os scaffolds constituem hipóteses mais abrangentes sobre a organização genômica, servindo como passo intermediário até a obtenção de montagens completas e anotadas.

As tecnologias de NGS diferem substancialmente no comprimento das leituras produzidas. Sequenciadores Illumina, por exemplo, geram leituras curtas de alta acurácia (50 a 300 pares de bases), que exigem algoritmos sofisticados para lidar com regiões repetitivas. Em contrapartida, tecnologias como PacBio e Oxford Nanopore produzem leituras longas, que podem alcançar milhares de pares de bases em uma única sequência, favorecendo a continuidade da montagem, ainda que a custos mais elevados e com menor rendimento global (EKBLOM; WOLF, 2014).

**Figura 2** — Comparação entre montagens SR e LR.
*Fonte: ORELLANA et al. (2023).*

A Figura 2 compara *metagenome-assembled genomes* (MAGs) recuperados a partir de tecnologias de *short reads* (Illumina) e *long reads* (PacBio HiFi). Genomas derivados de *long reads* apresentam maior continuidade (N50 mais elevado), menor fragmentação (menos contigs) e maior número de genes preditos em comparação aos derivados de *short reads* (ORELLANA et al., 2023) — evidenciando como a escolha da tecnologia influencia diretamente a qualidade da montagem. Enquanto os *short reads* permitem recuperar um número maior de genomas devido à profundidade de sequenciamento, eles produzem montagens mais fragmentadas. Já os *long reads* oferecem menor diversidade detectada, mas geram montagens mais contínuas e completas, preservando genes essenciais, como o 16S rRNA, e representando melhor regiões repetitivas.

O DNA mitocondrial (mtDNA) reúne um conjunto de particularidades que o tornam um alvo privilegiado para a montagem *de novo*. Em metazoários, ele corresponde a uma molécula circular de fita dupla, com tamanho geralmente entre 12 e 22 mil pares de bases, codificando 37 genes distribuídos em 13 proteínas associadas à fosforilação oxidativa, 22 RNAs transportadores (tRNAs) e 2 RNAs ribossômicos (rRNAs) (ASAKAWA et al., 1995; BOORE, 1999). A Figura 3, a seguir, apresenta a representação circular típica de um genoma mitocondrial de metazoários, evidenciando a organização gênica compacta, incluindo genes codificadores de proteínas, tRNAs e rRNAs.

**Figura 3** — Representação circular do genoma mitocondrial.
*Fonte: Do próprio autor (2024).*

A ausência de íntrons, a estrutura compacta e a elevada taxa de cópias por célula tornam o mtDNA especialmente acessível, facilitando sua recuperação mesmo em experimentos voltados ao genoma nuclear, nos quais ele aparece como subproduto. Além disso, características funcionais como a herança predominantemente materna, a ausência geral de recombinação e a taxa evolutiva superior à do DNA nuclear ampliam sua utilidade em estudos biológicos. Essas propriedades explicam sua ampla aplicação em análises de genética de populações, filogeografia, sistemática e conservação da biodiversidade.

**Figura 4** — Exemplo de arquivo FASTA.
*Fonte: Do próprio autor (2024).*

A Figura 4 representa um trecho de arquivo FASTA do genoma mitocondrial da *Anodorhynchus leari*, representando a forma textual dos dados brutos utilizados como entrada em pipelines de montagem genômica. Esse exemplo evidencia como os dados genômicos são inicialmente armazenados como cadeias de nucleotídeos em arquivos de texto simples, que posteriormente são processados e reconstruídos em estruturas circulares biologicamente coerentes. A contraposição entre o dado cru (FASTA) e a representação organizada (Figura 3) demonstra a importância dos algoritmos de montagem para transformar informação bruta em conhecimento biológico.

Para explorar essas vantagens, foram desenvolvidos algoritmos específicos de montagem. Uma das abordagens mais empregadas é a *seed-and-extend*, em que uma sequência inicial conhecida (*seed*) — que pode ser um gene ou até mesmo um fragmento de organela de espécie próxima — serve como ponto de partida para estender iterativamente a montagem até que a molécula circular seja reconstruída. O NOVOPlasty é um dos principais programas baseados nessa estratégia, destacando-se pela eficiência na recuperação de genomas mitocondriais e cloroplastidiais a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017). Avanços recentes, como o pipeline MitoHiFi, exploram leituras longas de alta fidelidade, evidenciando a constante evolução das ferramentas voltadas à montagem de organelas (ULIANO-SILVA et al., 2023).

### 2.2 Reprodutibilidade em Ciência Computacional

A reprodutibilidade é um pilar do método científico, consistindo na capacidade de um pesquisador independente replicar um experimento, com os mesmos materiais e métodos, e obter resultados consistentes. Em pesquisas computacionais, isso se traduz na possibilidade de executar o mesmo código e software sobre os mesmos dados e chegar ao mesmo resultado (PENG, 2011; SANDVE et al., 2013).

Contudo, a ciência moderna enfrenta o que tem sido amplamente denominado como uma "crise de reprodutibilidade", na qual um número significativo de estudos publicados se mostra difícil ou impossível de ser replicado (BAKER, 2016). No domínio da bioinformática, essa crise é particularmente acentuada devido à alta complexidade dos pipelines de análise. Um único pipeline pode envolver dezenas de ferramentas de software distintas, cada uma desenvolvida por diferentes grupos, em diferentes linguagens de programação e com um conjunto único e, por vezes, conflitante de dependências de bibliotecas e do próprio sistema operacional. A tarefa de instalar e configurar manualmente esse ecossistema de software para replicar uma análise é extremamente suscetível a erros, uma condição frequentemente descrita como *dependency hell* (inferno de dependências) (GRÜNING et al., 2018).

Pequenas variações na versão de uma ferramenta ou de uma biblioteca subjacente podem levar a resultados drasticamente diferentes, comprometendo a validade e a confiabilidade da pesquisa. Essa barreira técnica não apenas dificulta a verificação dos resultados, mas também impede a reutilização e a adaptação de métodos computacionais, contrariando os princípios FAIR (*Findable, Accessible, Interoperable, and Reusable*), que visam maximizar o valor dos dados e das análises científicas (SANDVE et al., 2013; WILKINSON et al., 2016; GRÜNING et al., 2018; BRITISH ECOLOGICAL SOCIETY, 2017).

### 2.3 Tecnologia de Conteinerização com Docker

Como solução para o problema do *dependency hell* e para garantir a reprodutibilidade computacional, a comunidade científica tem adotado práticas e tecnologias da engenharia de software, com destaque para a virtualização em nível de sistema operacional, também conhecida como containerização. A plataforma Docker se estabeleceu como a principal ferramenta desse paradigma, permitindo o empacotamento de uma aplicação e de todo o seu ambiente de execução — incluindo bibliotecas, códigos-fonte e dependências — em uma unidade padronizada e isolada chamada *container* (MERKEL, 2014; GRÜNING et al., 2018).

O funcionamento do Docker baseia-se em dois conceitos centrais: a imagem e o container. A imagem é um template estático e imutável, construído a partir de um arquivo de texto com instruções sequenciais chamado *Dockerfile*. Esse arquivo serve como uma "receita" que define o sistema operacional base, as dependências a serem instaladas e os comandos a serem executados, criando um ambiente de software completo e autocontido. O *container*, por sua vez, é uma instância executável e isolada de uma imagem. Ao contrário de máquinas virtuais tradicionais, os containers compartilham o kernel do sistema operacional hospedeiro, o que os torna leves e eficientes para iniciar (BRITISH ECOLOGICAL SOCIETY, 2017).

Essa arquitetura garante que um software containerizado se comporte de maneira idêntica em qualquer ambiente que possua o Docker instalado — seja em computadores pessoais, servidores de produção ou ambientes de nuvem — resolvendo de forma eficaz o problema da inconsistência de ambientes e sendo um passo fundamental para a criação de pipelines científicos portáveis e reprodutíveis (DI TOMMASO et al., 2017).

### 2.4 Sistemas de Gerenciamento de *Workflows* Científicos

A containerização resolve a questão do ambiente de software, mas não a orquestração das múltiplas etapas que compõem uma análise bioinformática. A execução manual e sequencial de cada container é ineficiente e propensa a erros. Para gerenciar essa complexidade, foram desenvolvidos Sistemas de Gerenciamento de Workflows (*Workflow Management Systems*), como o Nextflow (DI TOMMASO et al., 2017) e o Snakemake (KÖSTER; RAHMANN, 2012). Essas plataformas permitem que os pesquisadores definam pipelines computacionais complexos por meio de linguagens de alto nível.

Nesses sistemas, um workflow é decomposto em uma série de processos ou regras interconectadas. Cada processo define uma tarefa computacional discreta, especificando seus comandos, suas entradas (arquivos ou variáveis) e suas saídas. O gerenciador de workflow é então responsável por monitorar o estado dos dados e executar os processos na ordem correta, automatizando o fluxo de informação entre as etapas.

A principal vantagem dessas ferramentas reside na abstração da camada de execução. Elas integram-se nativamente com tecnologias como Docker, permitindo que cada processo seja executado em seu próprio container pré-configurado. Além disso, são compatíveis com diversos ambientes de execução, desde computadores pessoais até clusters de computação de alto desempenho (HPC) e plataformas de nuvem. Isso confere ao pipeline características essenciais para a ciência moderna: portabilidade, escalabilidade e retomada automática (*resumability*) em caso de falhas, evitando o reprocessamento de etapas já concluídas.

A combinação de um gerenciador de workflows com a containerização é hoje considerada o padrão-ouro para o desenvolvimento de pipelines de bioinformática, sendo a base de iniciativas comunitárias como o nf-core, que busca padronizar e compartilhar pipelines que seguem as melhores práticas de reprodutibilidade (EWELS et al., 2020).

O referencial teórico expôs os conceitos fundamentais que dão suporte à proposta deste trabalho, abordando desde os princípios de montagem *de novo* de genomas até as soluções contemporâneas de reprodutibilidade, como a containerização e os sistemas de gerenciamento de workflows. Esses conceitos fornecem a base conceitual sobre a qual se apoiam os trabalhos relacionados, discutidos no Capítulo 3, que permitem situar a presente proposta no estado da arte da bioinformática aplicada à montagem de mitogenomas.

---

## 3. Trabalhos Relacionados

Nesta seção, são apresentados e analisados trabalhos recentes da literatura que, assim como este projeto, propõem soluções computacionais para a montagem de genomas de organelas. O objetivo desta análise comparativa é contextualizar a presente proposta no estado da arte da área, identificando as abordagens metodológicas existentes e destacando os diferenciais e a contribuição original do pipeline a ser desenvolvido neste trabalho.

### 3.1 MitoHiFi: Montagem com Leituras Longas

Um trabalho de destaque na área é o MitoHiFi, um pipeline desenvolvido para a montagem de genomas mitocondriais a partir de dados de leituras longas de alta fidelidade (PacBio HiFi). Proposto por Uliano-Silva et al. (2023), o objetivo principal era criar um método que aproveitasse a extensão e a precisão deste tipo de dado para gerar mitogenomas completos e corretamente circularizados com maior eficiência. A metodologia do MitoHiFi integra ferramentas como o hifiasm-meta para montagem e o GetOrganelle para a etapa de finalização, sendo orquestrada por meio de scripts em Shell e utilizando Conda para o gerenciamento das dependências de software.

**Figura 4** — Página oficial do repositório MitoHiFi no GitHub.
*Fonte: Uliano-Silva et al. (2023).*

A abordagem do MitoHiFi é extremamente atual por focar em uma tecnologia de sequenciamento de ponta. No entanto, sua implementação via scripts e Conda, embora eficaz para seu propósito, pode apresentar limitações de portabilidade e reprodutibilidade quando comparada a soluções baseadas em containers e sistemas de workflow. O diferencial do presente TCC reside precisamente neste ponto: enquanto o MitoHiFi inova na aplicação de um novo tipo de dado biológico, nossa proposta foca em uma arquitetura de software mais robusta, utilizando Docker e Nextflow, para garantir que o pipeline de montagem (baseado em leituras curtas) seja universalmente reprodutível e portável.

### 3.2 MITObim: Estratégia de *Baiting and Iterative Mapping*

O MITObim, desenvolvido por Hahn et al. (2013), foi um dos primeiros métodos amplamente utilizados para montagem de organelas a partir de dados de NGS. Sua estratégia de *baiting and iterative mapping* consiste em iniciar o processo com uma sequência guia (um fragmento do genoma alvo ou de uma espécie próxima) e, a partir dela, recuperar leituras homólogas, refinando o alinhamento de forma iterativa até a montagem completa do mitogenoma.

**Figura 5** — Página oficial do repositório MITObim no GitHub.
*Fonte: Hahn et al. (2013).*

Embora tenha sido pioneiro na aplicação desta abordagem, o MITObim apresenta limitações em termos de eficiência computacional e escalabilidade, especialmente quando aplicado a conjuntos de dados maiores ou de maior complexidade. Apesar de sua relevância histórica, vem sendo gradualmente substituído por ferramentas mais modernas, como o NOVOPlasty e o GetOrganelle.

### 3.3 GetOrganelle: Montagem Baseada em Grafos

O GetOrganelle, proposto por Jin et al. (2020), é uma das ferramentas mais recentes e populares para montagem de organelas. Ele utiliza grafos de montagem derivados de dados de NGS para isolar e reconstruir genomas de organelas, oferecendo resultados robustos tanto para mitocôndrias quanto para cloroplastos. A ferramenta tem se destacado por sua precisão na recuperação de sequências circulares e por oferecer resultados consistentes mesmo em organismos com maior complexidade genômica.

**Figura 6** — Página oficial do repositório GetOrganelle no GitHub.
*Fonte: Jin et al. (2020).*

No entanto, assim como outras ferramentas tradicionais, o GetOrganelle depende de um ambiente de software adequado para ser executado, geralmente configurado manualmente pelo usuário ou via gerenciadores como Conda, o que ainda pode limitar sua reprodutibilidade plena em diferentes contextos computacionais.

### 3.4 MToolBox: Foco em Genomas Mitocondriais Humanos

O MToolBox, apresentado por Calabrese et al. (2014), foi uma das primeiras pipelines integradas voltadas especificamente para o estudo de genomas mitocondriais humanos. Ele combina montagem, anotação e análise funcional em um único fluxo, permitindo identificar variantes mitocondriais relevantes para estudos biomédicos e populacionais.

Apesar de sua contribuição para o campo, o MToolBox é fortemente orientado ao contexto humano e não é tão amplamente aplicável a outros organismos. Além disso, por ter sido desenvolvido antes da adoção massiva de containers e workflows, carece de recursos de portabilidade e reprodutibilidade presentes em soluções mais modernas.

**Figura 7** — Página oficial do repositório MToolBox no GitHub.
*Fonte: Calabrese et al. (2014).*

### 3.5 NOVOPlasty: Estratégia *Seed-and-Extend*

O NOVOPlasty, descrito por Dierckxsens et al. (2017), é atualmente uma das ferramentas mais utilizadas para montagem de genomas mitocondriais e cloroplastidiais a partir de dados de NGS de leituras curtas. Baseado na estratégia *seed-and-extend*, o programa inicia a montagem a partir de uma sequência semente fornecida pelo usuário e estende iterativamente o genoma até que a molécula circular esteja completa. Sua popularidade decorre da facilidade de uso e da capacidade de gerar montagens confiáveis mesmo a partir de dados de genoma total.

No entanto, o NOVOPlasty, como software isolado, ainda demanda configuração manual do ambiente e execução direta pelo usuário, o que pode dificultar sua integração em fluxos mais complexos e comprometer a reprodutibilidade quando diferentes versões ou configurações são utilizadas. É justamente neste ponto que se insere o presente TCC: ao invés de propor um novo montador, o objetivo é integrar o NOVOPlasty em um ecossistema de containers e workflow (Docker + Nextflow), assegurando que seu uso seja portável, automatizado e reprodutível em diferentes ambientes de execução.

**Figura 8** — Página oficial do repositório NOVOPlasty no GitHub.
*Fonte: Dierckxsens et al. (2017).*

A análise dos trabalhos relacionados permitiu identificar diferentes abordagens para a montagem de organelas, incluindo pipelines baseados em leituras longas, estratégias iterativas e montadores especializados como o NOVOPlasty. Embora cada solução apresente contribuições relevantes, permanecem lacunas quanto à padronização, portabilidade e reprodutibilidade dos fluxos de análise. Nesse contexto, o Capítulo 4 descreve a metodologia adotada neste trabalho, detalhando a concepção do pipeline proposto e as ferramentas que o compõem.

---

## 4. Metodologia

Nesta seção, é apresentado o plano metodológico para o desenvolvimento do pipeline de montagem de mitogenomas. A abordagem combina etapas de análise de bioinformática com práticas de engenharia de software para garantir a automação, a reprodutibilidade e a portabilidade da solução final.

### 4.1 Visão Geral do *Pipeline*

O pipeline proposto neste trabalho foi concebido para realizar a montagem de genomas mitocondriais a partir de dados de sequenciamento de nova geração (NGS), integrando ferramentas bioinformáticas consolidadas em um ambiente de execução automatizado, reprodutível e portável. A arquitetura segue uma estrutura modular, em que cada etapa corresponde a uma fase distinta do processamento dos dados, sendo encapsulada em *containers* Docker para garantir consistência na execução (MERKEL, 2014) e integrada por meio do gerenciador de workflows Nextflow (DI TOMMASO et al., 2017).

De forma geral, o pipeline é composto por cinco macroetapas principais. A primeira é a aquisição e preparo dos dados, responsável por obter *datasets* públicos em repositórios como o Sequence Read Archive (SRA) do NCBI, uma das principais bases de dados de sequenciamento disponíveis para a comunidade científica (LEINONEN; SUGAWARA; SHUMWAY, 2011). Em seguida, realiza-se o controle de qualidade e pré-processamento, que envolve a avaliação das leituras por meio do FastQC (ANDREWS, 2010) e a remoção de adaptadores e bases de baixa qualidade utilizando ferramentas como Cutadapt (MARTIN, 2011) ou Trim Galore (BABRAHAM BIOINFORMATICS, 2019).

A terceira etapa corresponde à montagem do mitogenoma, realizada com o NOVOPlasty, ferramenta especializada que utiliza a estratégia *seed-and-extend* para reconstruir a molécula circular completa a partir de leituras curtas (DIERCKXSENS; MARDULYN; SMITS, 2017). A quarta etapa, de pós-processamento e anotação, avalia a qualidade da montagem e realiza a anotação funcional do genoma, utilizando ferramentas como MITOS2 (DONATH et al., 2019) e tRNAscan-SE (LOWE; CHAN, 2016).

Por fim, a quinta etapa concentra-se na automação e reprodutibilidade, em que cada componente do pipeline é containerizado com Docker (MERKEL, 2014) e orquestrado via Nextflow (DI TOMMASO et al., 2017). Esse padrão metodológico segue práticas amplamente adotadas em iniciativas comunitárias, como o projeto nf-core, que disponibiliza workflows bioinformáticos containerizados e reprodutíveis (EWELS et al., 2020).

**Figura 9** — Visão geral do Pipeline proposto.
*Fonte: Do próprio autor (2025).*

A Figura 9 apresenta a visão geral do pipeline, representada em formato de fluxograma, ilustrando o encadeamento das etapas e destacando a integração entre análise bioinformática e práticas de engenharia de software.

### 4.2 Ferramentas e Tecnologias Utilizadas

A implementação de um pipeline bioinformático demanda não apenas a escolha de algoritmos adequados, mas também a definição criteriosa de ferramentas que sejam estáveis, bem documentadas e amplamente utilizadas pela comunidade científica. Nesse contexto, a seleção de softwares para este trabalho foi orientada por três princípios fundamentais: (i) funcionalidade biológica, garantindo que cada ferramenta seja apropriada para a etapa que executa, desde a aquisição dos dados até a anotação final do genoma; (ii) confiabilidade e validação prévia, privilegiando programas consolidados na literatura e em uso corrente em pesquisas genômicas; e (iii) compatibilidade com práticas modernas de engenharia de software, assegurando que cada aplicação possa ser encapsulada em containers e integrada a sistemas de gerenciamento de workflows.

Além de atender a esses critérios, as ferramentas escolhidas refletem o estado da arte em bioinformática aplicada ao estudo de organelas, cobrindo desde utilitários básicos, como o SRA-Toolkit para acesso a bancos de dados, até plataformas mais complexas, como o NOVOPlasty para montagem especializada de mitogenomas. Da mesma forma, tecnologias transversais como Docker e Nextflow foram incorporadas não pelo papel biológico que desempenham, mas pela sua relevância em assegurar portabilidade, escalabilidade e reprodutibilidade, elementos indispensáveis para a robustez metodológica do pipeline.

Assim, a seguir neste capítulo é descrito de forma integrada as principais ferramentas e tecnologias que compõem o pipeline de montagem de genomas mitocondriais. Cada software, desde os utilitários de aquisição de dados até as plataformas de anotação, foi selecionado com base em sua relevância científica, confiabilidade e compatibilidade com práticas modernas de reprodutibilidade computacional. A incorporação de tecnologias como Docker, Nextflow e GitHub complementa esse conjunto, garantindo que a execução do pipeline ocorra de maneira padronizada, portável e transparente, em consonância com as melhores práticas da bioinformática contemporânea.

#### 4.2.1 SRA-Toolkit

O SRA-Toolkit é o pacote oficial de utilitários desenvolvido e mantido pelo National Center for Biotechnology Information (NCBI) para acesso, download e manipulação de dados depositados no Sequence Read Archive (SRA), um dos maiores repositórios públicos de dados de sequenciamento do mundo (LEINONEN; SUGAWARA; SHUMWAY, 2011).

Do ponto de vista técnico, o SRA-Toolkit é composto por um conjunto de ferramentas de linha de comando com funções especializadas. As duas principais utilizadas neste pipeline são:

1. **prefetch** — responsável pelo download do arquivo binário no formato `.sra`, que corresponde à representação comprimida e proprietária dos dados de sequenciamento armazenados no NCBI. O download via prefetch é mais estável e eficiente do que o acesso direto via streaming, especialmente para datasets de grande volume.

2. **fasterq-dump** — ferramenta de conversão do formato `.sra` para o formato FASTQ, padrão na bioinformática moderna por reunir em um mesmo arquivo a sequência nucleotídica e as informações de qualidade (*scores* Phred) de cada base. O fasterq-dump é a versão otimizada e multi-threaded do antigo fastq-dump, oferecendo velocidades de conversão significativamente superiores.

Para dados *paired-end*, o fasterq-dump produz automaticamente dois arquivos FASTQ separados — um para as leituras *forward* (R1) e outro para as *reverse* (R2) — preservando a correspondência entre os pares de leituras. No caso do dataset de *A. leari* (SRR28399504), a conversão do arquivo `.sra` de 11 GB resultou em aproximadamente 87 GB de dados FASTQ descomprimidos (43,5 GB por arquivo), contendo 118,5 milhões de pares de leituras.

Uma consideração prática importante é que o fasterq-dump, diferentemente de seu antecessor fastq-dump, não suporta a flag `-X` para limitação do número de reads durante a conversão. Para viabilizar a truncagem do dataset (necessária para reduzir o custo computacional sem comprometer a montagem do mitogenoma), o pipeline implementa uma estratégia de pós-conversão: após a geração completa dos arquivos FASTQ, o comando `head` é utilizado para reter apenas os primeiros N reads, descartando o excedente. Essa abordagem é detalhada na Seção 4.3.

Entre as vantagens do SRA-Toolkit destacam-se: (i) a integração nativa com o banco SRA, garantindo acesso direto e confiável aos dados; (ii) a verificação automática de integridade dos arquivos baixados via checksums; e (iii) a ampla compatibilidade com formatos de saída utilizados em ferramentas subsequentes do pipeline. Como limitação, pode-se citar o volume de armazenamento temporário necessário durante a conversão, que pode alcançar três vezes o tamanho do arquivo `.sra` original.

A escolha pelo SRA-Toolkit neste trabalho se justifica por ser a ferramenta oficial e mais amplamente validada para acesso a dados públicos de sequenciamento, assegurando que a etapa de aquisição de dados seja realizada de forma padronizada, rastreável e compatível com os princípios de ciência aberta.

#### 4.2.2 FastQC

O FastQC é uma das ferramentas mais utilizadas na etapa de controle de qualidade de dados de sequenciamento de nova geração (NGS). Desenvolvido pelo Babraham Institute, tornou-se referência para avaliação inicial de leituras em formato FASTQ, sendo amplamente empregado em pipelines de bioinformática devido à sua rapidez, robustez e facilidade de interpretação (ANDREWS, 2010).

Do ponto de vista técnico, o FastQC processa os arquivos de entrada em formato FASTQ e gera como saída relatórios em HTML interativo e arquivos compactados (.zip), contendo gráficos e estatísticas sobre diferentes aspectos das sequências. Entre os módulos mais relevantes estão:

1. *Per base sequence quality*, que avalia a distribuição da qualidade das bases em cada posição do read;
2. *Per sequence quality scores*, que identifica leituras com baixa qualidade geral;
3. *Per base GC content*, que compara a distribuição de GC observada com a esperada;
4. *Overrepresented sequences*, que detecta possíveis contaminantes ou adaptadores;
5. *Sequence duplication levels*, que informa o grau de redundância entre leituras;
6. *Adapter content*, que indica a presença de adaptadores de sequenciamento não removidos.

Entre as principais vantagens do FastQC estão: (i) a rapidez de execução mesmo em grandes conjuntos de dados, (ii) a geração de relatórios de fácil interpretação, e (iii) a capacidade de identificar múltiplos problemas em uma única execução. No entanto, o programa possui limitações: por ser uma análise pré-alinhamento, não oferece informações sobre a qualidade do mapeamento ou sobre a distribuição de leituras no genoma de referência. Além disso, os resultados devem ser interpretados com cautela, já que pequenas variações em métricas, como o conteúdo GC, podem refletir características biológicas legítimas e não necessariamente artefatos técnicos.

A escolha pelo FastQC neste trabalho se justifica por sua ampla aceitação na comunidade científica como ferramenta de triagem inicial de qualidade. Sua aplicação permite identificar rapidamente problemas técnicos, orientar a necessidade de *trimming* ou filtragem, e garantir que apenas leituras de qualidade adequada sejam utilizadas nas etapas de montagem do mitogenoma. Dessa forma, o FastQC atua como um ponto de controle essencial para assegurar a confiabilidade das análises subsequentes.

#### 4.2.3 Cutadapt e Trim Galore

Após a avaliação inicial da qualidade das leituras, pode ser necessário realizar etapas de filtragem para remover adaptadores residuais e bases de baixa qualidade. Nesse contexto, duas ferramentas complementares são utilizadas no pipeline: Cutadapt e Trim Galore.

O Cutadapt é um programa desenvolvido por Martin (2011) especificamente para identificar e remover sequências de adaptadores em dados de alto rendimento. As leituras geradas por tecnologias como Illumina frequentemente contêm fragmentos de adaptadores ligados durante a preparação da biblioteca, que, se não forem removidos, podem comprometer a montagem ou o mapeamento subsequente. O Cutadapt utiliza algoritmos de busca eficiente para detectar esses fragmentos nas extremidades das leituras e removê-los, além de permitir o corte de bases de baixa qualidade de acordo com pontuações Phred. Como entrada, o programa recebe arquivos FASTQ, e como saída gera novos FASTQ já filtrados, prontos para análises subsequentes.

O Trim Galore é um *wrapper* que automatiza a execução do Cutadapt em conjunto com o FastQC, simplificando o processo de *trimming* e tornando-o mais acessível ao usuário (BABRAHAM BIOINFORMATICS, 2019). Com ele, é possível remover adaptadores mesmo quando a sequência exata não é previamente conhecida, o que é particularmente útil em experimentos com bibliotecas de origem diversa. Além disso, o Trim Galore aplica rotinas de filtragem por qualidade e integra relatórios do FastQC para avaliar o impacto do corte.

Entre as principais vantagens do Cutadapt e do Trim Galore estão: (i) a remoção precisa de adaptadores, prevenindo erros de montagem, (ii) a flexibilidade de configuração de parâmetros, e (iii) a integração automatizada que reduz a intervenção manual. Suas limitações incluem a necessidade de escolher cuidadosamente os limiares de qualidade, já que cortes excessivos podem resultar em perda de informação biológica relevante, enquanto cortes insuficientes podem deixar contaminantes residuais.

A escolha por essas ferramentas neste trabalho se justifica pelo seu amplo uso na comunidade científica e pela confiabilidade de seus resultados. A integração de Cutadapt e Trim Galore garante que as leituras que avançam para a montagem apresentem qualidade adequada, reduzindo o risco de artefatos e aumentando a robustez dos resultados obtidos pelo pipeline.

#### 4.2.4 NOVOPlasty

O NOVOPlasty é um montador especializado desenvolvido para reconstrução de genomas de organelas, como mitocôndrias e cloroplastos, a partir de dados de sequenciamento de genoma total. Publicado por Dierckxsens, Mardulyn e Smits (2017), tornou-se uma das ferramentas mais amplamente utilizadas para esse fim devido à sua eficiência e facilidade de uso.

Do ponto de vista técnico, o NOVOPlasty baseia-se na estratégia *seed-and-extend*. O processo de montagem é iniciado a partir de uma sequência fornecida pelo usuário — denominada *seed* — que pode ser um fragmento do genoma da organela em estudo, uma sequência de uma espécie próxima ou até mesmo uma sequência completa de organela de outro organismo. A partir dessa seed, o software realiza uma extensão bidirecional, identificando leituras sobrepostas e construindo gradualmente a sequência circular completa do genoma. Para acelerar o processo, o programa utiliza tabelas de hash, que permitem buscas rápidas por sobreposições entre as leituras.

O NOVOPlasty recebe como entrada arquivos FASTQ contendo leituras pareadas, bem como um arquivo FASTA com a sequência *seed*. A saída consiste em um arquivo FASTA representando o genoma circularizado, acompanhado de arquivos de log e, em alguns casos, múltiplas opções de montagem (quando diferentes caminhos são possíveis no grafo de montagem).

Entre as principais vantagens do NOVOPlasty estão: (i) sua alta eficiência na recuperação de genomas completos de organelas mesmo a partir de dados de genoma total, (ii) a baixa exigência de parâmetros complexos, o que o torna acessível a usuários com diferentes níveis de experiência, e (iii) a capacidade de lidar com diferentes tipos de sementes, permitindo flexibilidade em estudos de espécies com genomas pouco conhecidos. No entanto, como limitações, pode-se destacar: (i) a dependência da qualidade e da escolha adequada da *seed*, que influencia diretamente o sucesso da montagem, e (ii) o risco de circularizações incorretas em casos de regiões altamente repetitivas ou coberturas desbalanceadas.

A escolha do NOVOPlasty neste trabalho se justifica pelo fato de ser um dos montadores mais validados e aceitos na comunidade para organelas, além de sua compatibilidade com dados de leituras curtas, como os gerados por sequenciamento Illumina. Sua utilização como núcleo do pipeline permite combinar confiabilidade na montagem do mitogenoma com a inovação proposta neste TCC, que consiste em integrá-lo a uma arquitetura containerizada e orquestrada por workflows, elevando sua portabilidade e reprodutibilidade.

#### 4.2.5 Mitos2

O MITOS2 é uma ferramenta desenvolvida para a anotação funcional de genomas mitocondriais em metazoários. Apresentada por Donath et al. (2019), constitui a evolução da primeira versão do MITOS, oferecendo maior precisão na predição de genes, integração com bases de dados atualizadas e maior flexibilidade de execução.

Tecnicamente, o MITOS2 combina métodos de homologia e algoritmos de predição *ab initio* para identificar as principais classes de genes presentes em genomas mitocondriais. O sistema utiliza alinhamentos com modelos de substituição específicos para proteínas codificadas por mitocôndrias, além de perfis de RNA para identificar rRNAs e tRNAs. As entradas consistem em um arquivo FASTA contendo a sequência do genoma mitocondrial a ser anotado, e as saídas incluem: (i) arquivos tabulares com as coordenadas dos genes anotados, (ii) representações gráficas do genoma mostrando a organização gênica, e (iii) arquivos de anotação em formatos padrão, como GFF3.

Embora seja amplamente utilizado em sua versão web, o MITOS2 também pode ser executado por linha de comando em servidores Linux. Essa funcionalidade é particularmente relevante para a integração em pipelines automatizados, pois permite encapsular a ferramenta em containers Docker e orquestrar sua execução por meio de sistemas de workflow, como o Nextflow. Dessa forma, elimina-se a dependência da interface online e amplia-se a reprodutibilidade e a escalabilidade da análise.

Entre as principais vantagens do MITOS2 destacam-se: (i) a automação completa da etapa de anotação, reduzindo a necessidade de curadoria manual extensiva; (ii) a compatibilidade com diferentes grupos taxonômicos de metazoários; e (iii) a produção de relatórios gráficos de fácil interpretação. Por outro lado, como limitações, pode apresentar divergências em regiões com estruturas gênicas incomuns e requer validação complementar para garantir a anotação correta em casos atípicos.

A escolha do MITOS2 neste trabalho se justifica por sua ampla aceitação na comunidade científica e por sua capacidade de integrar de forma prática e padronizada a anotação de genomas mitocondriais reconstruídos. Sua utilização no pipeline proposto assegura que os genomas montados pelo NOVOPlasty possam ser caracterizados de forma robusta e reprodutível, viabilizando estudos comparativos de organização gênica e aplicações filogenéticas.

#### 4.2.6 tRNAscan-SE

O tRNAscan-SE é a principal ferramenta disponível para a identificação de genes de RNA transportador (tRNA) em genomas de procariotos e eucariotos. Originalmente desenvolvido por Lowe e Eddy nos anos 1990, foi posteriormente atualizado e consolidado em sua versão atual, o tRNAscan-SE 2.0, que apresenta melhorias significativas em termos de sensibilidade, especificidade e velocidade (LOWE; CHAN, 2016).

Do ponto de vista técnico, o tRNAscan-SE utiliza uma combinação de métodos heurísticos rápidos e algoritmos baseados em *covariance models*, que permitem detectar sequências de tRNA com alta precisão, mesmo em genomas complexos. O software recebe como entrada arquivos FASTA contendo sequências genômicas e gera como saída relatórios tabulares com as coordenadas dos genes de tRNA previstos, bem como informações sobre suas estruturas secundárias. Além disso, é capaz de distinguir entre tRNAs funcionais e pseudogenes de tRNA, uma funcionalidade importante para estudos de genomas mitocondriais, que frequentemente apresentam regiões com estruturas atípicas.

Entre as vantagens do tRNAscan-SE destacam-se: (i) a alta sensibilidade para detectar tRNAs genuínos, (ii) a ampla aceitação como padrão-ouro para anotação de tRNAs, e (iii) a integração com outros softwares de anotação, como o MITOS2. Por outro lado, como limitações, pode apresentar maior tempo de execução em genomas muito grandes e depender de ajustes de parâmetros quando aplicado a organismos com estruturas de tRNA altamente divergentes.

A escolha do tRNAscan-SE neste trabalho se justifica por sua confiabilidade e ampla utilização na comunidade científica. Sua inclusão no pipeline assegura que a anotação dos genomas mitocondriais reconstruídos seja completa, contemplando não apenas os genes codificadores de proteínas e rRNAs, mas também os tRNAs, que são essenciais para a funcionalidade do genoma e frequentemente utilizados em análises comparativas e filogenéticas. Dessa forma, o tRNAscan-SE atua como um componente complementar e indispensável para a anotação de mitogenomas.

#### 4.2.7 Docker

No pipeline proposto, cada ferramenta foi encapsulada em um container Docker, assegurando que todas as dependências fossem executadas em ambientes controlados e consistentes (MERKEL, 2014). Essa abordagem eliminou a necessidade de configurações manuais complexas e garantiu que softwares como FastQC, Cutadapt, NOVOPlasty e MITOS2 operassem de forma uniforme em diferentes sistemas operacionais. Com isso, a montagem e a anotação de mitogenomas tornam-se replicáveis em qualquer ambiente computacional que disponha do Docker instalado, ampliando a portabilidade do pipeline.

#### 4.2.8 Nextflow

A orquestração do pipeline foi implementada no Nextflow, que estruturou os processos em etapas independentes, conectando automaticamente suas entradas e saídas (DI TOMMASO et al., 2017). Essa organização permitiu automatizar a lógica iterativa da montagem, como a possibilidade de reaproveitar contigs parciais do NOVOPlasty como novas seeds em execuções subsequentes. Além disso, o Nextflow possibilitou a paralelização de tarefas, a retomada de execuções interrompidas e a integração direta com containers Docker. Dessa forma, o pipeline desenvolvido combina simplicidade de uso com escalabilidade e reprodutibilidade, adequando-se tanto a computadores locais quanto a servidores de alto desempenho.

#### 4.2.9 GitHub

O GitHub é uma plataforma de hospedagem de código baseada no sistema de controle de versão Git, amplamente utilizada para o desenvolvimento colaborativo de software e projetos científicos. Desde sua criação em 2008, tornou-se um dos principais repositórios públicos de código-fonte, reunindo milhões de projetos de diferentes áreas. Na ciência, seu uso tem crescido de forma expressiva, acompanhando o movimento de ciência aberta e reprodutibilidade computacional.

Do ponto de vista técnico, o GitHub oferece funcionalidades essenciais para o gerenciamento de projetos: (i) controle de versão, permitindo rastrear todas as modificações realizadas em arquivos de código ou documentação; (ii) ramificações (branches), que possibilitam o desenvolvimento paralelo de diferentes versões ou funcionalidades; e (iii) integração contínua, que permite automatizar testes e validações sempre que alterações são feitas. Além disso, a plataforma possui suporte para issues e pull requests, que facilitam o gerenciamento de contribuições externas e a discussão entre colaboradores.

Para pipelines de bioinformática, a disponibilização em um repositório GitHub garante não apenas transparência, mas também acessibilidade e manutenção a longo prazo. Repositórios podem incluir: (i) o código do workflow, (ii) arquivos de configuração e exemplos de execução, (iii) Dockerfiles das ferramentas utilizadas, e (iv) documentação detalhada, como tutoriais e guias de instalação.

As principais vantagens do GitHub incluem: (i) a difusão do projeto, já que o repositório pode ser acessado por qualquer pesquisador; (ii) a rastreabilidade, permitindo verificar quando e como cada alteração foi realizada; e (iii) a integração com iniciativas comunitárias como o nf-core, que adota o GitHub como plataforma padrão para o compartilhamento e padronização de pipelines bioinformáticos (EWELS et al., 2020). Como limitação, pode-se citar a necessidade de acesso à internet para sincronização dos repositórios, o que pode restringir seu uso em ambientes de HPC com restrições de rede.

Neste trabalho, o GitHub é utilizado como repositório público para disponibilizar o pipeline desenvolvido, assegurando sua acessibilidade, versionamento e documentação adequada. Essa prática se alinha às diretrizes de ciência aberta e reforça o compromisso com a reprodutibilidade e a transparência metodológica.

### 4.3 Aquisição e Preparo dos Dados

A primeira etapa prática do pipeline consiste na aquisição e organização dos dados de sequenciamento que servirão de insumo para a montagem dos genomas mitocondriais. Para este trabalho, foram utilizados datasets públicos provenientes do Sequence Read Archive (SRA), repositório mantido pelo NCBI e reconhecido como um dos principais bancos internacionais de dados de sequenciamento de alto rendimento (LEINONEN; SUGAWARA; SHUMWAY, 2011). A escolha por dados públicos justifica-se pela sua ampla disponibilidade, pela diversidade de organismos representados e pela possibilidade de reproduzir e validar resultados obtidos em trabalhos prévios.

O acesso aos dados é realizado por meio do SRA-Toolkit. O processo inicia-se com o comando `prefetch`, responsável pelo download dos arquivos brutos no formato `.sra`, seguido da conversão para o formato FASTQ utilizando o `fasterq-dump`.

Dois datasets foram selecionados para este trabalho, com propósitos complementares:

**a) Espécie de validação — *Diploprion bifasciatus* (peixe-sabão barrado):** O dataset SRR36182901, gerado por sequenciamento *paired-end* na plataforma Illumina NovaSeq X Plus (2×151 bp), foi selecionado como caso de teste por três razões: (i) seu mitogenoma (PZ143763.1, 16.805 bp) foi publicado a partir deste mesmo dataset, permitindo comparação direta; (ii) o volume de dados é compacto (~1,2 GB comprimido, 12 milhões de pares de leituras), viabilizando testes rápidos e iterativos durante o desenvolvimento do pipeline; e (iii) a espécie pertence a um grupo taxonômico distinto (Perciformes), demonstrando a generalidade do pipeline para além de aves.

**b) Espécie principal — *Anodorhynchus leari* (arara-azul-de-lear):** O dataset SRR28399504, gerado por sequenciamento *paired-end* na plataforma Illumina HiSeq X Ten (2×150 bp), corresponde ao caso de estudo central deste TCC. Depositado por Iridian Genomes em 2024, o dataset contém 118,5 milhões de pares de leituras (~10,7 GB comprimido / ~87 GB em FASTQ). Por se tratar de dados de genoma total (WGS), a fração mitocondrial corresponde a uma pequena proporção do total — porém, dada a abundância natural de cópias de mtDNA por célula, essa fração é suficiente para a montagem completa do mitogenoma.

**Tabela 1** – Datasets utilizados no pipeline.

| Espécie | Accession | Plataforma | Reads totais | Tamanho SRA | Finalidade |
|---|---|---|---|---|---|
| *D. bifasciatus* | SRR36182901 | NovaSeq X Plus (2×151 bp) | 12 M pares | ~1,2 GB | Validação |
| *A. leari* | SRR28399504 | HiSeq X Ten (2×150 bp) | 118,5 M pares | ~11 GB | Estudo principal |

Para ambas as espécies, sequências-semente (*seeds*) do gene COX1 foram obtidas de espécies filogeneticamente próximas ou da própria espécie a partir de mitogenomas depositados no GenBank. Para *A. leari*, utilizou-se o COX1 de *A. hyacinthinus* (NC_082165.1, posições 5359–6906, 1.548 bp), espécie-irmã do mesmo gênero. Para *D. bifasciatus*, utilizou-se o COX1 extraído do próprio mitogenoma publicado (PZ143763.1, 1.560 bp). As sementes foram armazenadas em formato FASTA no diretório `data/seeds/` do repositório.

**Estratégia de truncagem** — Devido ao grande volume do dataset de *A. leari* (118,5 M reads), a conversão completa para FASTQ gera aproximadamente 87 GB de dados. Para reduzir o custo computacional nas etapas subsequentes sem comprometer a montagem do mitogenoma, foi implementada uma estratégia de truncagem no módulo `sra_download.nf`: após a conversão completa pelo `fasterq-dump`, o comando `head -n (N × 4)` é utilizado para reter apenas os primeiros 20 milhões de pares de leituras (sendo 4 linhas por read no formato FASTQ). Essa abordagem reduziu o volume de cada arquivo de ~43,5 GB para ~7,4 GB. A escolha de 20 milhões de reads foi baseada na estimativa de cobertura:

> C = (N × f × 2 × L) / G

Onde N = 20.000.000 reads, f = fração mitocondrial estimada (~0,17%), L = 150 bp e G = 17.000 bp, resultando em uma cobertura esperada de aproximadamente 600×, valor amplamente suficiente para montagem *de novo* de mitogenomas.

A truncagem é realizada antes da trimagem (etapa seguinte), deliberadamente: como a qualidade das bases e a proporção de adaptadores são uniformemente distribuídas ao longo de uma corrida Illumina HiSeq X Ten, trimar primeiro os 118,5 M reads para depois truncar seria computacionalmente ~6× mais custoso sem ganho biológico. A ordem truncagem → trimagem é, portanto, a mais eficiente para este cenário.

Após a aquisição e truncagem, os dados são organizados em diretórios padronizados pelo próprio Nextflow (diretório `work/`), respeitando uma estrutura hierárquica que diferencia arquivos brutos, processados e resultados. Ao adotar essa prática, reforça-se o alinhamento do projeto com os princípios FAIR (WILKINSON et al., 2016).

### 4.4 Controle de Qualidade e Pré-processamento

Após a aquisição e organização dos dados de sequenciamento, a etapa seguinte do pipeline consiste na avaliação da qualidade das leituras e na aplicação de procedimentos de filtragem. Essa fase é essencial para assegurar que apenas dados confiáveis avancem para a montagem do genoma mitocondrial, evitando que artefatos técnicos comprometam a acurácia do resultado final.

O primeiro passo é a execução do FastQC (ANDREWS, 2010), que gera relatórios interativos em HTML contendo métricas fundamentais: (i) qualidade por posição nas sequências, (ii) distribuição de valores de qualidade ao longo de todas as leituras, (iii) conteúdo GC, (iv) presença de adaptadores residuais e (v) detecção de sequências sobre-representadas. No caso do dataset de *A. leari* (SRR28399504), o FastQC confirmou a presença significativa de sequências de adaptadores Illumina nas reads brutas, justificando a necessidade da etapa de trimagem.

Em seguida, foi aplicado o Trim Galore (BABRAHAM BIOINFORMATICS, 2019), wrapper que integra o Cutadapt (MARTIN, 2011) com o FastQC, automatizando a detecção e remoção de adaptadores e bases de baixa qualidade. O Trim Galore opera em modo *paired-end*, mantendo a sincronização entre R1 e R2 — se uma read é descartada por comprimento mínimo, a parceira correspondente também é removida.

Para o dataset de *A. leari* (20 milhões de reads truncados), o Trim Galore detectou automaticamente o adaptador Illumina TruSeq (`AGATCGGAAGAGC`) e produziu os seguintes resultados:

**Tabela 2** – Resultados do pré-processamento com Trim Galore para *A. leari* (SRR28399504).

| Métrica | R1 (forward) | R2 (reverse) |
|---|---|---|
| Reads de entrada | 20.000.000 | 20.000.000 |
| Reads de saída | 20.000.000 (100%) | 20.000.000 (100%) |
| Bases de entrada | 3,00 Gb | 3,00 Gb |
| Bases de saída | 2,15 Gb (71,8%) | 2,16 Gb (72,0%) |
| Removido por baixa qualidade | 30 Mb (1,0%) | 20 Mb (0,7%) |
| Removido por adaptadores | 816 Mb (27,2%) | 821 Mb (27,4%) |
| Reads com adaptador detectado | 16.285.736 (81,4%) | 16.157.912 (80,8%) |
| Adaptador identificado | AGATCGGAAGAGC (TruSeq) | AGATCGGAAGAGC |

A alta proporção de reads com adaptador (81%) é um valor esperado em bibliotecas de fragmentos curtos sequenciadas em plataformas que geram reads de 150 bp: quando o fragmento de DNA inserido é menor que o comprimento da read, o sequenciador avança até a região do adaptador, que é então capturado na sequência. A remoção dessas sequências artificiais é indispensável para evitar erros de montagem.

Apenas 1% das bases foi removido por baixa qualidade (Phred < 20), indicando que a corrida de sequenciamento foi de alta qualidade. Notavelmente, 100% das reads passaram pelos filtros de comprimento mínimo, sem perda de pares — evidência de que o dataset possui consistência e profundidade adequadas.

Além de aumentar a confiabilidade biológica das leituras, o pré-processamento contribuiu significativamente para a redução do volume de dados:

**Tabela 3** – Redução progressiva de volume de dados para *A. leari*.

| Etapa | Tamanho por arquivo | Reads (pares) | Redução |
|---|---|---|---|
| SRA comprimido (.sra) | 11,0 GB (total) | 118,5 M | — |
| FASTQ bruto (fasterq-dump) | ~43,5 GB | 118,5 M | +295% (descompressão) |
| FASTQ truncado (head 20M) | ~7,4 GB | 20 M | −83% |
| FASTQ trimado (Trim Galore) | ~5,3 GB | ~20 M | −28% |

Essa redução progressiva demonstra como cada etapa do pipeline contribui para tornar o processamento mais eficiente, sem comprometer a representatividade dos dados para a montagem do mitogenoma.

### 4.5 Montagem do Mitogenoma

A montagem do genoma mitocondrial constitui a etapa central do pipeline, pois é nela que se obtém a sequência circular completa que servirá de base para a anotação funcional e análises comparativas. Para essa finalidade foi utilizado o NOVOPlasty, um dos montadores mais consolidados para organelas, reconhecido por sua eficiência em reconstruir mitogenomas e cloroplastos a partir de dados de genoma total (DIERCKXSENS; MARDULYN; SMITS, 2017).

O NOVOPlasty adota a estratégia de *seed-and-extend*, em que a montagem é iniciada a partir de uma sequência semente (*seed*) fornecida pelo usuário. Essa seed pode ser um fragmento de gene mitocondrial (como *cox1*) proveniente da própria espécie em estudo ou de uma espécie próxima. A partir desse ponto inicial, o programa busca leituras sobrepostas no conjunto de dados e realiza a extensão bidirecional até completar a reconstrução do genoma. O processo resulta em uma molécula circular característica do DNA mitocondrial, armazenada em um arquivo FASTA, acompanhada de arquivos de log e relatórios de execução.

Do ponto de vista técnico, o NOVOPlasty recebe como entrada os arquivos FASTQ pareados, resultantes do pré-processamento, além da sequência seed e de um arquivo de configuração com parâmetros como tamanho esperado do genoma e valor de k-mer. Como saída, além do genoma circularizado, o software pode gerar múltiplas versões da montagem quando diferentes caminhos são possíveis no grafo de sobreposição das leituras.

Cabe destacar que a montagem de mitogenomas pode ser um processo iterativo. Nem sempre a primeira execução do NOVOPlasty resulta em um genoma circular completo; nesses casos, o contig mais longo obtido ou o resultado parcial da montagem pode ser utilizado como nova seed em uma segunda execução. Essa estratégia aumenta a chance de se alcançar a circularização correta e reflete uma prática consolidada em bioinformática de organelas.

Uma vez que a montagem pode exigir múltiplas execuções até alcançar a circularização, o pipeline foi projetado para contemplar essa necessidade de iteração. Quando a primeira execução do NOVOPlasty não resultar em um genoma completo, o próprio fluxo permite reutilizar automaticamente o contig mais longo obtido como nova seed em uma segunda rodada. Dessa forma, a lógica iterativa é incorporada ao workflow de maneira transparente, sem exigir que o usuário reinicie manualmente o processo.

Neste trabalho, a montagem foi conduzida para dois casos de estudo: (i) *Diploprion bifasciatus* como espécie de validação, e (ii) *Anodorhynchus leari* como espécie principal. Para *A. leari*, utilizou-se como semente o gene COX1 de *A. hyacinthinus* (NC_082165.1), espécie-irmã do mesmo gênero, com genoma mitocondrial de referência de 16.999 bp.

### 4.6 Pós-processamento e Anotação

Concluída a etapa de montagem, torna-se essencial realizar o pós-processamento do genoma obtido, tanto para verificar sua integridade quanto para proceder à anotação funcional. Essa fase garante que o produto gerado pelo NOVOPlasty represente de fato um genoma mitocondrial completo, circularizado e biologicamente plausível, apto a servir como base para análises comparativas, filogenéticas e de conservação.

A avaliação inicial da montagem é conduzida por meio da inspeção de métricas como o tamanho total do genoma, que em aves situa-se geralmente entre 16 e 20 kb, a verificação da circularidade da molécula e a presença dos genes esperados em genomas mitocondriais típicos, incluindo 13 genes codificadores de proteínas, 2 rRNAs e 22 tRNAs (BOORE, 1999).

A etapa seguinte consiste na anotação funcional do genoma reconstruído, para a qual será utilizado o MITOS2, plataforma amplamente validada para predição automática de genes mitocondriais em metazoários (DONATH et al., 2019). De maneira complementar, será empregado o tRNAscan-SE, considerado padrão-ouro para a detecção de genes de tRNA (LOWE; CHAN, 2016). A combinação de ambas as ferramentas confere maior robustez ao processo de anotação.

Os produtos finais dessa etapa incluirão: (i) um arquivo FASTA anotado contendo as sequências gênicas identificadas; (ii) um arquivo GFF3 com as coordenadas de cada gene; (iii) relatórios tabulares e representações gráficas do genoma circular; e (iv) arquivos de estruturas secundárias dos RNAs previstos.

### 4.7 Automação e Reprodutibilidade

A automação do pipeline e a padronização de sua execução representam os principais diferenciais desta proposta em relação a abordagens convencionais de montagem de genomas mitocondriais. Mais do que reunir ferramentas em sequência, busca-se disponibilizar um fluxo de trabalho que possa ser aplicado por diferentes grupos de pesquisa em distintos ambientes computacionais, assegurando consistência analítica, transparência metodológica e facilidade de reutilização.

Para atingir esse objetivo, cada etapa do pipeline foi encapsulada em containers Docker, garantindo que as ferramentas operem em ambientes controlados e reproduzíveis, independentemente do sistema em que forem executadas (MERKEL, 2014).

A orquestração é realizada pelo Nextflow, que automatiza a execução dos processos, gerencia dependências e possibilita paralelização quando aplicável. O sistema ainda permite a retomada da análise a partir do ponto de falha, evitando reprocessamento desnecessário (DI TOMMASO et al., 2017).

Além disso, todo o código-fonte do workflow é disponibilizado em repositório público no GitHub (https://github.com/matheus-sobreira/mitogenome-pipeline), acompanhado de documentação detalhada e exemplos de uso. A execução é simplificada por meio de perfis de configuração por espécie:

```
# Executar para A. leari (estudo principal)
nextflow run main.nf -profile a_leari

# Executar para D. bifasciatus (validação)
nextflow run main.nf -profile test
```

Essa prática de ciência aberta não apenas reforça a transparência metodológica, mas também incentiva a adaptação e o aprimoramento do pipeline por outros pesquisadores, em consonância com iniciativas comunitárias como o nf-core (EWELS et al., 2020).

---

## 5. Resultados e Discussão

Este capítulo apresenta os resultados obtidos na execução do pipeline de montagem de mitogenomas, aplicado a dois estudos de caso: (i) *Diploprion bifasciatus* como espécie de validação, cujo mitogenoma publicado permite comparação direta, e (ii) *Anodorhynchus leari* como espécie principal, objeto central deste TCC. Os resultados são descritos por etapa do pipeline, seguidos de uma discussão sobre a robustez da abordagem e sua relevância biológica e científica.

### 5.1 Resultados por Etapa do Pipeline

#### 5.1.1 Aquisição e Preparo dos Dados

A etapa de aquisição de dados foi executada com sucesso para ambos os datasets. O SRA-Toolkit (`prefetch` + `fasterq-dump`) realizou o download e a conversão dos arquivos `.sra` para FASTQ sem erros. Para o dataset de *A. leari* (SRR28399504), o `prefetch` obteve o arquivo `.sra` de 11 GB, e o `fasterq-dump` realizou a conversão completa em dois arquivos FASTQ totalizando ~87 GB.

Uma correção técnica relevante foi necessária durante o desenvolvimento: a flag `-X`, utilizada no antigo `fastq-dump` para limitar o número de reads, não é suportada pelo `fasterq-dump` (versão 3.0.10), gerando o erro *"Unknown argument '-X'"* com código de saída 3. A solução implementada no módulo `sra_download.nf` foi a conversão completa seguida de truncagem via `head -n (N × 4)`, onde N é o número de reads desejado. Essa abordagem, embora exija armazenamento temporário dos arquivos completos, é robusta e independente de versão da ferramenta.

Para *A. leari*, a truncagem de 118,5 M para 20 M reads reduziu cada arquivo de ~43,5 GB para ~7,4 GB, conforme planejado.

#### 5.1.2 Controle de Qualidade

O FastQC (versão 0.12.1) gerou relatórios de qualidade para ambos os datasets. Para *A. leari*, os relatórios confirmaram perfis de qualidade típicos de dados Illumina HiSeq X Ten, com *scores* Phred médios acima de Q30 em todas as posições, presença significativa de adaptadores Illumina TruSeq e níveis de duplicação moderados — padrão esperado em bibliotecas WGS de alta cobertura.

#### 5.1.3 Pré-processamento

O Trim Galore (versão 0.6.10, com Cutadapt 4.6) foi executado em modo *paired-end*, com detecção automática de adaptadores. Os resultados detalhados para *A. leari* estão apresentados na Tabela 2 (Seção 4.4). Em resumo: 81,4% das reads R1 e 80,8% das R2 continham sequências de adaptador Illumina TruSeq (`AGATCGGAAGAGC`), resultando na remoção de 27,2–27,4% das bases por adaptadores e 0,7–1,0% por baixa qualidade. Todas as 20 milhões de reads passaram nos filtros de comprimento mínimo (100% de retenção).

Esses valores indicam que a biblioteca de sequenciamento continha fragmentos de inserção relativamente curtos, resultando em *read-through* dos adaptadores — comportamento normal e esperado em bibliotecas WGS padrão. O resultado confirma a importância da etapa de trimagem como passo obrigatório antes da montagem.

#### 5.1.4 Montagem do Mitogenoma

O NOVOPlasty (versão 4.3.1) foi executado com sucesso para ambas as espécies, produzindo montagens circularizadas na primeira execução — sem necessidade da estratégia iterativa de re-seeding com contigs parciais.

**Tabela 4** – Resultados da montagem por NOVOPlasty.

| Métrica | *D. bifasciatus* (validação) | *A. leari* (principal) |
|---|---|---|
| Status da montagem | Circularizado ✓ | Circularizado ✓ |
| Número de contigs | 1 | 1 |
| Tamanho do contig | — | 16.986 bp |
| Referência mais próxima | PZ143763.1 (16.805 bp) | *A. hyacinthinus* (16.999 bp) |
| K-mer utilizado | 33 | 39 |
| Reads totais processadas | — | 20.285.220 |
| Reads alinhadas ao mtDNA | — | 34.644 (0,17%) |
| Cobertura média | — | 306× |
| Fração subsampleada pelo NOVOPlasty | — | 61,9% |
| Insert size observado | — | 213 bp |
| Semente utilizada | COX1 *D. bifasciatus* (PZ143763.1) | COX1 *A. hyacinthinus* (NC_082165.1) |

Para *A. leari*, o mitogenoma de **16.986 bp** difere em apenas 13 bp (< 0,1%) do mitogenoma de referência de *A. hyacinthinus* (16.999 bp), resultado biologicamente coerente para espécies do mesmo gênero. A semente COX1 interespecífica (*A. hyacinthinus* → *A. leari*) foi recuperada com sucesso, demonstrando que o NOVOPlasty tolera bem diferenças genéticas dentro do mesmo gênero.

A fração de reads mitocondriais (0,17% do total) produziu cobertura média de 306×, amplamente acima do limiar mínimo de 50× necessário para montagem confiável de organelas. Isso confirma a adequação da estratégia de truncagem para 20 M reads: mesmo utilizando apenas 16,9% do dataset original, foi possível obter cobertura elevada e montagem circular completa.

#### 5.1.5 Validação com *D. bifasciatus*

A espécie de validação (*D. bifasciatus*, SRR36182901) forneceu uma montagem circularizada com k-mer 33 na primeira execução. A comparação com o mitogenoma publicado (PZ143763.1, gerado a partir do mesmo dataset) permite validar a corretude do pipeline e confirmar que as etapas de download, trimagem e montagem preservam a integridade biológica dos dados. A montagem foi bem-sucedida, confirmando a funcionalidade do pipeline em organismos de diferentes grupos taxonômicos (peixes e aves).

### 5.2 Discussão sobre a Robustez da Abordagem

A robustez metodológica do pipeline foi demonstrada pela obtenção de montagens completas e circularizadas em ambos os estudos de caso, abrangendo organismos de grupos taxonômicos distintos (Perciformes e Psittaciformes) e plataformas de sequenciamento diferentes (NovaSeq X Plus e HiSeq X Ten).

Entre os aspectos positivos observados, destaca-se a **reprodutibilidade**, garantida pela combinação entre a containerização das ferramentas via Docker e a orquestração do fluxo pelo Nextflow. Cada execução do pipeline produz um diretório `work/` com hash único por processo, permitindo rastreabilidade completa de cada etapa. O uso do flag `-resume` permite retomar execuções interrompidas sem reprocessamento, característica que se mostrou essencial durante o desenvolvimento — quando falhas de ferramentas (como a incompatibilidade da flag `-X` no fasterq-dump) puderam ser corrigidas sem necessidade de reiniciar todo o pipeline.

A **estratégia de truncagem** por `head`, embora exija o armazenamento temporário dos arquivos FASTQ completos, mostrou-se robusta e produziu resultados de montagem indistinguíveis dos que seriam obtidos com o dataset completo. A justificativa reside na uniformidade de qualidade e composição ao longo das corridas Illumina, que torna a truncagem posicional estatisticamente equivalente a uma amostragem aleatória para fins de montagem mitogenômica.

O módulo de NOVOPlasty implementado no pipeline inclui uma **lógica iterativa** com múltiplos k-mers e re-seeding automático. Embora essa funcionalidade não tenha sido necessária nos dois casos testados (ambos circularizaram na primeira tentativa), ela confere resiliência ao pipeline para cenários mais desafiadores, como baixa cobertura ou regiões repetitivas complexas.

Como **limitação** observada na prática, o volume de armazenamento temporário durante a etapa de download/conversão pode ser significativo: para o dataset de *A. leari*, foram necessários ~11 GB (arquivo `.sra`) + ~87 GB (FASTQ bruto temporário) + ~15 GB (FASTQ truncado), totalizando ~113 GB de pico de armazenamento antes da limpeza automática. Esse requisito pode ser restritivo em ambientes com espaço limitado, embora o pipeline realize limpeza automática das reads brutas após cada etapa.

Em comparação com trabalhos relacionados, como o MitoHiFi (ULIANO-SILVA et al., 2023), que explora o potencial das leituras longas de alta fidelidade, o presente pipeline se diferencia por ser desenhado especificamente para dados de leituras curtas, ainda predominantes em muitas plataformas de sequenciamento. Essa escolha amplia sua aplicabilidade, dado que dados Illumina continuam sendo os mais acessíveis e economicamente viáveis para a maioria dos grupos de pesquisa.

### 5.3 Relevância Biológica e Científica

A relevância deste trabalho transcende a dimensão metodológica. Do ponto de vista biológico, o pipeline possibilitou a reconstrução reprodutível do mitogenoma de *Anodorhynchus leari*, espécie endêmica do Brasil classificada como "Em Perigo" de extinção. O mitogenoma obtido — **16.986 bp, circularizado, com 306× de cobertura** — constitui um recurso genômico validável e reprodutível, diferenciando-se da montagem anterior (realizada em contexto acadêmico sem documentação metodológica padronizada) por estar associado a:

(i) um pipeline público e versionado no GitHub;
(ii) containers Docker com versões fixas de cada ferramenta;
(iii) arquivos de configuração declarativos que permitem reprodução exata do experimento;
(iv) relatórios de qualidade e trimagem preservados como saída do pipeline.

A similaridade de 99,9% entre o mitogenoma montado (16.986 bp) e o de *A. hyacinthinus* (16.999 bp) é biologicamente coerente, dado que as duas espécies são irmãs dentro do gênero *Anodorhynchus*. Essa proximidade valida tanto a escolha da semente interespecífica quanto a qualidade da montagem obtida.

A validação cruzada com *D. bifasciatus* (grupo taxonômico distinto) demonstra que o pipeline não é específico para aves ou para um tipo particular de dados, ampliando sua aplicabilidade para a comunidade científica. A disponibilização pública do pipeline, acompanhada de perfis de execução para cada espécie (`nextflow run main.nf -profile a_leari` ou `-profile test`), permite que qualquer pesquisador reproduza integralmente os resultados aqui apresentados.

Em termos mais amplos, a disponibilização de um pipeline automatizado, containerizado e orquestrado com Nextflow se alinha às demandas contemporâneas da bioinformática por soluções reprodutíveis, portáveis e escaláveis. Ao tornar esse recurso acessível por meio de repositório público, acompanhado de documentação detalhada, o presente trabalho contribui para o fortalecimento de práticas de ciência aberta e colabora com a comunidade internacional que atua no desenvolvimento de workflows bioinformáticos.

---

## 6. Considerações Finais

O presente trabalho teve como objetivo desenvolver um pipeline de bioinformática automatizado, reprodutível e portável para a montagem *de novo* de genomas mitocondriais. A integração de ferramentas consolidadas — SRA-Toolkit, FastQC, Trim Galore/Cutadapt e NOVOPlasty — em um fluxo containerizado com Docker e orquestrado pelo Nextflow constituiu a principal contribuição metodológica desta proposta.

O pipeline foi validado com sucesso em dois cenários distintos. A espécie de validação, *Diploprion bifasciatus*, forneceu uma montagem circularizada comparável ao mitogenoma publicado para a mesma espécie. O estudo de caso principal, *Anodorhynchus leari*, resultou em um mitogenoma de **16.986 bp**, circularizado com cobertura média de **306×**, com diferença inferior a 0,1% em relação à espécie-irmã *A. hyacinthinus*, confirmando a validade biológica da montagem.

Os resultados demonstraram que:

a) O pipeline é capaz de reconstruir mitogenomas completos e circularizados a partir de dados públicos de WGS, utilizando sementes interespecíficas quando a sequência da espécie-alvo não está disponível;

b) A estratégia de truncagem para 20 milhões de reads (de um total de 118,5 milhões) foi suficiente para produzir cobertura de 306×, validando a abordagem de redução de volume como estratégia viável para economizar recursos computacionais;

c) A containerização individual das ferramentas eliminou problemas de incompatibilidade de versões e dependências, materializando o princípio de reprodutibilidade no próprio design do pipeline;

d) A modularidade do Nextflow DSL2, aliada a perfis de configuração por espécie (`-profile a_leari`, `-profile test`), permite que o pipeline seja reutilizado para qualquer organismo com ajustes mínimos — bastando fornecer um accession SRA e uma semente COX1.

Como perspectivas futuras, destacam-se: (i) a implementação da Fase 2 do pipeline — anotação funcional com MITOS2 e tRNAscan-SE, que identificará os 37 genes típicos do mitogenoma; (ii) a comparação detalhada do mitogenoma de *A. leari* com outras espécies de Psittaciformes para estudos filogenéticos; (iii) a publicação do mitogenoma no GenBank; e (iv) a adequação completa do documento ao modelo ABNT/UERN no Overleaf (LaTeX).

Em síntese, o pipeline desenvolvido representa uma contribuição metodológica concreta para a bioinformática de organelas, demonstrando que é possível unir rigor técnico, aplicabilidade biológica e boas práticas de reprodutibilidade computacional em uma solução acessível e reutilizável. A aplicação ao caso de uma espécie ameaçada da biodiversidade brasileira confere ao trabalho não apenas relevância acadêmica, mas também potencial impacto para a conservação.

---

## 7. Referências

ANDERSON, S.; et al. Sequence and organization of the human mitochondrial genome. *Nature*, v. 290, p. 457-465, 1981.

ANDREWS, S. FastQC: a quality control tool for high throughput sequence data. 2010. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Acesso em: 31 ago. 2025.

ASAKAWA, S.; HIMENO, H.; MIURA, K.; WATANABE, K. Nucleotide sequence and gene organization of the starfish *Asterina pectinifera* mitochondrial genome. *Genetics*, v. 140, n. 3, p. 1047-1060, 1995.

BABRAHAM BIOINFORMATICS. Trim Galore! 2019. Disponível em: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. Acesso em: 31 ago. 2025.

BAKER, M. 1,500 scientists lift the lid on reproducibility. *Nature*, v. 533, n. 7604, p. 452-454, 2016.

BOORE, J. L. Animal mitochondrial genomes. *Nucleic Acids Research*, v. 27, n. 8, p. 1767-1780, 1999.

BRITISH ECOLOGICAL SOCIETY. *A guide to reproducible code in ecology and evolution*. London: British Ecological Society, 2017.

CALABRESE, C.; et al. MToolBox: a bioinformatics pipeline for the characterization of human mitochondrial genome variants. *Bioinformatics*, v. 30, n. 21, p. 3115-3117, 2014. doi: 10.1093/bioinformatics/btu483.

DIERCKXSENS, N.; MARDULYN, P.; SMITS, G. NOVOPlasty: de novo assembly of organelle genomes from whole genome data. *Nucleic Acids Research*, v. 45, n. 4, p. e18, 2017.

DI TOMMASO, P.; et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology*, v. 35, n. 4, p. 316-319, 2017.

DONATH, A.; et al. Improved annotation of protein-coding gene boundaries in metazoan mitochondrial genomes. *Molecular Ecology Resources*, v. 19, n. 4, p. 609-615, 2019. doi: 10.1111/1755-0998.12985.

EKBLOM, R.; WOLF, J. B. W. A field guide to whole-genome sequencing, assembly and annotation. *Evolutionary Applications*, v. 7, n. 9, p. 1026-1042, 2014.

EWELS, P. A.; et al. The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*, v. 38, n. 3, p. 276-278, 2020.

GRÜNING, B.; et al. Practical computational reproducibility in the life sciences. *Cell Systems*, v. 6, n. 6, p. 631-635, 2018.

HAHN, C.; BACHMANN, L.; CHEVREUX, B. Reconstructing mitochondrial genomes directly from genomic next-generation sequencing reads — a baiting and iterative mapping approach. *Nucleic Acids Research*, v. 41, n. 13, p. e129, 2013. doi: 10.1093/nar/gkt371.

JIN, J.-J.; et al. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology*, v. 21, n. 1, p. 241, 2020. doi: 10.1186/s13059-020-02154-5.

KÖSTER, J.; RAHMANN, S. Snakemake—a scalable bioinformatics workflow engine. *Bioinformatics*, v. 28, n. 19, p. 2520-2522, 2012.

LEINONEN, R.; SUGAWARA, H.; SHUMWAY, M. The Sequence Read Archive. *Nucleic Acids Research*, v. 39, supl. 1, p. D19-D21, 2011. doi: 10.1093/nar/gkq1019.

LOWE, T. M.; CHAN, P. P. tRNAscan-SE On-line: integrating search and context for analysis of transfer RNA genes. *Nucleic Acids Research*, v. 44, supl. W1, p. W54-W57, 2016. doi: 10.1093/nar/gkw413.

MARTIN, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, v. 17, n. 1, p. 10-12, 2011. doi: 10.14806/ej.17.1.200.

MERKEL, D. Docker: lightweight Linux containers for consistent development and deployment. *Linux Journal*, n. 239, p. 2, 2014.

MILLER, J. R.; KOREN, S.; SUTTON, G. Assembly algorithms for next-generation sequencing data. *Genomics*, v. 95, n. 6, p. 315-327, 2010.

ORELLANA, L. H.; et al. Comparing genomes recovered from time-series metagenomes using long- and short-read sequencing technologies. *Microbiome*, v. 11, n. 1, p. 105, 2023.

PENG, R. D. Reproducible research in computational science. *Science*, v. 334, n. 6060, p. 1226-1227, 2011. DOI: https://doi.org/10.1126/science.1213847.

SANDVE, G. K.; et al. Ten simple rules for reproducible computational research. *PLoS Computational Biology*, v. 9, n. 10, p. e1003285, 2013.

ULIANO-SILVA, M.; et al. MitoHiFi: a pipeline for mitochondrial genome assembly from PacBio HiFi reads. *BMC Bioinformatics*, v. 24, n. 288, 2023. DOI: 10.1186/s12859-023-05385-y.

WILKINSON, M. D.; et al. The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data*, v. 3, p. 160018, 2016.
