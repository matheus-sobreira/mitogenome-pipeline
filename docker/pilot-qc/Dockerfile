# ──────────────────────────────────────────────────────────────────────────────
# Pilot QC — Análise da amostra piloto (qualidade + fração mitocondrial)
# Versões: BWA 0.7.18 | Samtools 1.19.2 | Seqtk 1.4
# Base: ubuntu:22.04
#
# Separada da imagem sra-tools de propósito (ver DEC-03 no vault): a etapa de
# download e a de análise têm dependências distintas, e o split permite iterar
# na fórmula de scoring reaproveitando a amostra já baixada via -resume.
# ──────────────────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="BWA + Samtools + Seqtk para análise da amostra piloto"
LABEL version="1.0"

ARG DEBIAN_FRONTEND=noninteractive
ARG BWA_VERSION=0.7.18
ARG SAMTOOLS_VERSION=1.19.2
ARG SEQTK_VERSION=1.4

# Dependências de compilação e de runtime
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        wget \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instala BWA (mapeamento da amostra piloto contra a referência mitocondrial)
RUN wget -q \
    "https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz" && \
    tar -xzf "v${BWA_VERSION}.tar.gz" && \
    make -C "bwa-${BWA_VERSION}" && \
    mkdir -p /opt/bwa && \
    mv "bwa-${BWA_VERSION}/bwa" /opt/bwa/ && \
    rm -rf "bwa-${BWA_VERSION}" "v${BWA_VERSION}.tar.gz"

# Instala Samtools (profundidade e estatísticas de mapeamento)
# --without-curses: só usamos depth/flagstat/sort/view, não a interface tview
RUN wget -q \
    "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xjf "samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd "samtools-${SAMTOOLS_VERSION}" && \
    ./configure --prefix=/opt/samtools --without-curses && \
    make && \
    make install && \
    cd .. && \
    rm -rf "samtools-${SAMTOOLS_VERSION}" "samtools-${SAMTOOLS_VERSION}.tar.bz2"

# Instala Seqtk (subamostragem e manipulação de FASTQ)
RUN wget -q \
    "https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz" && \
    tar -xzf "v${SEQTK_VERSION}.tar.gz" && \
    make -C "seqtk-${SEQTK_VERSION}" && \
    mkdir -p /opt/seqtk && \
    mv "seqtk-${SEQTK_VERSION}/seqtk" /opt/seqtk/ && \
    rm -rf "seqtk-${SEQTK_VERSION}" "v${SEQTK_VERSION}.tar.gz"

ENV PATH="/opt/bwa:/opt/samtools/bin:/opt/seqtk:$PATH"

WORKDIR /data

CMD ["bwa"]
