# ──────────────────────────────────────────────────────────────────────────────
# FastQC — Controle de qualidade de dados NGS
# Versão: 0.12.1  |  Base: ubuntu:22.04
# ──────────────────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="FastQC para avaliação de qualidade de leituras NGS"
LABEL version="0.12.1"

ARG DEBIAN_FRONTEND=noninteractive
ARG FASTQC_VERSION=0.12.1

# Java (obrigatório para FastQC) + wget + unzip
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
        wget \
        unzip \
        perl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instala FastQC
RUN wget -q \
    "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" && \
    unzip "fastqc_v${FASTQC_VERSION}.zip" && \
    mv FastQC /opt/fastqc && \
    chmod +x /opt/fastqc/fastqc && \
    rm "fastqc_v${FASTQC_VERSION}.zip"

ENV PATH="/opt/fastqc:$PATH"

WORKDIR /data

CMD ["fastqc", "--version"]
