# ──────────────────────────────────────────────────────────────────────────────
# Trim Galore — Remoção de adaptadores e trimming de qualidade
# Versão: 0.6.10  |  Dependências: Cutadapt 4.6 + FastQC 0.12.1
# Base: ubuntu:22.04
# ──────────────────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="Trim Galore para remoção de adaptadores e controle de qualidade"
LABEL version="0.6.10"

ARG DEBIAN_FRONTEND=noninteractive
ARG FASTQC_VERSION=0.12.1
ARG TRIMGALORE_VERSION=0.6.10

# Dependências: Java (FastQC), Python3 (Cutadapt), Perl, pigz (compressão paralela)
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
        python3 \
        python3-pip \
        python3-setuptools \
        perl \
        wget \
        unzip \
        pigz && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instala Cutadapt (dependência do Trim Galore)
RUN pip3 install --no-cache-dir cutadapt==4.6

# Instala FastQC (necessário para relatórios do Trim Galore)
RUN wget -q \
    "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" && \
    unzip "fastqc_v${FASTQC_VERSION}.zip" && \
    mv FastQC /opt/fastqc && \
    chmod +x /opt/fastqc/fastqc && \
    rm "fastqc_v${FASTQC_VERSION}.zip"

# Instala Trim Galore
RUN wget -q \
    "https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/${TRIMGALORE_VERSION}.tar.gz" && \
    tar -xzf "${TRIMGALORE_VERSION}.tar.gz" && \
    mv "TrimGalore-${TRIMGALORE_VERSION}" /opt/trim-galore && \
    chmod +x /opt/trim-galore/trim_galore && \
    rm "${TRIMGALORE_VERSION}.tar.gz"

ENV PATH="/opt/fastqc:/opt/trim-galore:$PATH"

WORKDIR /data

CMD ["trim_galore", "--version"]
