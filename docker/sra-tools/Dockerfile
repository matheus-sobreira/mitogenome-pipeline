# ──────────────────────────────────────────────────────────────────────────────
# SRA-Toolkit — Download de dados do NCBI SRA
# Versão: 3.0.10  |  Base: ubuntu:22.04
# ──────────────────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="SRA-Toolkit para download de dados do NCBI SRA"
LABEL version="3.0.10"

ARG DEBIAN_FRONTEND=noninteractive
ARG SRA_VERSION=3.0.10

# Dependências mínimas
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        wget \
        ca-certificates \
        libxml2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instala SRA-Toolkit
RUN wget -q \
    "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" && \
    tar -xzf "sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" && \
    mv "sratoolkit.${SRA_VERSION}-ubuntu64" /opt/sra-tools && \
    rm "sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz"

ENV PATH="/opt/sra-tools/bin:$PATH"

# Configuração inicial para evitar prompts interativos
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "d3b07384-d9a3-4f3c-b7e2-4c5c1d3abe42"\n' \
        > /root/.ncbi/user-settings.mkfg

WORKDIR /data

CMD ["fasterq-dump", "--version"]
