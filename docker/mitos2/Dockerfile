# ──────────────────────────────────────────────────────────────────────────────
# MITOS2 — Anotação funcional de genomas mitocondriais
# Versão: 2.1.9  |  Base: continuumio/miniconda3:24.1.2-0
#
# MITOS2 anota genes proteicos, rRNAs e tRNAs de mitogenomas.
# Dependências instaladas via bioconda: blast, hmmer, infernal, viennarna, mafft.
# ──────────────────────────────────────────────────────────────────────────────

FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="MITOS2 para anotação funcional de genomas mitocondriais"
LABEL version="2.1.9"

ARG MITOS_VERSION=2.1.9

ARG DEBIAN_FRONTEND=noninteractive

# Instala mamba para resolução de dependências mais rápida
RUN conda install -y -n base -c conda-forge mamba=1.5.8 && \
    conda clean -afy

# Instala MITOS2 e todas as suas dependências via bioconda
RUN mamba install -y -n base \
        -c bioconda \
        -c conda-forge \
        mitos=${MITOS_VERSION} && \
    conda clean -afy && \
    find /opt/conda -name "*.pyc" -delete && \
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

WORKDIR /data

CMD ["runmitos.py", "--help"]
