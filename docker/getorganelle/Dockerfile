# ──────────────────────────────────────────────────────────────────────────────
# GetOrganelle — Montador de genomas de organela (sondagem, Eixo F)
# Versão: 1.7.7.1  |  Base: continuumio/miniconda3
#
# ⚠️ IMAGEM DE SONDAGEM, não integrada ao pipeline. Existe para responder três
# perguntas (ver DEC-17 no vault):
#   1. o GetOrganelle tem janela de intake limitada por memória, como o
#      NOVOPlasty (1,691 M reads/GB)?
#   2. qual o limiar de cobertura dele, comparado aos ~36–48× do NOVOPlasty?
#   3. a seed COX1 e o mitogenoma congênere funcionam como `-s`?
#
# Motivo do interesse: licença GPL-3.0 (redistribuível), enquanto o NOVOPlasty
# proíbe distribuição de derivados e bloqueia o nf-core. Faz SÓ montagem, então
# é substituição direta preservando o MITOS2 e o pós-processamento.
# ──────────────────────────────────────────────────────────────────────────────

FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="GetOrganelle para sondagem de montador alternativo ao NOVOPlasty"
LABEL version="1.7.7.1"

ARG GETORGANELLE_VERSION=1.7.7.1

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

RUN conda install -y -c bioconda -c conda-forge \
        "getorganelle=${GETORGANELLE_VERSION}" && \
    conda clean -afy

# Banco de seeds/labels do modo animal_mt. Sem ele o GetOrganelle não classifica
# os contigs, mesmo quando a seed é fornecida via -s.
RUN get_organelle_config.py --add animal_mt

WORKDIR /data

CMD ["get_organelle_from_reads.py", "--help"]
