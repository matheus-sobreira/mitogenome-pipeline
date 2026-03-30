# ──────────────────────────────────────────────────────────────────────────────
# NOVOPlasty — Montagem de genomas de organelas (seed-and-extend)
# Versão: 4.3.1  |  Base: ubuntu:22.04
# ──────────────────────────────────────────────────────────────────────────────

FROM ubuntu:22.04

LABEL maintainer="Matheus Sobreira Benevides"
LABEL description="NOVOPlasty para montagem de genomas mitocondriais"
LABEL version="4.3.1"

ARG DEBIAN_FRONTEND=noninteractive
ARG NOVOPLASTY_VERSION=4.3.1

# NOVOPlasty é um script Perl — somente Perl é necessário
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        perl \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instala NOVOPlasty
# O padrão de tag no GitHub é 'NOVOPlasty<versão>' (ex: NOVOPlasty4.3.1)
RUN wget -q \
    "https://github.com/ndierckx/NOVOPlasty/archive/refs/tags/NOVOPlasty${NOVOPLASTY_VERSION}.tar.gz" \
    -O novoplasty.tar.gz && \
    tar -xzf novoplasty.tar.gz && \
    mv NOVOPlasty-NOVOPlasty${NOVOPLASTY_VERSION} /opt/novoplasty && \
    rm novoplasty.tar.gz && \
    chmod +x /opt/novoplasty/NOVOPlasty4.3.pl

ENV PATH="/opt/novoplasty:$PATH"

WORKDIR /data

CMD ["perl", "/opt/novoplasty/NOVOPlasty4.3.pl"]
