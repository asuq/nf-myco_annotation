FROM python:3.12-slim

ARG CODETTA_REF=863359ed326276602d44e48227b6003ac6ffd266
ARG CODETTA_VERSION=v2.0
ARG CODETTA_SOURCE_BRANCH=main

LABEL org.opencontainers.image.title="Codetta"
LABEL org.opencontainers.image.version="${CODETTA_VERSION}"
LABEL org.opencontainers.image.revision="${CODETTA_REF}"
LABEL org.opencontainers.image.source="https://github.com/kshulgina/codetta"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gzip \
        make \
        procps \
        tar \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir numpy scipy "setuptools<82" \
    && python3 -c "import pkg_resources"

WORKDIR /opt

RUN curl -fsSL "https://github.com/kshulgina/codetta/archive/${CODETTA_REF}.tar.gz" -o codetta.tar.gz \
    && tar -xzf codetta.tar.gz \
    && mv "codetta-${CODETTA_REF}" codetta \
    && rm codetta.tar.gz

WORKDIR /opt/codetta

RUN bash setup.sh \
    && chmod +x /opt/codetta/codetta.py /opt/codetta/codetta_align.py /opt/codetta/codetta_summary.py /opt/codetta/codetta_infer.py \
    && ln -s /opt/codetta/codetta.py /usr/local/bin/codetta.py \
    && ln -s /opt/codetta/codetta_align.py /usr/local/bin/codetta_align.py \
    && ln -s /opt/codetta/codetta_summary.py /usr/local/bin/codetta_summary.py \
    && ln -s /opt/codetta/codetta_infer.py /usr/local/bin/codetta_infer.py \
    && printf '%s\n' "${CODETTA_VERSION}" > /opt/codetta/.codetta_version \
    && printf '%s\n' "${CODETTA_REF}" > /opt/codetta/.codetta_source_commit \
    && printf '%s\n' "${CODETTA_SOURCE_BRANCH}" > /opt/codetta/.codetta_source_branch

ENV CODETTA_HOME=/opt/codetta
