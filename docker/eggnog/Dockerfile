FROM quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2

USER root

RUN apt-get update \
    && apt-get install --yes --no-install-recommends pigz \
    && rm -rf /var/lib/apt/lists/*

COPY docker/eggnog/gunzip /usr/local/bin/gunzip

RUN chmod +x /usr/local/bin/gunzip
