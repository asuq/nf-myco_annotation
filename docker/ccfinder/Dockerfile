ARG CCFINDER_VERSION=4.3.2
ARG CCFINDER_REF=release-4.3.2
ARG MACSYFINDER_VERSION=2.0.*
ARG CASFINDER_MODEL=CASFinder==3.1.0
ARG VMATCH_URL=http://vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz

FROM ubuntu:22.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG CCFINDER_VERSION
ARG CCFINDER_REF
ARG MACSYFINDER_VERSION
ARG CASFINDER_MODEL
ARG VMATCH_URL

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN set -eux; \
    test "$(dpkg --print-architecture)" = "amd64"; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        ca-certificates \
        cpanminus \
        curl \
        default-jre \
        emboss \
        emboss-lib \
        git \
        hmmer \
        libbio-perl-perl \
        libdatetime-perl \
        libdigest-md5-perl \
        libxml-simple-perl \
        muscle \
        ncbi-blast+ \
        parallel \
        perl \
        prodigal \
        python3 \
        python3-pip \
        python3-venv \
        wget; \
    rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    git clone --branch "${CCFINDER_REF}" --depth 1 https://github.com/dcouvin/CRISPRCasFinder.git /usr/local/CRISPRCasFinder; \
    rm -rf /usr/local/CRISPRCasFinder/.git; \
    chmod 755 /usr/local/CRISPRCasFinder /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl

RUN set -eux; \
    curl -fsSL "${VMATCH_URL}" -o /tmp/vmatch.tar.gz; \
    tar -xzf /tmp/vmatch.tar.gz -C /tmp; \
    VMATCH_DIR="$(find /tmp -maxdepth 1 -type d -name 'vmatch-*' | head -n 1)"; \
    test -n "${VMATCH_DIR}"; \
    install -d /opt/vmatch/bin; \
    gcc -Wall -Werror -fPIC -O3 -shared "${VMATCH_DIR}/SELECT/sel392.c" -o /usr/local/CRISPRCasFinder/sel392v2.so; \
    install -m 0755 "${VMATCH_DIR}/vmatch" /opt/vmatch/bin/vmatch; \
    install -m 0755 "${VMATCH_DIR}/mkvtree" /opt/vmatch/bin/mkvtree; \
    install -m 0755 "${VMATCH_DIR}/vsubseqselect" /opt/vmatch/bin/vsubseqselect; \
    rm -rf /tmp/vmatch.tar.gz "${VMATCH_DIR}"

RUN set -eux; \
    python3 -m venv /opt/macsyfinder-venv; \
    install -d /opt/macsyfinder-venv/share/macsyfinder/models; \
    /opt/macsyfinder-venv/bin/pip install --no-cache-dir --upgrade pip setuptools wheel; \
    /opt/macsyfinder-venv/bin/pip install --no-cache-dir "macsyfinder==${MACSYFINDER_VERSION}" "gitpython>=3,<4"; \
    /opt/macsyfinder-venv/bin/macsydata install --target /opt/macsyfinder-venv/share/macsyfinder/models "${CASFINDER_MODEL}"

RUN set -eux; \
    install -d /opt/perl5; \
    ensure_perl_module() { \
        local perl_module="$1"; \
        shift; \
        if ! perl -M"${perl_module}" -e1 >/dev/null 2>&1; then \
            cpanm --notest -L /opt/perl5 "$@"; \
        fi; \
    }; \
    ensure_perl_module "Bio::AlignIO" "Bio::AlignIO"; \
    ensure_perl_module "Bio::SeqIO" "Bio::SeqIO"; \
    ensure_perl_module "Bio::DB::Fasta" "Bio::DB::Fasta"; \
    ensure_perl_module "XML::Simple" "XML::Simple"; \
    ensure_perl_module "Digest::MD5" "Digest::MD5"

FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG CCFINDER_VERSION
ARG CCFINDER_REF
ARG MACSYFINDER_VERSION
ARG CASFINDER_MODEL
ARG VMATCH_URL

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

LABEL org.opencontainers.image.title="CRISPRCasFinder" \
      org.opencontainers.image.description="CRISPRCasFinder v4.3.2 runtime for nf-myco_annotation" \
      org.opencontainers.image.version="${CCFINDER_VERSION}" \
      org.opencontainers.image.source="https://github.com/dcouvin/CRISPRCasFinder" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.base.name="ubuntu:22.04" \
      io.nf_myco_annotation.ccfinder_ref="${CCFINDER_REF}" \
      io.nf_myco_annotation.macsyfinder_version="${MACSYFINDER_VERSION}" \
      io.nf_myco_annotation.casfinder_model="${CASFINDER_MODEL}" \
      io.nf_myco_annotation.vmatch_url="${VMATCH_URL}"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH=/opt/perl5/bin:/opt/macsyfinder-venv/bin:/usr/local/bin:${PATH}
ENV PERL5LIB=/opt/perl5/lib/perl5
WORKDIR /work

RUN set -eux; \
    test "$(dpkg --print-architecture)" = "amd64"; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        default-jre \
        emboss \
        emboss-lib \
        hmmer \
        libbio-perl-perl \
        libdatetime-perl \
        libdigest-md5-perl \
        libxml-simple-perl \
        muscle \
        ncbi-blast+ \
        parallel \
        perl \
        prodigal \
        python3; \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/CRISPRCasFinder /usr/local/CRISPRCasFinder
COPY --from=builder /opt/macsyfinder-venv /opt/macsyfinder-venv
COPY --from=builder /opt/perl5 /opt/perl5
COPY --from=builder /opt/vmatch /opt/vmatch

RUN set -eux; \
    ln -sf /opt/vmatch/bin/vmatch /usr/local/bin/vmatch; \
    ln -sf /opt/vmatch/bin/mkvtree /usr/local/bin/mkvtree; \
    ln -sf /opt/vmatch/bin/vsubseqselect /usr/local/bin/vsubseqselect; \
    test -f /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; \
    test -f /usr/local/CRISPRCasFinder/sel392v2.so; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/CRISPR_crisprdb.csv; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv; \
    test -d /usr/local/CRISPRCasFinder/CasFinder-2.0.3; \
    find /opt/macsyfinder-venv/share/macsyfinder/models -maxdepth 4 -type d -name CASFinder | grep -q .; \
    command -v vmatch; \
    command -v mkvtree; \
    command -v vsubseqselect; \
    command -v muscle; \
    command -v prodigal; \
    command -v hmmsearch; \
    command -v blastn; \
    command -v macsyfinder; \
    perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -v >/tmp/ccfinder-version.txt
