ARG CCFINDER_VERSION=4.2.30
ARG CCFINDER_REF=release-4.2.20
ARG MACSYFINDER_VERSION=1.0.5
ARG VMATCH_URL=http://vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
ARG SINGULARITY_INSPECTOR_IMAGE=quay.io/singularity/singularity:v3.8.2-slim@sha256:2558c93c63aa8f721aae106ff1052dff5442171d7c00aed3e81d0fe62172b118

FROM ubuntu:18.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG CCFINDER_VERSION
ARG CCFINDER_REF
ARG MACSYFINDER_VERSION
ARG VMATCH_URL
ARG SINGULARITY_INSPECTOR_IMAGE

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

COPY CRISPRCasFinder.runtime.patch /tmp/CRISPRCasFinder.runtime.patch
COPY macsyfinder-1.0.5.tar.gz /tmp/macsyfinder-1.0.5.tar.gz

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
        libdate-calc-perl \
        libdatetime-perl \
        libdigest-md5-perl \
        libexpat1-dev \
        libxml-parser-perl \
        libxml-sax-expat-perl \
        libxml-simple-perl \
        muscle \
        ncbi-blast+ \
        parallel \
        patch \
        perl \
        prodigal \
        python \
        wget; \
    rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    git clone --branch "${CCFINDER_REF}" --depth 1 https://github.com/dcouvin/CRISPRCasFinder.git /usr/local/CRISPRCasFinder; \
    patch -d /usr/local/CRISPRCasFinder -p0 < /tmp/CRISPRCasFinder.runtime.patch; \
    install -d /usr/local/CRISPRCasFinder/bin /usr/local/CRISPRCasFinder/CasFinder-2.0.4; \
    ln -sf ../CasFinder-2.0.3/CASprofiles-2.0.3 /usr/local/CRISPRCasFinder/CasFinder-2.0.4/CASprofiles-2.0.4; \
    ln -sf ../CasFinder-2.0.3/DEF-Class-2.0.3 /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-Class-2.0.4; \
    ln -sf ../CasFinder-2.0.3/DEF-SubTyping-2.0.3 /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-SubTyping-2.0.4; \
    ln -sf ../CasFinder-2.0.3/DEF-Typing-2.0.3 /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-Typing-2.0.4; \
    rm -rf /usr/local/CRISPRCasFinder/.git; \
    chmod 755 /usr/local/CRISPRCasFinder /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl

RUN set -eux; \
    curl -fsSL "${VMATCH_URL}" -o /tmp/vmatch.tar.gz; \
    tar -xzf /tmp/vmatch.tar.gz -C /tmp; \
    VMATCH_DIR="$(find /tmp -maxdepth 1 -type d -name 'vmatch-*' | head -n 1)"; \
    test -n "${VMATCH_DIR}"; \
    gcc -Wall -Werror -fPIC -O3 -shared "${VMATCH_DIR}/SELECT/sel392.c" -o /usr/local/CRISPRCasFinder/sel392v2.so; \
    install -m 0755 "${VMATCH_DIR}/vmatch" /usr/local/CRISPRCasFinder/bin/vmatch2; \
    install -m 0755 "${VMATCH_DIR}/mkvtree" /usr/local/CRISPRCasFinder/bin/mkvtree2; \
    install -m 0755 "${VMATCH_DIR}/vsubseqselect" /usr/local/CRISPRCasFinder/bin/vsubseqselect2; \
    ln -sf vmatch2 /usr/local/CRISPRCasFinder/bin/vmatch; \
    ln -sf mkvtree2 /usr/local/CRISPRCasFinder/bin/mkvtree; \
    ln -sf vsubseqselect2 /usr/local/CRISPRCasFinder/bin/vsubseqselect; \
    rm -rf /tmp/vmatch.tar.gz "${VMATCH_DIR}"

RUN set -eux; \
    tar -xzf /tmp/macsyfinder-1.0.5.tar.gz -C /tmp; \
    cd /tmp/macsyfinder-1.0.5; \
    python setup.py build; \
    python setup.py install --no-viewer; \
    rm -rf /tmp/macsyfinder-1.0.5 /tmp/macsyfinder-1.0.5.tar.gz

RUN set -eux; \
    cpanm --notest \
        Bio::DB::Fasta \
        Bio::Seq \
        Bio::SeqIO \
        Bio::AlignIO; \
    install -d /opt/perl5; \
    ensure_perl_module() { \
        local perl_module="$1"; \
        shift; \
        if ! perl -M"${perl_module}" -e1 >/dev/null 2>&1; then \
            cpanm --notest -L /opt/perl5 "$@"; \
        fi; \
    }; \
    ensure_perl_module "Date::Calc" "Date::Calc"; \
    ensure_perl_module "Digest::MD5" "Digest::MD5"; \
    ensure_perl_module "JSON::Parse" "JSON::Parse"; \
    ensure_perl_module "Try::Tiny" "Try::Tiny"; \
    ensure_perl_module "Unix::Sysexits" "Unix::Sysexits"; \
    ensure_perl_module "XML::Simple" "XML::Simple"; \
    install -d /usr/local/share/perl/5.26.1 /usr/local/lib/x86_64-linux-gnu/perl/5.26.1

FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG CCFINDER_VERSION
ARG CCFINDER_REF
ARG MACSYFINDER_VERSION
ARG VMATCH_URL
ARG SINGULARITY_INSPECTOR_IMAGE

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

LABEL org.opencontainers.image.title="CRISPRCasFinder" \
      org.opencontainers.image.description="CRISPRCasFinder runtime aligned to the inspected working Singularity image contract" \
      org.opencontainers.image.version="${CCFINDER_VERSION}" \
      org.opencontainers.image.source="https://github.com/dcouvin/CRISPRCasFinder" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.base.name="ubuntu:18.04" \
      io.nf_myco_annotation.ccfinder_ref="${CCFINDER_REF}" \
      io.nf_myco_annotation.macsyfinder_version="${MACSYFINDER_VERSION}" \
      io.nf_myco_annotation.vmatch_url="${VMATCH_URL}" \
      io.nf_myco_annotation.reference_runtime="CRISPRCasFinder.simg" \
      io.nf_myco_annotation.reference_runtime_version="4.2.30" \
      io.nf_myco_annotation.reference_runtime_inspector="${SINGULARITY_INSPECTOR_IMAGE}"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH=/usr/local/CRISPRCasFinder/bin:/opt/perl5/bin:/usr/local/bin:${PATH}
ENV PERL5LIB=/opt/perl5/lib/perl5
WORKDIR /work

COPY macsyfinder-wrapper.sh /tmp/macsyfinder-wrapper.sh
COPY ccfinder-smoke-test.sh /tmp/ccfinder-smoke-test.sh

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
        libdate-calc-perl \
        libdatetime-perl \
        libdigest-md5-perl \
        libxml-parser-perl \
        libxml-sax-expat-perl \
        libxml-simple-perl \
        muscle \
        ncbi-blast+ \
        parallel \
        perl \
        prodigal \
        python; \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/CRISPRCasFinder /usr/local/CRISPRCasFinder
COPY --from=builder /usr/local/bin/macsyfinder /usr/local/bin/macsyfinder
COPY --from=builder /usr/local/lib/python2.7 /usr/local/lib/python2.7
COPY --from=builder /usr/local/lib/x86_64-linux-gnu/perl/5.26.1 /usr/local/lib/x86_64-linux-gnu/perl/5.26.1
COPY --from=builder /usr/local/share/macsyfinder /usr/local/share/macsyfinder
COPY --from=builder /usr/local/share/perl/5.26.1 /usr/local/share/perl/5.26.1
COPY --from=builder /etc/macsyfinder /etc/macsyfinder
COPY --from=builder /opt/perl5 /opt/perl5

RUN set -eux; \
    install -m 0755 /tmp/macsyfinder-wrapper.sh /usr/local/CRISPRCasFinder/macsyfinder; \
    install -m 0755 /tmp/ccfinder-smoke-test.sh /usr/local/bin/ccfinder-smoke-test; \
    ln -sf /usr/local/CRISPRCasFinder/bin/vmatch2 /usr/local/bin/vmatch2; \
    ln -sf /usr/local/CRISPRCasFinder/bin/vmatch /usr/local/bin/vmatch; \
    ln -sf /usr/local/CRISPRCasFinder/bin/mkvtree2 /usr/local/bin/mkvtree2; \
    ln -sf /usr/local/CRISPRCasFinder/bin/mkvtree /usr/local/bin/mkvtree; \
    ln -sf /usr/local/CRISPRCasFinder/bin/vsubseqselect2 /usr/local/bin/vsubseqselect2; \
    ln -sf /usr/local/CRISPRCasFinder/bin/vsubseqselect /usr/local/bin/vsubseqselect; \
    test -f /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; \
    test -f /usr/local/CRISPRCasFinder/sel392v2.so; \
    test -d /usr/local/CRISPRCasFinder/CasFinder-2.0.4; \
    test -L /usr/local/CRISPRCasFinder/CasFinder-2.0.4/CASprofiles-2.0.4; \
    test -L /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-Class-2.0.4; \
    test -L /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-SubTyping-2.0.4; \
    test -L /usr/local/CRISPRCasFinder/CasFinder-2.0.4/DEF-Typing-2.0.4; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/CRISPR_crisprdb.csv; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv; \
    test -f /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv; \
    test -x /usr/local/CRISPRCasFinder/macsyfinder; \
    test -x /usr/local/bin/ccfinder-smoke-test; \
    test -f /usr/local/share/perl/5.26.1/Bio/DB/IndexedBase.pm; \
    test -f /usr/local/share/perl/5.26.1/Bio/DB/Fasta.pm; \
    command -v python; \
    command -v vmatch2; \
    command -v mkvtree2; \
    command -v vsubseqselect2; \
    command -v muscle; \
    command -v prodigal; \
    command -v macsyfinder; \
    grep -F 'my $version = "4.2.30";' /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; \
    grep -F 'my $PathVmatch2 = "/usr/local/CRISPRCasFinder/bin/vmatch2";' /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; \
    grep -F 'my $PathMacsyfinder = "/usr/local/CRISPRCasFinder/macsyfinder";' /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; \
    test -z "$(perl -MBio::DB::IndexedBase -e 1 2>&1)"; \
    python -c 'import macsypy'; \
    perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -v >/tmp/ccfinder-version.txt; \
    /usr/local/bin/ccfinder-smoke-test
