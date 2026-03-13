FROM quay.io/biocontainers/prokka:1.15.6--pl5321hdfd78af_0

RUN MINCED_BIN="$(command -v minced)" && \
    mv "${MINCED_BIN}" "${MINCED_BIN}.real"

COPY docker/prokka/minced-wrapper.sh /usr/local/bin/minced

RUN chmod +x /usr/local/bin/minced && \
    prokka --version >/dev/null && \
    first_line="$(minced --version 2>&1 | sed -n '1p')" && \
    test "${first_line}" = "minced 0.4.2" && \
    ! printf '%s\n' "${first_line}" | grep -F 'hsperfdata' >/dev/null && \
    test -L /usr/local/bin/minced.real
