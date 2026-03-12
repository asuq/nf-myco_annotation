FROM quay.io/biocontainers/padloc:2.0.0--hdfd78af_1

RUN PADLOC_BIN="$(command -v padloc)" && \
    sed -i 's#mkdir -p "${SRC_DIR}/../data"#mkdir -p "/tmp/padloc-launcher-data"#' "${PADLOC_BIN}" && \
    sed -i 's#DATA=$(normpath "${SRC_DIR}/../data")#DATA=$(normpath "/tmp/padloc-launcher-data")#' "${PADLOC_BIN}" && \
    mkdir -p /tmp/padloc-launcher-data

RUN PADLOC_BIN="$(command -v padloc)" && \
    grep -F 'mkdir -p "/tmp/padloc-launcher-data"' "${PADLOC_BIN}" >/dev/null && \
    grep -F 'DATA=$(normpath "/tmp/padloc-launcher-data")' "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'mkdir -p "${SRC_DIR}/../data"' "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'DATA=$(normpath "${SRC_DIR}/../data")' "${PADLOC_BIN}" >/dev/null
