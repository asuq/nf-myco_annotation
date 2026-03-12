FROM quay.io/biocontainers/padloc:2.0.0--hdfd78af_1

ENV PADLOC_DATA_DIR=/opt/padloc-data

RUN PADLOC_BIN="$(command -v padloc)" && \
    sed -i "s#mkdir -p \"\${SRC_DIR}/../data\"#mkdir -p \"${PADLOC_DATA_DIR}\"#" "${PADLOC_BIN}" && \
    sed -i "s#DATA=\$(normpath \"\${SRC_DIR}/../data\")#DATA=\$(normpath \"${PADLOC_DATA_DIR}\")#" "${PADLOC_BIN}" && \
    mkdir -p "${PADLOC_DATA_DIR}" && \
    padloc --data "${PADLOC_DATA_DIR}" --db-update

RUN PADLOC_BIN="$(command -v padloc)" && \
    grep -F "mkdir -p \"${PADLOC_DATA_DIR}\"" "${PADLOC_BIN}" >/dev/null && \
    grep -F "DATA=\$(normpath \"${PADLOC_DATA_DIR}\")" "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'mkdir -p "${SRC_DIR}/../data"' "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'DATA=$(normpath "${SRC_DIR}/../data")' "${PADLOC_BIN}" >/dev/null && \
    test -f "${PADLOC_DATA_DIR}/hmm/padlocdb.hmm" && \
    padloc --version >/dev/null
