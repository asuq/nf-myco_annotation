FROM quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2

RUN set -eu; \
    script_path="$(command -v download_eggnog_data.py)"; \
    python - "${script_path}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
source = path.read_text(encoding="utf-8")
source = source.replace(
    "http://eggnogdb.embl.de/download/emapperdb-",
    "http://eggnog5.embl.de/download/emapperdb-",
)
source = source.replace(
    "http://eggnogdb.embl.de/download/novel_fams-",
    "http://eggnog5.embl.de/download/novel_fams-",
)
path.write_text(source, encoding="utf-8")
PY

RUN set -eu; \
    script_path="$(command -v download_eggnog_data.py)"; \
    ! grep -F 'http://eggnogdb.embl.de/download/emapperdb-' "${script_path}"; \
    ! grep -F 'http://eggnogdb.embl.de/download/novel_fams-' "${script_path}"; \
    grep -F 'http://eggnog5.embl.de/download/emapperdb-' "${script_path}" >/dev/null; \
    grep -F 'http://eggnog5.embl.de/download/novel_fams-' "${script_path}" >/dev/null
