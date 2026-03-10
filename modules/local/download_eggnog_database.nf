/*
 * Reuse or download the eggNOG database into its final directory.
 */
process DOWNLOAD_EGGNOG_DATABASE {
    tag "eggnog"
    label 'process_high'

    input:
    tuple val(destination), val(download_enabled), val(force)

    output:
    tuple val('eggnog'), val(destination), path('prep_mode.txt'), val('native_download'), path('lineages.txt'), emit: finalise_input
    path 'versions.yml', emit: versions

    script:
    """
    destination_path="${destination}"
    mode="prepared"
    expected_files=("eggnog.db" "eggnog_proteins.dmnd")

    has_expected_files() {
        for name in "\${expected_files[@]}"; do
            if [[ ! -f "\$1/\${name}" ]]; then
                return 1
            fi
        done
        return 0
    }

    if [[ -e "\${destination_path}" && ! -d "\${destination_path}" ]]; then
        echo "eggNOG destination must be a directory: \${destination_path}" >&2
        exit 1
    fi

    if [[ -d "\${destination_path}" && -f "\${destination_path}/.nf_myco_ready.json" ]] && has_expected_files "\${destination_path}"; then
        mode="reuse"
    else
        if [[ -d "\${destination_path}" ]] && has_expected_files "\${destination_path}"; then
            mode="prepared"
        else
            if [[ -d "\${destination_path}" ]]; then
                if [[ "${force}" != "true" ]]; then
                    echo "eggNOG destination exists but is not ready: \${destination_path}" >&2
                    exit 1
                fi
                rm -rf "\${destination_path}"
            fi

            if [[ "${download_enabled}" != "true" ]]; then
                echo "eggNOG destination is missing and download is disabled: \${destination_path}" >&2
                exit 1
            fi

            mkdir -p "\${destination_path}"

            script_path="\$(command -v download_eggnog_data.py)"
            patched_script="${PWD}/download_eggnog_data.py.patched"

            python - "\${script_path}" "\${patched_script}" <<'PY'
from pathlib import Path
import sys

source = Path(sys.argv[1]).read_text(encoding="utf-8")
source = source.replace(
    "http://eggnogdb.embl.de/download/emapperdb-",
    "http://eggnog5.embl.de/download/emapperdb-",
)
source = source.replace(
    "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level",
    "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level",
)
source = source.replace(
    "http://eggnogdb.embl.de/download/novel_fams-",
    "http://eggnog5.embl.de/download/novel_fams-",
)
Path(sys.argv[2]).write_text(source, encoding="utf-8")
PY
            chmod +x "\${patched_script}"
            python "\${patched_script}" --data_dir "\${destination_path}" -y ${force ? '-f' : ''}

            if ! has_expected_files "\${destination_path}"; then
                echo "eggNOG download did not populate the required files: \${destination_path}" >&2
                find "\${destination_path}" -maxdepth 2 -type f | sort >&2
                exit 1
            fi
        fi
    fi

    printf '%s\n' "\${mode}" > prep_mode.txt
    : > lineages.txt

    {
        printf '"%s":\n' "${task.process}"
        printf '  eggnog_mapper: "%s"\n' "\$(emapper.py --version 2>&1 | head -n 1 || echo NA)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p "${destination}"
    : > "${destination}/eggnog.db"
    : > "${destination}/eggnog_proteins.dmnd"
    : > "${destination}/.nf_myco_ready.json"
    printf '%s\n' 'reuse' > prep_mode.txt
    : > lineages.txt
    {
        printf '"%s":\n' "${task.process}"
        printf '  eggnog_mapper: "%s"\n' 'stub'
    } > versions.yml
    """
}
