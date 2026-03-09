#!/usr/bin/env bash
set -euo pipefail

tmpdir="$(mktemp -d)"
cleanup() {
    rm -rf "${tmpdir}"
}
trap cleanup EXIT

task_root="${tmpdir}"
run_root="${task_root}/ccfinder_run"
tool_output_root="${task_root}/ccfinder_raw"
input_source="/usr/local/CRISPRCasFinder/install_test/sequence.fasta"
staged_input="${run_root}/sequence.fasta"
run_log="${task_root}/ccfinder-smoke.log"

mkdir -p "${run_root}"
cp "${input_source}" "${staged_input}"

awk -v output_root="${task_root}" '
    function flush_record(    output_path, start, chunk) {
        if (contig_id == "") {
            return
        }
        output_path = output_root "/" contig_id ".fna"
        print ">" header > output_path
        for (start = 1; start <= length(sequence); start += 60) {
            chunk = substr(sequence, start, 60)
            print chunk >> output_path
        }
        close(output_path)
    }
    /^>/ {
        flush_record()
        header = substr($0, 2)
        contig_id = $1
        sub(/^>/, "", contig_id)
        sub(/\.[0-9]+$/, "", contig_id)
        sequence = ""
        next
    }
    {
        gsub(/[[:space:]]/, "", $0)
        sequence = sequence $0
    }
    END {
        flush_record()
    }
' "${staged_input}"

(
    cd "${run_root}"
    : > index.html
    perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl \
        -in "$(basename "${staged_input}")" \
        -outdir "${tool_output_root}" \
        -soFile /usr/local/CRISPRCasFinder/sel392v2.so \
        -DBcrispr /usr/local/CRISPRCasFinder/supplementary_files/CRISPR_crisprdb.csv \
        -repeats /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv \
        -DIRrepeat /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
        -cpuMacSyFinder 1 \
        -cpuProkka 1 \
        -log \
        -html \
        -levelMin 2 \
        -cas \
        -ccvRep \
        -getSummaryCasfinder \
        -gcode 11 \
        -quiet
) >"${run_log}" 2>&1

find "${tool_output_root}" -type f -name 'CRISPR-Cas_summary.tsv' | grep -q .
find "${tool_output_root}" -type f -name 'Casfinder_summary_*.tsv' | grep -q .
find "${tool_output_root}" -type f -name 'rawCas.fna' | grep -q .
find "${tool_output_root}" -type f -name 'result.json' | grep -q .
