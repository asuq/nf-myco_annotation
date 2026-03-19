"""Regression tests for parser-safe Nextflow module scripts."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULES_DIR = ROOT / "modules" / "local"
RAW_COMMAND_SUBSTITUTION = re.compile(r'(?<!\\)\$\(')
SOFT_FAIL_RETRY_MARKERS = {
    "barrnap.nf": "retrying_barrnap",
    "busco.nf": "retrying_busco",
    "ccfinder.nf": "retrying_ccfinder",
    "checkm2.nf": "retrying_checkm2",
    "codetta.nf": "retrying_codetta",
    "eggnog.nf": "retrying_eggnog",
    "padloc.nf": "retrying_padloc",
    "prokka.nf": "retrying_prokka",
}


def find_raw_command_substitutions(path: Path) -> list[str]:
    """Return line-numbered raw shell substitutions that would break parsing."""
    matches: list[str] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if RAW_COMMAND_SUBSTITUTION.search(line):
            matches.append(f"{path.relative_to(ROOT)}:{line_number}:{line.strip()}")
    return matches


class NextflowModuleSyntaxTestCase(unittest.TestCase):
    """Protect local Nextflow modules from raw shell substitution syntax."""

    def test_local_modules_escape_shell_command_substitutions(self) -> None:
        """Require parser-safe escaping for shell command substitution tokens."""
        matches: list[str] = []
        for path in sorted(MODULES_DIR.glob("*.nf")):
            matches.extend(find_raw_command_substitutions(path))

        self.assertEqual(matches, [])

    def test_assign_gcode_and_qc_stages_unique_checkm2_input_names(self) -> None:
        """Require unique staged names for paired CheckM2 reports."""
        module_path = MODULES_DIR / "assign_gcode_and_qc.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn("name: 'checkm2_gcode4_report.tsv'", module_text)
        self.assertIn("name: 'checkm2_gcode11_report.tsv'", module_text)

    def test_checkm2_resolves_directory_parameters_to_dmnd_files(self) -> None:
        """Require CheckM2 to accept one DMND directory and retry transient failures."""
        module_path = MODULES_DIR / "checkm2.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn("tuple val(meta), path(genome), path(checkm2_db)", module_text)
        self.assertIn('database_path="${checkm2_db}"', module_text)
        self.assertIn('dmnd_candidates=("\\${database_path}"/*.dmnd)', module_text)
        self.assertIn('--database_path "\\${database_path}"', module_text)
        self.assertIn('max_attempts="${params.soft_fail_attempts}"', module_text)
        self.assertIn('while (( attempt <= max_attempts ))', module_text)
        self.assertIn("retrying_checkm2", module_text)

    def test_runtime_database_dirs_are_staged_into_checkm2_codetta_and_eggnog(self) -> None:
        """Require external database directories to enter runtime modules as path inputs."""
        main_text = (ROOT / "main.nf").read_text(encoding="utf-8")
        per_sample_qc_text = (
            ROOT / "subworkflows" / "local" / "per_sample_qc.nf"
        ).read_text(encoding="utf-8")
        per_sample_annotation_text = (
            ROOT / "subworkflows" / "local" / "per_sample_annotation.nf"
        ).read_text(encoding="utf-8")
        codetta_text = (MODULES_DIR / "codetta.nf").read_text(encoding="utf-8")
        eggnog_text = (MODULES_DIR / "eggnog.nf").read_text(encoding="utf-8")

        self.assertIn("checkm2Db = Channel.fromPath(params.checkm2_db, checkIfExists: true)", main_text)
        self.assertIn("codettaDb = Channel.fromPath(params.codetta_db, checkIfExists: true)", main_text)
        self.assertIn("eggnogDb = Channel.fromPath(params.eggnog_db, checkIfExists: true)", main_text)
        self.assertIn("PER_SAMPLE_QC(", main_text)
        self.assertIn("checkm2Db,", main_text)
        self.assertIn("PER_SAMPLE_ANNOTATION(", main_text)
        self.assertIn("codettaDb,", main_text)
        self.assertIn("eggnogDb,", main_text)
        self.assertIn("take:\n    sample_genomes\n    checkm2_db\n    busco_datasets", per_sample_qc_text)
        self.assertIn(".combine(checkm2_db)", per_sample_qc_text)
        self.assertIn(
            "take:\n    sample_genomes\n    gcode_summaries\n    codetta_db\n    eggnog_db",
            per_sample_annotation_text,
        )
        self.assertIn("CODETTA(sample_genomes.combine(codetta_db))", per_sample_annotation_text)
        self.assertIn("SUMMARISE_CODETTA(CODETTA.out.summary_input)", per_sample_annotation_text)
        self.assertIn("EGGNOG(eggnog_inputs.combine(eggnog_db))", per_sample_annotation_text)
        self.assertIn("tuple val(meta), path(genome), path(codetta_db)", codetta_text)
        self.assertIn('{ "${params.outdir}/samples/${meta.accession}" }', codetta_text)
        self.assertIn('cp -R "${codetta_db}/". "\\${resource_directory}/"', codetta_text)
        self.assertIn("tuple val(meta), path(faa), path(eggnog_db)", eggnog_text)
        self.assertIn('--data_dir "${eggnog_db}"', eggnog_text)

    def test_summarise_busco_emits_lineage_specific_summary_names(self) -> None:
        """Require unique BUSCO summary filenames per lineage."""
        module_path = MODULES_DIR / "summarise_busco.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('path("busco_summary_${lineage}.tsv")', module_text)
        self.assertIn('--output "busco_summary_${lineage}.tsv"', module_text)

    def test_download_busco_dataset_preserves_lineage_directory_name(self) -> None:
        """Require downloaded BUSCO datasets to be staged under their lineage names."""
        module_path = MODULES_DIR / "download_busco_dataset.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('tuple val(lineage), path("${lineage}"), emit: dataset', module_text)
        self.assertIn('ln -s "\\${dataset_dir}" "${lineage}"', module_text)

    def test_busco_stages_offline_datasets_under_download_root(self) -> None:
        """Require BUSCO offline runs to stage lineage datasets and retry transient failures."""
        module_path = MODULES_DIR / "busco.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('busco_download_root="busco_downloads"', module_text)
        self.assertIn('staged_lineage_dir="\\${busco_download_root}/lineages/${lineage}"', module_text)
        self.assertIn('dataset_source="\\$(cd "${dataset_dir}" && pwd)"', module_text)
        self.assertIn('ln -s "\\${dataset_source}" "\\${staged_lineage_dir}"', module_text)
        self.assertIn('--download_path "\\${busco_download_root}"', module_text)
        self.assertIn('--lineage_dataset "${lineage}"', module_text)
        self.assertIn('max_attempts="${params.soft_fail_attempts}"', module_text)
        self.assertIn('while (( attempt <= max_attempts ))', module_text)
        self.assertIn('retrying_busco', module_text)

    def test_soft_fail_modules_use_internal_retry_knob(self) -> None:
        """Require soft-fail runtime modules to use the shared internal retry knob."""
        for module_name, retry_marker in SOFT_FAIL_RETRY_MARKERS.items():
            module_text = (MODULES_DIR / module_name).read_text(encoding="utf-8")

            self.assertIn('max_attempts="${params.soft_fail_attempts}"', module_text, module_name)
            self.assertIn('while (( attempt <= max_attempts ))', module_text, module_name)
            self.assertIn(retry_marker, module_text, module_name)

    def test_busco_consumers_stage_unique_summary_names(self) -> None:
        """Require BUSCO summary consumers to stage input files uniquely."""
        expected_snippet = "path busco_tables, name: 'busco_tables/busco_table??.tsv'"
        for module_name in (
            "build_fastani_inputs.nf",
            "build_master_table.nf",
            "write_sample_status.nf",
        ):
            module_path = MODULES_DIR / module_name
            module_text = module_path.read_text(encoding="utf-8")
            self.assertIn(expected_snippet, module_text, module_name)

    def test_assembly_stats_flow_is_staged_into_ani_and_final_outputs(self) -> None:
        """Require the in-house assembly-stats TSV to feed all downstream consumers."""
        cohort_workflow = (ROOT / "subworkflows" / "local" / "cohort_ani.nf").read_text(
            encoding="utf-8"
        )
        final_outputs_workflow = (
            ROOT / "subworkflows" / "local" / "final_outputs.nf"
        ).read_text(encoding="utf-8")
        build_fastani_module = (MODULES_DIR / "build_fastani_inputs.nf").read_text(
            encoding="utf-8"
        )
        build_master_module = (MODULES_DIR / "build_master_table.nf").read_text(
            encoding="utf-8"
        )
        write_status_module = (MODULES_DIR / "write_sample_status.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("CALCULATE_ASSEMBLY_STATS(staged_manifest, staged_fasta_files)", cohort_workflow)
        self.assertIn("CALCULATE_ASSEMBLY_STATS.out.stats", cohort_workflow)
        self.assertIn("path assembly_stats", build_fastani_module)
        self.assertIn("--assembly-stats \"${assembly_stats}\"", build_fastani_module)
        self.assertIn("path assembly_stats", build_master_module)
        self.assertIn("--assembly-stats \"${assembly_stats}\"", build_master_module)
        self.assertIn("path assembly_stats", write_status_module)
        self.assertIn("--assembly-stats \"${assembly_stats}\"", write_status_module)
        self.assertIn("assembly_stats", final_outputs_workflow)

    def test_write_sample_status_runs_helper_via_python3(self) -> None:
        """Require the shared helper image to launch sample-status via Python."""
        module_path = MODULES_DIR / "write_sample_status.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('script_path="\\$(command -v build_sample_status.py)"', module_text)
        self.assertIn('python3 "\\${script_path}"', module_text)

    def test_collect_versions_runs_helper_via_python3(self) -> None:
        """Require version collection to resolve the staged helper explicitly."""
        module_path = MODULES_DIR / "collect_versions.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('script_path="\\$(command -v collect_versions.py)"', module_text)
        self.assertIn('python3 "\\${script_path}"', module_text)

    def test_summarise_codetta_runs_helper_via_python3(self) -> None:
        """Require the Codetta helper summariser to resolve the staged script explicitly."""
        module_path = MODULES_DIR / "summarise_codetta.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('script_path="\\$(command -v summarise_codetta.py)"', module_text)
        self.assertIn('python3 "\\${script_path}"', module_text)

    def test_validate_inputs_stages_sample_status_columns_asset(self) -> None:
        """Require input validation to stage the sample-status asset into the container."""
        module_text = (MODULES_DIR / "validate_inputs.nf").read_text(encoding="utf-8")
        workflow_text = (
            ROOT / "subworkflows" / "local" / "input_validation_and_staging.nf"
        ).read_text(encoding="utf-8")
        main_text = (ROOT / "main.nf").read_text(encoding="utf-8")

        self.assertIn("path sample_status_columns", module_text)
        self.assertIn('--sample-status-columns "${sample_status_columns}"', module_text)
        self.assertIn("--defer-genome-fasta-check", module_text)
        self.assertIn("sample_status_columns", workflow_text)
        self.assertIn("sampleStatusColumns = Channel.value", main_text)
        self.assertIn("INPUT_VALIDATION_AND_STAGING(sampleCsv, metadata, sampleStatusColumns)", main_text)

    def test_runtime_tool_modules_write_versions_without_indented_headers(self) -> None:
        """Require runtime tool modules to emit versions via printf, not heredoc indentation."""
        expected_modules = (
            "barrnap.nf",
            "busco.nf",
            "calculate_assembly_stats.nf",
            "checkm2.nf",
            "codetta.nf",
            "eggnog.nf",
            "padloc.nf",
            "prokka.nf",
            "stage_inputs.nf",
            "summarise_codetta.nf",
        )

        for module_name in expected_modules:
            module_text = (MODULES_DIR / module_name).read_text(encoding="utf-8")
            self.assertIn("printf ", module_text, module_name)
            self.assertNotIn("cat <<EOF > versions.yml", module_text, module_name)

    def test_stage_inputs_extracts_seqtk_version_from_version_line(self) -> None:
        """Require seqtk version probing to ignore the multi-line usage banner."""
        module_text = (MODULES_DIR / "stage_inputs.nf").read_text(encoding="utf-8")

        self.assertIn("/^Version:/ { print \\$2; exit }", module_text)
        self.assertIn('seqtk_version="\\${seqtk_version:-NA}"', module_text)

    def test_calculate_assembly_stats_extracts_seqtk_version_from_version_line(self) -> None:
        """Require assembly-stats provenance to use the same stable seqtk probe."""
        module_text = (MODULES_DIR / "calculate_assembly_stats.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("/^Version:/ { print \\$2; exit }", module_text)
        self.assertIn('seqtk_version="\\${seqtk_version:-NA}"', module_text)

    def test_barrnap_uses_the_last_non_blank_version_line(self) -> None:
        """Require the Barrnap probe to avoid its leading program-name line."""
        module_text = (MODULES_DIR / "barrnap.nf").read_text(encoding="utf-8")

        self.assertIn("awk 'NF { value=\\$0 } END { if (value) print value }'", module_text)

    def test_busco_uses_the_last_non_blank_version_line(self) -> None:
        """Require the BUSCO probe to avoid its leading program-name line."""
        module_text = (MODULES_DIR / "busco.nf").read_text(encoding="utf-8")

        self.assertIn("awk 'NF { value=\\$0 } END { if (value) print value }'", module_text)
        self.assertIn('busco_version="\\${busco_version:-NA}"', module_text)

    def test_prokka_uses_the_last_non_blank_version_line(self) -> None:
        """Require the Prokka probe to avoid its leading program-name line."""
        module_text = (MODULES_DIR / "prokka.nf").read_text(encoding="utf-8")

        self.assertIn("awk 'NF { value=\\$0 } END { if (value) print value }'", module_text)
        self.assertIn('prokka_version="\\${prokka_version:-NA}"', module_text)

    def test_padloc_uses_the_version_flag_for_provenance(self) -> None:
        """Require PADLOC provenance to capture the version string, not the banner."""
        module_text = (MODULES_DIR / "padloc.nf").read_text(encoding="utf-8")

        self.assertIn('padloc_version="\\$(padloc --version 2>&1 | awk \'NF { print; exit }\' || echo NA)"', module_text)
        self.assertNotIn('padloc --help 2>&1 | awk \'NF { print; exit }\'', module_text)

    def test_fastani_materialises_inputs_and_fails_early_on_empty_matrix(self) -> None:
        """Require FastANI to copy-stage inputs and fail before cluster parsing on broken output."""
        module_text = (MODULES_DIR / "fastani.nf").read_text(encoding="utf-8")

        self.assertIn("stageInMode 'copy'", module_text)
        self.assertIn("grep -q 'Could not open ' fastani.log", module_text)
        self.assertIn('if [[ -s "${fastani_paths}" && ! -s fastani.matrix ]]; then', module_text)
        self.assertIn('FastANI did not produce a matrix for a non-empty input list.', module_text)

    def test_ccfinder_extracts_a_numeric_version_from_verbose_output(self) -> None:
        """Require CRISPRCasFinder provenance to store a clean version token."""
        module_text = (MODULES_DIR / "ccfinder.nf").read_text(encoding="utf-8")

        self.assertIn(
            "sed -n 's/.*version \\\\([^,[:space:]]*\\\\).*/\\\\1/p' | head -n 1",
            module_text,
        )
        self.assertIn('ccfinder_version="\\${ccfinder_version:-NA}"', module_text)

    def test_ccfinder_does_not_pass_removed_casfinder_path_flags(self) -> None:
        """Require the pinned CRISPRCasFinder invocation to omit dead path flags."""
        module_path = MODULES_DIR / "ccfinder.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertNotIn('-cf "\\${ccfinder_root}/CasFinder-2.0.3" \\', module_text)
        self.assertNotIn('-CASFinder "\\${ccfinder_root}/CasFinder-2.0.3" \\', module_text)

    def test_ccfinder_relies_on_image_provided_tool_compatibility(self) -> None:
        """Require the module to rely on the container's bundled tool compatibility."""
        module_path = MODULES_DIR / "ccfinder.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertNotIn('real_muscle="\\$(command -v muscle || true)"', module_text)
        self.assertNotIn('real_mv="\\$(command -v mv || true)"', module_text)
        self.assertNotIn('tool_bin="\\${task_root}/tool_bin"', module_text)
        self.assertNotIn("translated_args+=(-in", module_text)
        self.assertNotIn("translated_args+=(-out", module_text)
        self.assertNotIn("protected=(%q %q %q)", module_text)
        self.assertNotIn('export PATH="\\${tool_bin}:\\$PATH"', module_text)

    def test_ccfinder_uses_isolated_run_root_outside_tool_output_tree(self) -> None:
        """Require the CRISPRCasFinder run root to stay separate from its output tree."""
        module_text = (MODULES_DIR / "ccfinder.nf").read_text(encoding="utf-8")

        self.assertIn('task_root="\\$PWD"', module_text)
        self.assertIn('run_root="\\${task_root}/ccfinder_run"', module_text)
        self.assertIn('tool_output_root="\\${task_root}/ccfinder_raw"', module_text)
        self.assertIn('pushd "\\${run_root}" >/dev/null', module_text)
        self.assertIn('-outdir "\\${tool_output_root}" \\', module_text)
        self.assertNotIn('-outdir "\\$PWD/ccfinder" \\', module_text)

    def test_ccfinder_stages_local_fasta_and_collects_run_root_outputs(self) -> None:
        """Require the isolated run to reopen the staged FASTA by local basename."""
        module_text = (MODULES_DIR / "ccfinder.nf").read_text(encoding="utf-8")

        self.assertIn('genome_path="\\$(cd "\\$(dirname "${genome}")" && pwd)/\\$(basename "${genome}")"', module_text)
        self.assertIn('genome_name="\\$(basename "${genome}")"', module_text)
        self.assertIn("awk -v output_root=\"\\${task_root}\" '", module_text)
        self.assertIn('output_path = output_root "/" contig_id ".fna"', module_text)
        self.assertIn('sub(/\\\\.[0-9]+\\$/, "", contig_id)', module_text)
        self.assertIn('print \\$0 > output_path', module_text)
        self.assertIn('while (length(sequence_line) > 60) {', module_text)
        self.assertIn('print substr(sequence_line, 1, 60) >> output_path', module_text)
        self.assertIn('cp "\\${genome_path}" "\\${genome_name}"', module_text)
        self.assertIn(': > index.html', module_text)
        self.assertIn('perl "\\${ccfinder_root}/CRISPRCasFinder.pl" -in "\\${genome_name}" \\', module_text)
        self.assertIn(
            "result_json_path=\\$(find \"\\${tool_output_root}\" \"\\${run_root}\" -type f -name 'result.json' | head -n 1 || true)",
            module_text,
        )
        self.assertIn(
            "perl -0pi -e 's/:\\\\s*(?=,|\\\\}|\\\\])/: null/g' \"\\${result_json_path}\"",
            module_text,
        )
        self.assertIn(
            "internal_log_path=\\$(find \"\\${tool_output_root}\" \"\\${run_root}\" -type f \\\\( -name 'ccfinder.log' -o -name 'logFile_*' \\\\) | head -n 1 || true)",
            module_text,
        )

    def test_summarise_ccfinder_consumes_strict_json_only(self) -> None:
        """Require CRISPR summarisation to depend only on the emitted JSON artefact."""
        module_text = (MODULES_DIR / "summarise_ccfinder.nf").read_text(encoding="utf-8")
        workflow_text = (
            ROOT / "subworkflows" / "local" / "per_sample_annotation.nf"
        ).read_text(encoding="utf-8")

        self.assertIn("tuple val(meta), path(result_json)", module_text)
        self.assertNotIn("--ccfinder-dir", module_text)
        self.assertIn("SUMMARISE_CCFINDER(CCFINDER.out.result_json)", workflow_text)
        self.assertNotIn("CCFINDER.out.results.map", workflow_text)

    def test_cohort_ani_combines_busco_summaries_per_lineage(self) -> None:
        """Require the ANI workflow to aggregate BUSCO rows into lineage tables."""
        workflow_path = ROOT / "subworkflows" / "local" / "cohort_ani.nf"
        workflow_text = workflow_path.read_text(encoding="utf-8")

        self.assertIn('tuple("busco_summary_${lineage}.tsv", summary)', workflow_text)
        self.assertIn(".collectFile(", workflow_text)
        self.assertIn("parsed_busco = busco_tables", workflow_text)

    def test_cohort_ani_preserves_staged_manifest_header_order(self) -> None:
        """Require the staged ANI manifest header to stay first in the file."""
        workflow_path = ROOT / "subworkflows" / "local" / "cohort_ani.nf"
        workflow_text = workflow_path.read_text(encoding="utf-8")

        self.assertIn("collectFile(name: 'staged_genomes.tsv', newLine: true, sort: false)", workflow_text)

    def test_final_outputs_consumes_combined_busco_tables(self) -> None:
        """Require final outputs to consume combined BUSCO lineage tables directly."""
        workflow_path = ROOT / "subworkflows" / "local" / "final_outputs.nf"
        workflow_text = workflow_path.read_text(encoding="utf-8")

        self.assertIn("collected_busco = busco_summaries", workflow_text)
        self.assertNotIn(".map { meta, lineage, summary -> summary }", workflow_text)

    def test_final_outputs_keep_manifest_headers_first(self) -> None:
        """Require collected annotation manifests to preserve their header rows."""
        workflow_path = ROOT / "subworkflows" / "local" / "final_outputs.nf"
        workflow_text = workflow_path.read_text(encoding="utf-8")

        self.assertIn("collectFile(name: 'prokka_manifest.tsv', newLine: true, sort: false)", workflow_text)
        self.assertIn("collectFile(name: 'padloc_manifest.tsv', newLine: true, sort: false)", workflow_text)
        self.assertIn("collectFile(name: 'eggnog_manifest.tsv', newLine: true, sort: false)", workflow_text)
        self.assertIn("accession\\tstatus\\twarnings\\texit_code\\tannotations_size\\tresult_file_count", workflow_text)

    def test_final_outputs_calls_manifest_helpers_as_closures(self) -> None:
        """Require closure helpers in final outputs to use explicit `.call(...)`."""
        workflow_path = ROOT / "subworkflows" / "local" / "final_outputs.nf"
        workflow_text = workflow_path.read_text(encoding="utf-8")

        self.assertIn("extractExitCode.call(log)", workflow_text)
        self.assertIn("countTopLevelFiles.call(padlocDir)", workflow_text)
        self.assertIn("countTopLevelFiles.call(eggnogDir)", workflow_text)
        self.assertNotIn("extractExitCode(log)", workflow_text)
        self.assertNotIn("countTopLevelFiles(padlocDir)", workflow_text)
        self.assertNotIn("countTopLevelFiles(eggnogDir)", workflow_text)

    def test_padloc_creates_output_directory_before_running_tool(self) -> None:
        """Require PADLOC output directory creation before invoking the tool."""
        module_text = (MODULES_DIR / "padloc.nf").read_text(encoding="utf-8")

        mkdir_index = module_text.index("mkdir -p padloc")
        run_index = module_text.index('padloc --faa padloc_input.faa')
        self.assertLess(mkdir_index, run_index)

    def test_padloc_uses_the_bundled_database_image(self) -> None:
        """Require PADLOC to rely on the bundled database in its fixed image."""
        module_text = (MODULES_DIR / "padloc.nf").read_text(encoding="utf-8")
        workflow_text = (
            ROOT / "subworkflows" / "local" / "per_sample_annotation.nf"
        ).read_text(encoding="utf-8")
        main_text = (ROOT / "main.nf").read_text(encoding="utf-8")
        collect_versions_text = (MODULES_DIR / "collect_versions.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("tuple val(meta), path(gff), path(faa)", module_text)
        self.assertNotIn("path(padloc_db)", module_text)
        self.assertNotIn("padloc_db_dir=", module_text)
        self.assertNotIn("--data", module_text)
        self.assertIn("PADLOC(PROKKA.out.padloc_inputs)", workflow_text)
        self.assertNotIn("padlocDb = params.padloc_db", main_text)
        self.assertNotIn("--padloc-db", collect_versions_text)

    def test_merge_runtime_database_reports_uses_path_python_lookup(self) -> None:
        """Require merge helper execution to use PATH-resolved Python and helper binaries."""
        module_text = (MODULES_DIR / "merge_runtime_database_reports.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn('script_path="\\$(command -v merge_runtime_database_reports.py)"', module_text)
        self.assertIn('python_path="\\$(command -v python3)"', module_text)
        self.assertNotIn('script_path="/usr/local/bin/merge_runtime_database_reports.py"', module_text)
        self.assertNotIn('/usr/local/bin/python3', module_text)

    def test_eggnog_short_circuit_filters_only_eggnog_inputs(self) -> None:
        """Require acceptance eggNOG short-circuiting to filter only eggNOG jobs."""
        workflow_text = (
            ROOT / "subworkflows" / "local" / "per_sample_annotation.nf"
        ).read_text(encoding="utf-8")
        final_outputs_text = (
            ROOT / "subworkflows" / "local" / "final_outputs.nf"
        ).read_text(encoding="utf-8")
        config_text = (ROOT / "nextflow.config").read_text(encoding="utf-8")

        self.assertIn("eggnog_only_accessions = null", config_text)
        self.assertIn("params.eggnog_only_accessions", workflow_text)
        self.assertIn("EGGNOG(eggnog_inputs.combine(eggnog_db))", workflow_text)
        self.assertNotIn("PROKKA(annotation_candidates.filter", workflow_text)
        self.assertNotIn("CCFINDER(annotation_candidates.filter", workflow_text)
        self.assertNotIn("PADLOC(PROKKA.out.padloc_inputs.filter", workflow_text)
        self.assertIn("configuredEggnogOnlyAccessions = parseConfiguredAccessions.call(params.eggnog_only_accessions)", final_outputs_text)
        self.assertIn("eggnog_short_circuit", final_outputs_text)

    def test_collect_versions_stages_unique_input_names(self) -> None:
        """Require unique staged names for collected version files."""
        module_path = MODULES_DIR / "collect_versions.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn("path version_files, name: 'version_files/versions??.yml'", module_text)

    def test_runtime_database_prep_modules_escape_shell_command_substitutions(self) -> None:
        """Require parser-safe escaping in the new runtime DB prep modules."""
        matches: list[str] = []
        for path in (
            MODULES_DIR / "prepare_runtime_database.nf",
            MODULES_DIR / "download_checkm2_database.nf",
            MODULES_DIR / "download_busco_databases.nf",
            MODULES_DIR / "download_eggnog_database.nf",
            MODULES_DIR / "finalise_runtime_database.nf",
            MODULES_DIR / "merge_runtime_database_reports.nf",
        ):
            matches.extend(find_raw_command_substitutions(path))

        self.assertEqual(matches, [])

    def test_prepare_databases_requires_a_database_destination_and_busco_lineages(self) -> None:
        """Require the separate prep entry point to validate its core params."""
        workflow_text = (ROOT / "prepare_databases.nf").read_text(encoding="utf-8")

        self.assertIn("At least one database destination must be set for prepare_databases.nf.", workflow_text)
        self.assertIn("params.busco_lineages must be a non-empty list.", workflow_text)
        self.assertNotIn("nextflow.enable.configProcessNamesValidation = false", workflow_text)

    def test_prepare_databases_registers_main_subworkflows_for_shared_config_selectors(self) -> None:
        """Require the prep entry point to register main workflow process names."""
        workflow_text = (ROOT / "prepare_databases.nf").read_text(encoding="utf-8")

        self.assertIn(
            "include { BUSCO_DATASET_PREP as UNUSED_BUSCO_DATASET_PREP }",
            workflow_text,
        )
        self.assertIn("include { COHORT_16S as UNUSED_COHORT_16S }", workflow_text)
        self.assertIn("include { COHORT_ANI as UNUSED_COHORT_ANI }", workflow_text)
        self.assertIn(
            "include { COHORT_TAXONOMY as UNUSED_COHORT_TAXONOMY }",
            workflow_text,
        )
        self.assertIn("include { FINAL_OUTPUTS as UNUSED_FINAL_OUTPUTS }", workflow_text)
        self.assertIn(
            "include { INPUT_VALIDATION_AND_STAGING as UNUSED_INPUT_VALIDATION_AND_STAGING }",
            workflow_text,
        )
        self.assertIn(
            "include { PER_SAMPLE_ANNOTATION as UNUSED_PER_SAMPLE_ANNOTATION }",
            workflow_text,
        )
        self.assertIn("include { PER_SAMPLE_QC as UNUSED_PER_SAMPLE_QC }", workflow_text)

    def test_prepare_databases_uses_main_database_params_directly(self) -> None:
        """Require the prep workflow to reuse the main DB params directly."""
        workflow_text = (ROOT / "prepare_databases.nf").read_text(encoding="utf-8")
        merge_module_text = (MODULES_DIR / "merge_runtime_database_reports.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("taxdump   : normaliseDestination.call(params.taxdump)", workflow_text)
        self.assertIn("checkm2   : normaliseDestination.call(params.checkm2_db)", workflow_text)
        self.assertIn("busco_root: normaliseDestination.call(params.busco_db)", workflow_text)
        self.assertIn("codetta   : normaliseDestination.call(params.codetta_db)", workflow_text)
        self.assertIn("eggnog    : normaliseDestination.call(params.eggnog_db)", workflow_text)
        self.assertIn("canonicalFile.absolutePath", workflow_text)
        self.assertIn("def buildMountedDestination = { destination ->", workflow_text)
        self.assertIn("def buildMountedScratchRoot = { destinationParent ->", workflow_text)
        self.assertIn("file(parentFile.absolutePath)", workflow_text)
        self.assertIn("def mountedDestinations = destinations.collectEntries", workflow_text)
        self.assertIn("taxdumpRequest = destinations.taxdump", workflow_text)
        self.assertIn("checkm2Request = destinations.checkm2", workflow_text)
        self.assertIn("buscoRequest = destinations.busco_root", workflow_text)
        self.assertIn("codettaRequest = destinations.codetta", workflow_text)
        self.assertIn("eggnogRequest = destinations.eggnog", workflow_text)
        self.assertNotIn("params.padloc_db", workflow_text)
        self.assertIn('script_path="\\$(command -v merge_runtime_database_reports.py)"', merge_module_text)
        self.assertIn('python_path="\\$(command -v python3)"', merge_module_text)
        self.assertIn('"\\${python_path}" "\\${script_path}" \\', merge_module_text)

    def test_runtime_database_prep_modules_delegate_to_python_helpers(self) -> None:
        """Require the prep workflow to keep validation logic out of Groovy."""
        prepare_module_text = (MODULES_DIR / "prepare_runtime_database.nf").read_text(
            encoding="utf-8"
        )
        busco_module_text = (MODULES_DIR / "download_busco_databases.nf").read_text(
            encoding="utf-8"
        )
        checkm2_module_text = (MODULES_DIR / "download_checkm2_database.nf").read_text(
            encoding="utf-8"
        )
        eggnog_module_text = (MODULES_DIR / "download_eggnog_database.nf").read_text(
            encoding="utf-8"
        )
        finalise_module_text = (MODULES_DIR / "finalise_runtime_database.nf").read_text(
            encoding="utf-8"
        )
        merge_module_text = (MODULES_DIR / "merge_runtime_database_reports.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("stageInMode 'symlink'", prepare_module_text)
        self.assertIn(
            "tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(version), path(scratch_root_path), val(force)",
            prepare_module_text,
        )
        self.assertIn('script_path="\\$(command -v prepare_runtime_databases.py)"', prepare_module_text)
        self.assertIn('python_path="\\$(command -v python3)"', prepare_module_text)
        self.assertIn('destination_path="${destination_parent}/${destination_name}"', prepare_module_text)
        self.assertIn("--taxdump-dest", prepare_module_text)
        self.assertIn("--codetta-dest", prepare_module_text)
        self.assertIn('helper_args+=(--scratch-root "${scratch_root_path}")', prepare_module_text)
        self.assertIn('python_version="\\$("\\${python_path}" --version 2>&1 | sed \'s/^Python //\')"', prepare_module_text)
        self.assertIn('"\\${python_path}" "\\${script_path}" "\\${helper_args[@]}"', prepare_module_text)
        self.assertIn("stageInMode 'symlink'", checkm2_module_text)
        self.assertIn(
            "tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(force)",
            checkm2_module_text,
        )
        self.assertIn('destination_path="${destination_parent}/${destination_name}"', checkm2_module_text)
        self.assertIn('normalise_download_layout() {', checkm2_module_text)
        self.assertIn('nested_root="\\$1/CheckM2_database"', checkm2_module_text)
        self.assertIn("stageInMode 'symlink'", busco_module_text)
        self.assertIn(
            "tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(lineages), val(force)",
            busco_module_text,
        )
        self.assertIn('destination_path="${destination_parent}/${destination_name}"', busco_module_text)
        self.assertIn('busco --download_path "\\${destination_path}" --download "\\${lineage}"', busco_module_text)
        self.assertIn("stageInMode 'symlink'", eggnog_module_text)
        self.assertIn(
            "tuple val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(force)",
            eggnog_module_text,
        )
        self.assertIn('destination_path="${destination_parent}/${destination_name}"', eggnog_module_text)
        self.assertIn('script_path="\\$(command -v download_eggnog_data.py)"', eggnog_module_text)
        self.assertIn('python "\\${script_path}" --data_dir "\\${destination_path}" -y', eggnog_module_text)
        self.assertNotIn('download_eggnog_data.py.patched', eggnog_module_text)
        self.assertNotIn('http://eggnog5.embl.de/download/emapperdb-', eggnog_module_text)
        self.assertIn("stageInMode 'symlink'", finalise_module_text)
        self.assertIn(
            "tuple val(component), val(destination), path(destination_parent), val(destination_name), path(mode_file), val(source_label), path(lineages_file)",
            finalise_module_text,
        )
        self.assertIn('script_path="\\$(command -v finalise_runtime_database.py)"', finalise_module_text)
        self.assertIn('python_path="\\$(command -v python3)"', finalise_module_text)
        self.assertIn('destination_path="${destination_parent}/${destination_name}"', finalise_module_text)
        self.assertIn('python_version="\\$("\\${python_path}" --version 2>&1 | sed \'s/^Python //\')"', finalise_module_text)
        self.assertIn('"\\${python_path}" "\\${script_path}" "\\${helper_args[@]}"', finalise_module_text)
        self.assertIn('script_path="\\$(command -v merge_runtime_database_reports.py)"', merge_module_text)
        self.assertIn('python_path="\\$(command -v python3)"', merge_module_text)
        self.assertIn('"\\${python_path}" "\\${script_path}" \\', merge_module_text)

    def test_runtime_database_prep_respects_configured_busco_lineages(self) -> None:
        """Require the prep workflow to forward params.busco_lineages to BUSCO prep."""
        workflow_text = (ROOT / "prepare_databases.nf").read_text(encoding="utf-8")
        busco_module_text = (MODULES_DIR / "download_busco_databases.nf").read_text(
            encoding="utf-8"
        )
        finalise_module_text = (MODULES_DIR / "finalise_runtime_database.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("params.busco_lineages as List<String>", workflow_text)
        self.assertIn("""printf '%s\\n' "\\${lineages[@]}" > lineages.txt""", busco_module_text)
        self.assertIn('helper_args+=(--busco-lineage "\\${lineage}")', finalise_module_text)

    def test_busco_prep_accepts_valid_roots_without_ready_markers(self) -> None:
        """Require BUSCO prep to finalise valid lineage roots that lack a marker."""
        module_text = (MODULES_DIR / "download_busco_databases.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn('has_all_lineages=true', module_text)
        self.assertIn(
            'if [[ "\\${has_all_lineages}" == "true" && -f "\\${destination_path}/.nf_myco_ready.json" ]]; then',
            module_text,
        )
        self.assertIn('elif [[ "\\${has_all_lineages}" == "true" ]]; then', module_text)
        self.assertIn(
            'echo "BUSCO destination exists but is not ready: \\${destination_path}" >&2',
            module_text,
        )


if __name__ == "__main__":
    unittest.main()
