"""Regression tests for parser-safe Nextflow module scripts."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULES_DIR = ROOT / "modules" / "local"
RAW_COMMAND_SUBSTITUTION = re.compile(r'(?<!\\)\$\(')


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
        """Require CheckM2 to accept a database directory containing one `.dmnd` file."""
        module_path = MODULES_DIR / "checkm2.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('database_path="${checkm2Db}"', module_text)
        self.assertIn('dmnd_candidates=("\\${database_path}"/*.dmnd)', module_text)
        self.assertIn('--database_path "\\${database_path}"', module_text)

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
        """Require BUSCO offline runs to stage lineage datasets in BUSCO's expected layout."""
        module_path = MODULES_DIR / "busco.nf"
        module_text = module_path.read_text(encoding="utf-8")

        self.assertIn('busco_download_root="busco_downloads"', module_text)
        self.assertIn('staged_lineage_dir="\\${busco_download_root}/lineages/${lineage}"', module_text)
        self.assertIn('dataset_source="\\$(cd "${dataset_dir}" && pwd)"', module_text)
        self.assertIn('ln -s "\\${dataset_source}" "\\${staged_lineage_dir}"', module_text)
        self.assertIn('--download_path "\\${busco_download_root}"', module_text)
        self.assertIn('--lineage_dataset "${lineage}"', module_text)

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

    def test_padloc_requires_staged_database_directory(self) -> None:
        """Require PADLOC to use a staged database directory rather than image defaults."""
        module_text = (MODULES_DIR / "padloc.nf").read_text(encoding="utf-8")
        workflow_text = (
            ROOT / "subworkflows" / "local" / "per_sample_annotation.nf"
        ).read_text(encoding="utf-8")
        main_text = (ROOT / "main.nf").read_text(encoding="utf-8")
        collect_versions_text = (MODULES_DIR / "collect_versions.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("tuple val(meta), path(gff), path(faa), path(padloc_db)", module_text)
        self.assertIn('padloc_db_dir="\\$(cd "${padloc_db}" && pwd)"', module_text)
        self.assertIn('if [[ ! -f "\\${padloc_db_dir}/hmm/padlocdb.hmm" ]]; then', module_text)
        self.assertIn('--data "\\${padloc_db_dir}"', module_text)
        self.assertIn("PADLOC(PROKKA.out.padloc_inputs.combine(padloc_db))", workflow_text)
        self.assertIn("padlocDb = params.padloc_db", main_text)
        self.assertIn('--padloc-db "${params.padloc_db ?: \'NA\'}" \\', collect_versions_text)

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
        self.assertIn("EGGNOG(eggnog_inputs)", workflow_text)
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
            MODULES_DIR / "download_padloc_database.nf",
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

    def test_prepare_databases_uses_main_database_params_directly(self) -> None:
        """Require the prep workflow to reuse the main DB params directly."""
        workflow_text = (ROOT / "prepare_databases.nf").read_text(encoding="utf-8")
        merge_module_text = (MODULES_DIR / "merge_runtime_database_reports.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("taxdump   : normaliseDestination.call(params.taxdump)", workflow_text)
        self.assertIn("checkm2   : normaliseDestination.call(params.checkm2_db)", workflow_text)
        self.assertIn("busco_root: normaliseDestination.call(params.busco_db)", workflow_text)
        self.assertIn("eggnog    : normaliseDestination.call(params.eggnog_db)", workflow_text)
        self.assertIn("padloc    : normaliseDestination.call(params.padloc_db)", workflow_text)
        self.assertIn("taxdumpRequest = destinations.taxdump", workflow_text)
        self.assertIn("checkm2Request = destinations.checkm2", workflow_text)
        self.assertIn("buscoRequest = destinations.busco_root", workflow_text)
        self.assertIn("eggnogRequest = destinations.eggnog", workflow_text)
        self.assertIn("padlocRequest = destinations.padloc", workflow_text)
        self.assertIn('"\\${script_path}" \\', merge_module_text)

    def test_runtime_database_prep_modules_delegate_to_python_helpers(self) -> None:
        """Require the prep workflow to keep validation logic out of Groovy."""
        prepare_module_text = (MODULES_DIR / "prepare_runtime_database.nf").read_text(
            encoding="utf-8"
        )
        busco_module_text = (MODULES_DIR / "download_busco_databases.nf").read_text(
            encoding="utf-8"
        )
        finalise_module_text = (MODULES_DIR / "finalise_runtime_database.nf").read_text(
            encoding="utf-8"
        )
        merge_module_text = (MODULES_DIR / "merge_runtime_database_reports.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn('script_path="/usr/local/bin/prepare_runtime_databases.py"', prepare_module_text)
        self.assertIn("--taxdump-dest", prepare_module_text)
        self.assertIn('"\\${script_path}" "\\${helper_args[@]}"', prepare_module_text)
        self.assertIn('busco --download_path "\\${destination_path}" --download "\\${lineage}"', busco_module_text)
        self.assertIn('script_path="/usr/local/bin/finalise_runtime_database.py"', finalise_module_text)
        self.assertIn('"\\${script_path}" "\\${helper_args[@]}"', finalise_module_text)
        self.assertIn('script_path="/usr/local/bin/merge_runtime_database_reports.py"', merge_module_text)
        self.assertIn('"\\${script_path}" \\', merge_module_text)

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


if __name__ == "__main__":
    unittest.main()
