"""Regression checks for Nextflow runtime configuration contracts."""

from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_CONFIG = ROOT / "nextflow.config"
DOCKER_CONFIG = ROOT / "conf" / "docker.config"
SINGULARITY_CONFIG = ROOT / "conf" / "singularity.config"
SLURM_CONFIG = ROOT / "conf" / "slurm.config"
OIST_CONFIG = ROOT / "conf" / "oist.config"
LOCAL_CONFIG = ROOT / "conf" / "local.config"
DEBUG_CONFIG = ROOT / "conf" / "debug.config"
BASE_CONFIG = ROOT / "conf" / "base.config"


class NextflowConfigContractsTestCase(unittest.TestCase):
    """Protect configuration needed for containerised local runs."""

    def test_nextflow_defaults_define_runtime_params(self) -> None:
        """Keep warning-prone runtime params defined in the root config."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("ani_score_profile = 'default'", config_text)
        self.assertIn("padloc_db = null", config_text)
        self.assertIn("db_root = null", config_text)
        self.assertIn("download_missing_databases = false", config_text)
        self.assertIn("force_runtime_database_rebuild = false", config_text)
        self.assertIn("runtime_db_link_mode = 'copy'", config_text)
        self.assertIn("runtime_db_scratch_root = null", config_text)
        self.assertIn("taxdump_source = null", config_text)
        self.assertIn("taxdump_version = null", config_text)
        self.assertIn("checkm2_source = null", config_text)
        self.assertIn("checkm2_version = null", config_text)
        self.assertIn("busco_source_root = null", config_text)
        self.assertIn("busco_version = null", config_text)
        self.assertIn("eggnog_source = null", config_text)
        self.assertIn("eggnog_version = null", config_text)
        self.assertIn("padloc_source = null", config_text)
        self.assertIn("padloc_version = null", config_text)
        self.assertIn("eggnog_only_accessions = null", config_text)
        self.assertIn("singularity_cache_dir = null", config_text)
        self.assertIn("singularity_run_options = ''", config_text)
        self.assertIn("includeConfig 'conf/debug.config'", config_text)
        self.assertIn("includeConfig 'conf/oist.config'", config_text)
        self.assertNotIn("use_biocontainers", config_text)
        self.assertNotIn("runtime_db_helper_container", config_text)
        self.assertNotIn("apptainer_cache_dir", config_text)
        self.assertNotIn("apptainer_run_options", config_text)

    def test_singularity_runtime_surface_uses_renamed_config_contract(self) -> None:
        """Keep the renamed Singularity profile, include, and config file wired in."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        singularity_text = SINGULARITY_CONFIG.read_text(encoding="utf-8")
        docker_text = DOCKER_CONFIG.read_text(encoding="utf-8")
        slurm_text = SLURM_CONFIG.read_text(encoding="utf-8")

        self.assertTrue(SINGULARITY_CONFIG.is_file())
        self.assertFalse((ROOT / "conf" / "apptainer.config").exists())
        self.assertIn("includeConfig 'conf/singularity.config'", config_text)
        self.assertNotIn("includeConfig 'conf/apptainer.config'", config_text)
        self.assertIn("singularity {", config_text)
        self.assertNotIn("apptainer {", config_text)
        self.assertIn("profiles {", singularity_text)
        self.assertIn("singularity.enabled = true", singularity_text)
        self.assertIn("singularity.autoMounts = true", singularity_text)
        self.assertIn("params.singularity_cache_dir", singularity_text)
        self.assertIn("params.singularity_run_options", singularity_text)
        self.assertIn("singularity.enabled = false", docker_text)
        self.assertIn('params.slurm_account ? "--account=${params.slurm_account}" : null', slurm_text)
        self.assertNotIn("apptainer.enabled = false", docker_text)

    def test_oist_profile_is_available_as_one_standalone_runtime(self) -> None:
        """Keep the OIST profile self-contained for SLURM plus Singularity runs."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        oist_text = OIST_CONFIG.read_text(encoding="utf-8")

        self.assertTrue(OIST_CONFIG.is_file())
        self.assertIn("includeConfig 'conf/oist.config'", config_text)
        self.assertIn("profiles {", oist_text)
        self.assertIn("oist {", oist_text)
        self.assertIn("process.executor = 'slurm'", oist_text)
        self.assertIn("executor.queueSize = params.executor_queue_size", oist_text)
        self.assertIn("process.queue = params.slurm_queue ?: 'short'", oist_text)
        self.assertIn("params.slurm_cluster_options ?: null", oist_text)
        self.assertIn("singularity.enabled = true", oist_text)
        self.assertIn("singularity.autoMounts = true", oist_text)
        self.assertIn("singularity.cacheDir = params.singularity_cache_dir ?: null", oist_text)
        self.assertIn("docker.enabled = false", oist_text)
        self.assertIn("process.containerOptions = params.singularity_run_options ?: ''", oist_text)
        self.assertIn("process.resourceLimits = [", oist_text)
        self.assertIn("withLabel: process_medium", oist_text)
        self.assertIn("withLabel: process_high", oist_text)
        self.assertNotIn("params.slurm_account", oist_text)
        self.assertNotIn("beforeScript", oist_text)

    def test_debug_profile_sets_default_eggnog_smoke_accession(self) -> None:
        """Provide one composable debug profile for single-sample eggNOG runs."""
        config_text = DEBUG_CONFIG.read_text(encoding="utf-8")

        self.assertIn("profiles {", config_text)
        self.assertIn("debug {", config_text)
        self.assertIn("eggnog_only_accessions = 'GCA_000027325.1'", config_text)

    def test_python_container_uses_shared_repo_owned_helper_image(self) -> None:
        """Use one shared helper image that carries the ANI scientific stack."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("python_container = 'quay.io/asuq1617/python-scipy:3.12'", config_text)
        self.assertNotIn("python_container = 'python:3.12'", config_text)

    def test_ccfinder_container_points_at_the_clean_runtime_tag(self) -> None:
        """Use the cleaned CRISPRCasFinder image tag by default."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("ccfinder_container = 'quay.io/asuq1617/ccfinder:4.2.30'", config_text)
        self.assertNotIn("ccfinder_container = 'quay.io/asuq1617/ccfinder:4.3.2'", config_text)

    def test_python_helper_processes_are_pinned_explicitly(self) -> None:
        """Keep shared helper processes on the single Python helper image."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn("withName: CLUSTER_ANI {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = params.python_container", config_text)
        self.assertIn("withName: SELECT_ANI_REPRESENTATIVES {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = params.python_container", config_text)
        self.assertIn("withName: WRITE_SAMPLE_STATUS {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = params.python_container", config_text)

    def test_seqtk_helper_processes_are_pinned_explicitly(self) -> None:
        """Keep seqtk-backed helper processes on the shared seqtk image."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn("withName: STAGE_INPUTS {\n        container = params.seqtk_container", config_text)
        self.assertIn("withName: CALCULATE_ASSEMBLY_STATS {\n        container = params.seqtk_container", config_text)

    def test_runtime_database_prep_processes_use_the_dedicated_helper_image(self) -> None:
        """Keep the prep entry point on the dedicated runtime DB helper image."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn(
            "withName: PREP_RUNTIME_DATABASE {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = 'quay.io/asuq1617/nf-myco-db:0.1'",
            config_text,
        )
        self.assertIn(
            "withName: PREP_BUSCO_DATABASES {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = 'quay.io/asuq1617/nf-myco-db:0.1'",
            config_text,
        )
        self.assertIn(
            "withName: MERGE_RUNTIME_DATABASE_REPORTS {\n        errorStrategy = 'terminate'\n        maxRetries = 0\n        container = 'quay.io/asuq1617/nf-myco-db:0.1'",
            config_text,
        )

    def test_docker_profile_forces_amd64_ccfinder_on_arm_hosts(self) -> None:
        """Allow Apple Silicon Docker runs to pull the pinned CCFINDER image."""
        config_text = DOCKER_CONFIG.read_text(encoding="utf-8")

        self.assertIn("withName: CCFINDER", config_text)
        self.assertIn("container = params.ccfinder_container", config_text)
        self.assertIn("--platform linux/amd64", config_text)
        self.assertIn("['aarch64', 'arm64']", config_text)

    def test_local_profile_caps_memory_for_workstation_runs(self) -> None:
        """Keep local acceptance runs within a typical workstation memory budget."""
        config_text = LOCAL_CONFIG.read_text(encoding="utf-8")

        self.assertIn("params.max_memory = 16.GB", config_text)
        self.assertIn("withName: EGGNOG", config_text)
        self.assertIn("container = params.eggnog_container", config_text)
        self.assertIn("cpus = 4", config_text)
        self.assertIn("memory = 8.GB", config_text)

    def test_reporting_assets_overwrite_on_resumed_runs(self) -> None:
        """Allow resumed local runs to rewrite pipeline info artefacts safely."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn("timeline {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/timeline.html\"\n    overwrite = true", config_text)
        self.assertIn("report {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/report.html\"\n    overwrite = true", config_text)
        self.assertIn("trace {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/trace.tsv\"\n    overwrite = true", config_text)
        self.assertIn("dag {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/dag.html\"\n    overwrite = true", config_text)


if __name__ == "__main__":
    unittest.main()
