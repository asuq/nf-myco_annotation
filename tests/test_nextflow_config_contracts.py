"""Regression checks for Nextflow runtime configuration contracts."""

from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_CONFIG = ROOT / "nextflow.config"
DOCKER_CONFIG = ROOT / "conf" / "docker.config"
LOCAL_CONFIG = ROOT / "conf" / "local.config"
BASE_CONFIG = ROOT / "conf" / "base.config"


class NextflowConfigContractsTestCase(unittest.TestCase):
    """Protect configuration needed for containerised local runs."""

    def test_nextflow_defaults_define_runtime_params(self) -> None:
        """Keep warning-prone runtime params defined in the root config."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("ani_score_profile = 'default'", config_text)
        self.assertIn("use_biocontainers = true", config_text)

    def test_python_container_uses_shared_repo_owned_helper_image(self) -> None:
        """Use one shared helper image that carries the ANI scientific stack."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("python_container = 'quay.io/asuq1617/python-scipy:3.12'", config_text)
        self.assertNotIn("python_container = 'python:3.12'", config_text)

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
