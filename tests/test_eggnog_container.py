"""Container contract tests for the custom eggNOG runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "eggnog" / "Dockerfile"
IMAGE = "codex-eggnog:test"
RUN_DOCKER_TESTS = os.environ.get("RUN_DOCKER_TESTS") == "1"


def run_command(args: list[str], *, timeout: int | None = None) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and return the completed process."""
    return subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
        timeout=timeout,
    )


class EggnogContainerContractTestCase(unittest.TestCase):
    """Lock the custom eggNOG image contract."""

    def test_dockerfile_patches_the_downloader_at_build_time(self) -> None:
        """Require the custom image to fix the stale download hosts."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn(
            "FROM quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2",
            dockerfile_text,
        )
        self.assertIn('script_path="$(command -v download_eggnog_data.py)"', dockerfile_text)
        self.assertIn("http://eggnogdb.embl.de/download/emapperdb-", dockerfile_text)
        self.assertIn("http://eggnog5.embl.de/download/emapperdb-", dockerfile_text)
        self.assertIn("http://eggnogdb.embl.de/download/novel_fams-", dockerfile_text)
        self.assertIn("http://eggnog5.embl.de/download/novel_fams-", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the eggNOG container image contract test.",
    )
    def test_custom_image_installs_the_fixed_downloader(self) -> None:
        """Build the custom image and require the downloader URLs to be fixed."""
        build_result = run_command(
            [
                "docker",
                "build",
                "--platform",
                "linux/amd64",
                "-f",
                str(DOCKERFILE),
                "-t",
                IMAGE,
                ".",
            ]
        )
        self.assertEqual(
            build_result.returncode,
            0,
            msg=f"Docker build failed.\nSTDOUT:\n{build_result.stdout}\nSTDERR:\n{build_result.stderr}",
        )

        result = run_command(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                IMAGE,
                "bash",
                "-lc",
                (
                    "set -euo pipefail; "
                    'SCRIPT_PATH="$(command -v download_eggnog_data.py)"; '
                    "test -f \"$SCRIPT_PATH\"; "
                    "! grep -F 'http://eggnogdb.embl.de/download/emapperdb-' \"$SCRIPT_PATH\"; "
                    "! grep -F 'http://eggnogdb.embl.de/download/novel_fams-' \"$SCRIPT_PATH\"; "
                    "grep -F 'http://eggnog5.embl.de/download/emapperdb-' \"$SCRIPT_PATH\" >/dev/null; "
                    "grep -F 'http://eggnog5.embl.de/download/novel_fams-' \"$SCRIPT_PATH\" >/dev/null; "
                    "python -c \"import importlib.metadata as m; print(m.version('eggnog-mapper'))\""
                ),
            ],
            timeout=30,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"Custom eggNOG image contract failed.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )
        self.assertIn("2.1.13", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
