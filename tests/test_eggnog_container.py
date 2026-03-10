"""Container contract tests for the custom eggNOG image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "eggnog" / "Dockerfile"
GUNZIP_WRAPPER = ROOT / "docker" / "eggnog" / "gunzip"
IMAGE_TAG = "codex-eggnog-pigz:test"
RUN_DOCKER_TESTS = os.environ.get("RUN_DOCKER_TESTS") == "1"


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and return the completed process."""
    return subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


class EggnogContainerContractTestCase(unittest.TestCase):
    """Lock the custom eggNOG image contract."""

    def test_dockerfile_installs_pigz_and_overrides_gunzip(self) -> None:
        """Require the image definition to add pigz and the wrapper script."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")
        wrapper_text = GUNZIP_WRAPPER.read_text(encoding="utf-8")

        self.assertIn("FROM quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2", dockerfile_text)
        self.assertIn("apt-get install --yes --no-install-recommends pigz", dockerfile_text)
        self.assertIn("COPY docker/eggnog/gunzip /usr/local/bin/gunzip", dockerfile_text)
        self.assertEqual(wrapper_text, "#!/bin/sh\nexec pigz -d \"$@\"\n")

    @unittest.skipUnless(RUN_DOCKER_TESTS, "Set RUN_DOCKER_TESTS=1 to build the eggNOG image.")
    def test_built_image_exposes_pigz_backed_gunzip(self) -> None:
        """Build the image and require gunzip to dispatch through pigz."""
        build_result = run_command(
            [
                "docker",
                "build",
                "--platform",
                "linux/amd64",
                "-f",
                str(DOCKERFILE),
                "-t",
                IMAGE_TAG,
                ".",
            ]
        )
        self.assertEqual(
            build_result.returncode,
            0,
            msg=f"Docker build failed.\nSTDOUT:\n{build_result.stdout}\nSTDERR:\n{build_result.stderr}",
        )

        run_result = run_command(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                IMAGE_TAG,
                "bash",
                "-lc",
                (
                    "set -euo pipefail; "
                    "command -v pigz >/dev/null; "
                    "test \"$(command -v gunzip)\" = \"/usr/local/bin/gunzip\"; "
                    "gunzip --help >/dev/null"
                ),
            ]
        )
        self.assertEqual(
            run_result.returncode,
            0,
            msg=f"eggNOG image contract failed.\nSTDOUT:\n{run_result.stdout}\nSTDERR:\n{run_result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
