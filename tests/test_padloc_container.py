"""Container contract tests for the custom PADLOC runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "padloc" / "Dockerfile"
IMAGE = "codex-padloc:test"
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


class PadlocContainerContractTestCase(unittest.TestCase):
    """Lock the custom PADLOC image contract."""

    def test_dockerfile_patches_the_launcher_at_build_time(self) -> None:
        """Require the PADLOC Dockerfile to patch the broken bundled data path."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn("FROM quay.io/biocontainers/padloc:2.0.0--hdfd78af_1", dockerfile_text)
        self.assertIn('mkdir -p "/tmp/padloc-launcher-data"', dockerfile_text)
        self.assertIn('DATA=$(normpath "/tmp/padloc-launcher-data")', dockerfile_text)
        self.assertNotIn("patch_padloc_launcher.py", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the PADLOC container image contract test.",
    )
    def test_custom_image_avoids_the_missing_bundled_data_dir(self) -> None:
        """Require the custom image to avoid the broken bundled data path."""
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

        try:
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
                        "PADLOC_BIN=$(command -v padloc); "
                        "grep -F 'mkdir -p \"/tmp/padloc-launcher-data\"' \"$PADLOC_BIN\" >/dev/null; "
                        "grep -F 'DATA=$(normpath \"/tmp/padloc-launcher-data\")' \"$PADLOC_BIN\" >/dev/null; "
                        "set +e; "
                        "padloc --data /tmp/padloc_dest --db-update > /tmp/padloc.out 2>&1; "
                        "status=$?; "
                        "set -e; "
                        "cat /tmp/padloc.out; "
                        "exit \"$status\""
                    ),
                ],
                timeout=30,
            )
            output = result.stdout + result.stderr
        except subprocess.TimeoutExpired as exc:
            output = (exc.stdout or "") + (exc.stderr or "")

        self.assertNotIn("padloc_bin/../data", output)
        self.assertNotIn('DATA=$(normpath "${SRC_DIR}/../data")', output)


if __name__ == "__main__":
    unittest.main()
