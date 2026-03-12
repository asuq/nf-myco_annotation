"""Container contract tests for the PADLOC launcher patch."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IMAGE = "quay.io/biocontainers/padloc:2.0.0--hdfd78af_1"
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
    """Lock the task-local PADLOC launcher bootstrap patch."""

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the PADLOC container launcher contract test.",
    )
    def test_patched_launcher_does_not_use_missing_bundled_data_dir(self) -> None:
        """Require the patched launcher to bypass the missing bundled data directory."""
        shell_script = r"""
set -euo pipefail
mkdir -p /tmp/padloc_bin /tmp/padloc_bootstrap_data /tmp/padloc_dest
cp "$(command -v padloc)" /tmp/padloc_bin/padloc.real
cp /workspace/bin/patch_padloc_launcher.py /tmp/patch_padloc_launcher.py
python3 /tmp/patch_padloc_launcher.py /tmp/padloc_bin/padloc.real
cat <<'EOF' > /tmp/padloc_bin/padloc
#!/usr/bin/env bash
set -euo pipefail
exec bash "$PADLOC_WRAPPER_REAL" "$@"
EOF
chmod +x /tmp/padloc_bin/padloc.real /tmp/padloc_bin/padloc
export PADLOC_WRAPPER_REAL=/tmp/padloc_bin/padloc.real
export PADLOC_BOOTSTRAP_DATA=/tmp/padloc_bootstrap_data
export PATH=/tmp/padloc_bin:$PATH
set +e
padloc --data /tmp/padloc_dest --db-update > /tmp/padloc.out 2>&1
status=$?
set -e
cat /tmp/padloc.out
exit "$status"
"""

        try:
            result = run_command(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--platform",
                    "linux/amd64",
                    "-v",
                    f"{ROOT}:/workspace:ro",
                    IMAGE,
                    "bash",
                    "-lc",
                    shell_script,
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
