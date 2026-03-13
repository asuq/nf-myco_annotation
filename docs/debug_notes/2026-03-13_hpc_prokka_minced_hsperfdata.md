## Failure

- HPC tracked real run reported one failed sample: `GCA_000027325.1`.
- The published sample outputs showed:
  - empty `prokka.faa`
  - empty `prokka.gff`
  - downstream `padloc` failure on an empty FAA
  - downstream `eggnog` failure on an empty FAA

## Root cause

- Prokka aborted during its `minced` version check.
- The first line returned to Prokka was a Java `hsperfdata` warning rather than the
  expected `minced 0.4.2` banner.
- Prokka could not parse the version string and exited with `exit_code=2`.

## Fix

- Stop using the upstream Prokka BioContainer directly.
- Add a fixed custom Prokka image that wraps the bioinformatics `minced` launcher.
- The wrapper injects `-XX:+PerfDisableSharedMem` directly into the `java`
  invocation before executing the real `minced.jar`.
- Switch the default `prokka_container` to the fixed image.

## Expected impact

- `minced --version` has a deterministic first line on HPC.
- Prokka no longer aborts during tool discovery for this sample.
- PADLOC and eggNOG receive a non-empty FAA again without needing their own changes.
