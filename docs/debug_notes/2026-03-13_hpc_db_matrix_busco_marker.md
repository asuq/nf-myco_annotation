# HPC BUSCO valid-without-marker mismatch

## Failure

The OIST matrix run failed in `db3_valid_without_marker` for BUSCO even though
the seeded root already contained both configured lineage directories and their
`dataset.cfg` files.

The failing path was:

`db_cases/db3_valid_without_marker/busco`

## Cause

`DOWNLOAD_BUSCO_DATABASES` treated a BUSCO root as ready only when both of
these were true:

- all configured lineage directories were present with `dataset.cfg`
- `.nf_myco_ready.json` already existed

That differed from CheckM2 and eggNOG, which already accept valid existing
directories without markers as `prepared` and let the finalise step write the
marker.

## Fix

Align BUSCO with the intended matrix contract:

- valid lineage root plus marker -> `reuse`
- valid lineage root without marker -> `prepared`
- partial or broken lineage root -> fail unless force rebuild is enabled

## Validation

- `python3 -m unittest tests.test_run_oist_hpc_matrix tests.test_nextflow_module_syntax`
- `python3 -m unittest discover -s tests -p 'test_*.py'`
- `bash bin/run_oist_hpc_matrix.sh --dry-run --hpc-root /tmp/hpc_matrix all`
