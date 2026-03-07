# Test Data

`stub/` contains the minimal fixture set used by `-profile test -stub-run`.

- `sample_sheet.csv`: one-sample manifest with a local FASTA path
- `metadata.tsv`: matching metadata row with the ANI-required fields
- `taxdump/`: pinned miniature `names.dmp` and `nodes.dmp`
- `busco/`: empty lineage directories that satisfy offline path checks
- `checkm2_db/` and `eggnog_db/`: placeholder resource roots for config wiring

`headers/` contains header-only TSV templates used to build stable combined
tables even when a gated branch yields zero sample rows.
