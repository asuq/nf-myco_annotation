include { DOWNLOAD_BUSCO_DATASET } from '../../modules/local/download_busco_dataset'

/*
 * Resolve BUSCO lineage datasets for offline per-sample BUSCO runs. Reuse a
 * staged dataset directory by default and only download when explicitly asked.
 */
workflow BUSCO_DATASET_PREP {
    take:
    lineages

    main:
    def parseBoolean = { value ->
        if (value instanceof Boolean) {
            return value
        }
        if (value == null) {
            return false
        }
        return value.toString().toBoolean()
    }
    def downloadEnabled = parseBoolean.call(params.prepare_busco_datasets)
    def buscoDbRoot = params.busco_db ? new File(params.busco_db.toString()) : null
    def missingDatasetError = { lineage, datasetPath ->
        error "BUSCO lineage dataset is missing: ${lineage}. Expected path: ${datasetPath}. " +
            "Create it with prepare_databases.nf using the same --busco_db and " +
            "--busco_lineages values, or rerun main.nf with --prepare_busco_datasets true."
    }

    if (downloadEnabled) {
        if (params.busco_db) {
            reusableLineages = lineages.filter { lineage ->
                new File(params.busco_db.toString(), lineage.toString()).exists()
            }
            downloadLineages = lineages.filter { lineage ->
                !new File(params.busco_db.toString(), lineage.toString()).exists()
            }
            reusedDatasets = reusableLineages.map { lineage ->
                def datasetPath = new File(params.busco_db.toString(), lineage.toString()).absolutePath
                tuple(
                    lineage,
                    file(datasetPath, checkIfExists: true),
                )
            }
        } else {
            downloadLineages = lineages
            reusedDatasets = Channel.empty()
        }

        DOWNLOAD_BUSCO_DATASET(downloadLineages)
        datasets = reusedDatasets.mix(DOWNLOAD_BUSCO_DATASET.out.dataset)
        logs = DOWNLOAD_BUSCO_DATASET.out.log
        versions = DOWNLOAD_BUSCO_DATASET.out.versions
    } else {
        if (!params.busco_db) {
            error "params.busco_db is required unless params.prepare_busco_datasets=true"
        }

        def stubDatasetFallback = null
        if (workflow.stubRun) {
            stubDatasetFallback = buscoDbRoot
                .listFiles()
                ?.findAll { it.isDirectory() }
                ?.sort { it.name }
                ?.findResult { datasetDir -> datasetDir.absolutePath }
        }

        datasets = lineages.map { lineage ->
            def datasetPath = new File(params.busco_db.toString(), lineage.toString()).absolutePath
            if (
                workflow.stubRun
                && !new File(datasetPath).exists()
                && stubDatasetFallback != null
            ) {
                datasetPath = stubDatasetFallback
            }
            if (!new File(datasetPath).exists()) {
                missingDatasetError.call(lineage, datasetPath)
            }

            tuple(
                lineage,
                file(datasetPath, checkIfExists: true),
            )
        }
        logs = Channel.empty()
        versions = Channel.empty()
    }

    emit:
    datasets = datasets
    logs = logs
    versions = versions
}
