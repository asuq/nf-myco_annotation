include { DOWNLOAD_BUSCO_DATASET } from '../../modules/local/download_busco_dataset'

/*
 * Resolve BUSCO lineage datasets for offline per-sample BUSCO runs. Reuse a
 * staged dataset directory by default and only download when explicitly asked.
 */
workflow BUSCO_DATASET_PREP {
    take:
    lineages

    main:
    if (params.prepare_busco_datasets) {
        DOWNLOAD_BUSCO_DATASET(lineages)
        datasets = DOWNLOAD_BUSCO_DATASET.out.dataset
        logs = DOWNLOAD_BUSCO_DATASET.out.log
        versions = DOWNLOAD_BUSCO_DATASET.out.versions
    } else {
        if (!params.busco_download_dir) {
            error "params.busco_download_dir is required unless params.prepare_busco_datasets=true"
        }

        datasets = lineages.map { lineage ->
            tuple(
                lineage,
                file("${params.busco_download_dir}/${lineage}", checkIfExists: true),
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
