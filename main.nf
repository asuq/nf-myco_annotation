#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PER_SAMPLE_ANNOTATION } from './subworkflows/local/per_sample_annotation'
include { PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'
include { COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI } from './subworkflows/local/cohort_ani'

workflow {
    log.warn 'PADLOC outputs are produced but not merged into the master table in v1.'
    log.warn 'This repository currently contains the pipeline skeleton only.'
}
