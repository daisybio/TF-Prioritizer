if (params.rnaseq_samplesheet) { ch_rnaseq_samplesheet = file(params.rnaseq_samplesheet) } else { exit 1, 'RNAseq samplesheet not specified!' }

include { RNASEQ } from '../subworkflows/local/rnaseq'
include { COUNT_PREPROCESSING } from '../modules/local/count_preprocessing'

workflow RNASEQ {
    ch_versions = Channel.empty()

    if (params.rnaseq_counts) {
        ch_count = file(params.rnaseq_counts)
    } else {
        RNASEQ ()
        ch_count = RNASEQ.out.counts
        ch_bigwig = RNASEQ.out.bigwig
        ch_versions = ch_versions.mix(RNASEQ.out.versions)
    }

    COUNT_PREPROCESSING (ch_count, ch_rnaseq_samplesheet)
}