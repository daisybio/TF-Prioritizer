if (params.rnaseq_samplesheet) { ch_rnaseq_samplesheet = Channel.value(file(params.rnaseq_samplesheet)) } else { exit 1, 'RNAseq samplesheet not specified!' }

include { RNASEQ } from './nf-core_rnaseq'
include { COUNT_NORMALIZATION } from '../../modules/local/counts/count_preprocessing'
include { FILTER_COUNTS } from '../../modules/local/counts/count_preprocessing'
include { GROUP_COUNTS } from '../../modules/local/counts/count_preprocessing'
include { DESEQ2 } from '../../modules/local/counts/deseq'
include { ENSG_MAP_CREATION } from '../../modules/local/counts/ensg_mapping'
include { ENSG_MAPPING } from '../../modules/local/counts/ensg_mapping'

workflow RNASEQ {
    ch_versions = Channel.empty()

    if (params.rnaseq_counts) {
        ch_count = file(params.rnaseq_counts)
    } else {
        RNASEQ ()
        ch_count = RNASEQ.out.counts
        ch_bigwig = RNASEQ.out.bigwig
        ch_versions = ch_versions.mix(RNASEQ.out.versions)
        ch_tpm = RNASEQ.out.tpm
    }

    ch_map = ENSG_MAP_CREATION (ch_count, Channel.value(params.taxonomy))

    FILTER_COUNTS (ch_count, params.min_count, params.min_tpm)
    COUNT_NORMALIZATION (FILTER_COUNTS.out.count, ch_rnaseq_samplesheet)
    ENSG_MAPPING (COUNT_NORMALIZATION.out, ch_map)

    ch_groups = ch_rnaseq_samplesheet
        .splitCsv(header: true, sep: ',')
        .map { it.group }
        .unique()
    
    ch_pairings = ch_groups.combine(ch_groups)
        .filter { it[0] < it[1] }

    DESEQ2 (ENSG_MAPPING.out, ch_rnaseq_samplesheet, ch_pairings)
    
    emit:
    deseq2 = DESEQ2.out
    versions = ch_versions
    ensg_map = ch_map
}