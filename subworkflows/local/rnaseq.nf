if (params.rnaseq_samplesheet) { ch_rnaseq_samplesheet = Channel.value(file(params.rnaseq_samplesheet)) } else { exit 1, 'RNAseq samplesheet not specified!' }

include { COUNT_NORMALIZATION } from '../../modules/local/counts/count_preprocessing'
include { FILTER_COUNTS } from '../../modules/local/counts/count_preprocessing'
include { DESEQ2 } from '../../modules/local/counts/deseq'
include { ENSG_MAP_CREATION } from '../../modules/local/counts/ensg_mapping'
include { ENSG_MAPPING as COUNTS_ENSG_MAPPING } from '../../modules/local/counts/ensg_mapping'
include { ENSG_MAPPING as TPM_ENSG_MAPPING } from '../../modules/local/counts/ensg_mapping'
include { GROUP_COUNTS as COUNT_GROUP } from '../../modules/local/counts/count_preprocessing'
include { GROUP_COUNTS as TPM_GROUP } from '../../modules/local/counts/count_preprocessing'
include { DESEQ2_GROUP } from '../../modules/local/counts/deseq'

workflow RNASEQ {
    ch_versions = Channel.empty()

    ch_count = file(params.rnaseq_counts)

    ch_map = ENSG_MAP_CREATION (ch_count, Channel.value(file(params.tepic_gtf)))

    FILTER_COUNTS (ch_count, params.min_count, params.min_tpm)
    COUNT_NORMALIZATION (FILTER_COUNTS.out.count, ch_rnaseq_samplesheet)
    COUNTS_ENSG_MAPPING (COUNT_NORMALIZATION.out, ch_map)
    TPM_ENSG_MAPPING (FILTER_COUNTS.out.tpm, ch_map)

    ch_groups = ch_rnaseq_samplesheet
        .splitCsv(header: true, sep: ',')
        .map { it.group }
        .unique()
    
    ch_pairings = ch_groups.combine(ch_groups)
        .filter { it[0] < it[1] }

    DESEQ2 (COUNTS_ENSG_MAPPING.out, ch_rnaseq_samplesheet, ch_pairings) // group1, group2, deseq-file

    COUNT_GROUP (COUNTS_ENSG_MAPPING.out, ch_rnaseq_samplesheet)
    TPM_GROUP (TPM_ENSG_MAPPING.out, ch_rnaseq_samplesheet)
    DESEQ2_GROUP ( DESEQ2.out
                       .map { it[2] }
                       .collect()
                 )
    
    emit:
    deseq2 = DESEQ2.out
    versions = ch_versions
    ensg_map = ch_map
    deseq2_grouped = DESEQ2_GROUP.out
    count = COUNT_GROUP.out
    count_per_sample = COUNTS_ENSG_MAPPING.out
    tpm = TPM_GROUP.out
    samplesheet = ch_rnaseq_samplesheet
}