include { CAT_CAT } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE } from '../../modules/nf-core/bedtools/merge/main'
include { GAWK as FILTER_MIN_OCCURRENCE } from '../../modules/nf-core/gawk/main'
include { GAWK as CLEAN_BED } from '../../modules/nf-core/gawk/main'

workflow MERGE_PEAKS {
    take:
        ch_peaks

    main:
        ch_grouped = ch_peaks.map { meta, file -> return [meta.state, meta.antibody, file]}
                .groupTuple(by: [0,1])
                .map { state, antibody, files -> return [[id: state + "_" + antibody, state: state, antibody: antibody], files] }

        CAT_CAT(ch_grouped)
        BEDTOOLS_SORT(CAT_CAT.out.file_out, [])
        BEDTOOLS_MERGE(BEDTOOLS_SORT.out.sorted)
        FILTER_MIN_OCCURRENCE(BEDTOOLS_MERGE.out.bed, [])
        CLEAN_BED(FILTER_MIN_OCCURRENCE.out.output, [])

    emit:
        peaks = CLEAN_BED.out.output
}