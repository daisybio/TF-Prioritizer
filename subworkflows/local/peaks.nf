include { CLEAN_BED } from "../../modules/local/clean_bed"
include { BEDTOOLS_SORT } from "../../modules/nf-core/bedtools/sort/main"
include { MERGE_PEAKS } from "./merge_peaks"

workflow PEAKS {
    take:
        ch_peaks
    
    main:
        CLEAN_BED(ch_peaks)

        if (params.merge_peaks) {
            MERGE_PEAKS(CLEAN_BED.out)
        } else {
            BEDTOOLS_SORT(CLEAN_BED.out, []).out.sorted
        }
}