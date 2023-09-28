include { CLEAN_BED } from "../../modules/local/clean_bed"
include { BEDTOOLS_SORT as SORT_PEAKS } from "../../modules/nf-core/bedtools/sort/main"
include { BEDTOOLS_SORT as SORT_FOOTPRINTS } from "../../modules/nf-core/bedtools/sort/main"
include { MERGE_PEAKS } from "./merge_peaks"
include { BEDTOOLS_CLOSEST } from "../../modules/nf-core/bedtools/closest/main"
include { GAWK as FILTER_CLOSEST } from "../../modules/nf-core/gawk/main"
include { CAT_CAT } from "../../modules/nf-core/cat/cat/main"
include { BEDTOOLS_MERGE } from "../../modules/nf-core/bedtools/merge/main"
include { BEDTOOLS_SUBTRACT } from "../../modules/nf-core/bedtools/subtract/main"
include { TEPIC } from "../../modules/local/tepic"

workflow PEAKS {
    take:
        ch_peaks
    
    main:
        CLEAN_BED(ch_peaks)

        ch_sorted = params.merge_peaks ? 
                        MERGE_PEAKS(CLEAN_BED.out)    .peaks : 
                        SORT_PEAKS (CLEAN_BED.out, []).out.sorted

        if (params.peak_search_type == "inside") {
            ch_footprints = ch_sorted
        } else {
            BEDTOOLS_CLOSEST(ch_sorted.map{ meta, bed_file -> 
                copy_name = ".temp/" + bed_file.getName() + ".copy"
                bed_file.mklink(copy_name, overwrite: true)
                copy_file = file(copy_name)
                return [meta, bed_file, copy_file]
            }, [])
        
            FILTER_CLOSEST(BEDTOOLS_CLOSEST.out.output, [])

            if (params.peak_search_type == "incl_between") {
                ch_joined = FILTER_CLOSEST.out.output
                    .map{ meta, bed_file -> [meta.state, meta.antibody, bed_file]}
                    .join(
                        ch_sorted.map{ meta, bed_file -> [meta.state, meta.antibody, bed_file]},
                        by: [0, 1]
                    ).map{
                        state, antibody, closest, bed_file -> 
                            [[id: state + "_" + antibody, state: state, antibody: antibody], [bed_file, closest]]
                    }

                CAT_CAT(ch_joined)
                SORT_FOOTPRINTS(CAT_CAT.out.file_out, [])
                BEDTOOLS_MERGE(SORT_FOOTPRINTS.out.sorted)
                ch_footprints = BEDTOOLS_MERGE.out.bed
            } else {
                ch_footprints = FILTER_CLOSEST.out.output
            }
        }

   

        if (params.blacklist) {
            BEDTOOLS_SUBTRACT(ch_footprints.map{ meta, bed_file -> [meta, bed_file, params.blacklist]})
            ch_blacklisted = BEDTOOLS_SUBTRACT.out.bed
        } else {
            ch_blacklisted = ch_footprints
        }

        /*
        TEPIC (
            ch_peaks,
            params.pwm,
            params.gtf,
            params.fasta,
            params.tepic_window_size,
            params.tepic_loop_windows,
            params.tepic_exponential_decay,
            params.tepic_normalize_peak_length,
            params.tepic_max_minutes_per_chr,
            params.tepic_original_scaling,
            params.tepic_p_value
        )
        */
}