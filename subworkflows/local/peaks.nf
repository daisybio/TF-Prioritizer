include { BEDTOOLS_SORT as SORT_PEAKS } from "../../modules/nf-core/bedtools/sort/main"
include { BEDTOOLS_SORT as SORT_FOOTPRINTS } from "../../modules/nf-core/bedtools/sort/main"
include { MERGE_PEAKS } from "./merge_peaks"
include { BEDTOOLS_CLOSEST } from "../../modules/nf-core/bedtools/closest/main"
include { GAWK as FILTER_CLOSEST } from "../../modules/nf-core/gawk/main"
include { GAWK as CLEAN_BED_ANNOTATE_SAMPLES } from "../../modules/nf-core/gawk/main"
include { CAT_CAT as MERGE_ORIGINAL } from "../../modules/nf-core/cat/cat/main"
include { BEDTOOLS_MERGE } from "../../modules/nf-core/bedtools/merge/main"
include { BEDTOOLS_SUBTRACT as SUBTRACT_BLACKLIST } from "../../modules/nf-core/bedtools/subtract/main"
include { STARE } from "../../modules/local/tepic/stare"
include { REMOVE_VERSION_ANNOTATIONS } from "../../modules/local/tepic/remove_version_annotations"
include { COMBINE_TABLES as AFFINITY_MEAN } from "../../modules/local/combine_tables"
include { COMBINE_TABLES as AFFINITY_SUM } from "../../modules/local/combine_tables"
include { COMBINE_TABLES as AFFINITY_RATIO } from "../../modules/local/combine_tables"


workflow PEAKS {
    take:
        ch_peaks
    
    main:
        CLEAN_BED_ANNOTATE_SAMPLES(ch_peaks, [])

        if (params.merge_peaks) {
            ch_sorted = MERGE_PEAKS(CLEAN_BED_ANNOTATE_SAMPLES.out.output).peaks
        } else {
            SORT_PEAKS (CLEAN_BED_ANNOTATE_SAMPLES.out.output, [])
            ch_sorted = SORT_PEAKS.out.sorted
        }

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
                    .join(ch_sorted)
                    .map {
                        meta, closest, bed_file -> 
                            [meta, [closest, bed_file]]
                    }

                MERGE_ORIGINAL(ch_joined)
                SORT_FOOTPRINTS(MERGE_ORIGINAL.out.file_out, [])
                BEDTOOLS_MERGE(SORT_FOOTPRINTS.out.sorted)
                ch_footprints = BEDTOOLS_MERGE.out.bed
            } else {
                ch_footprints = FILTER_CLOSEST.out.output
            }
        }

   

        if (params.blacklist) {
            SUBTRACT_BLACKLIST(ch_footprints.map{ meta, bed_file -> [meta, bed_file, params.blacklist]})
            ch_blacklisted = SUBTRACT_BLACKLIST.out.bed
        } else {
            ch_blacklisted = ch_footprints
        }

        STARE (
            ch_footprints,
            params.pwm,
            params.gtf,
            params.fasta,
            params.stare_window_size,
            params.blacklist,
            params.stare_decay
        )

        ch_affinities = STARE.out.affinities

        REMOVE_VERSION_ANNOTATIONS(ch_affinities)

        if (params.merge_peaks)
        {
            ch_mean_affinities = REMOVE_VERSION_ANNOTATIONS.out
        } else {
            AFFINITY_MEAN(
                REMOVE_VERSION_ANNOTATIONS.out
                    .map{ meta, affinity_file -> 
                            [
                                [
                                    id: meta.state + "_" + meta.antibody,
                                    state: meta.state,
                                    antibody: meta.antibody
                                ], 
                                affinity_file
                            ]
                        }
                    .groupTuple(),
                "mean"
            )
            ch_mean_affinities = AFFINITY_MEAN.out
        }

        ch_affinity_antibody = ch_mean_affinities
            .map{ meta, affinity_file -> [meta, meta.antibody, affinity_file]}

        ch_affinity_pairings = ch_affinity_antibody
            .combine(
                ch_affinity_antibody,
                by: 1
            )
            .filter{
                antibody, meta1, affinites1, meta2, affinites2 -> 
                    meta1.state < meta2.state
            }
            .map{
                antibody, meta1, affinites1, meta2, affinites2 -> 
                    [[  id: meta1.state + ":" + meta2.state + "_" + antibody, 
                        state1: meta1.state, 
                        state2: meta2.state, 
                        antibody: antibody
                    ],
                    [affinites1, affinites2]]
            }

        AFFINITY_RATIO(ch_affinity_pairings, "ratio")
        AFFINITY_SUM(ch_affinity_pairings, "sum")

    emit:
        affinity_ratio = AFFINITY_RATIO.out
        affinity_sum = AFFINITY_SUM.out
        peaks = ch_blacklisted
}