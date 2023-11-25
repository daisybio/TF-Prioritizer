include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"
include { APPLY_MODEL } from "../../modules/local/ehmm/apply_model"
include { CAT_CAT as CAT_BACKGROUND } from "../../modules/nf-core/cat/cat/main"
include { CAT_CAT as MERGE_ASSAYS } from "../../modules/nf-core/cat/cat/main"
include { COPY as RENAME_BAMS } from "../../modules/local/copy"

workflow EHMM {
    take:
        background
        enhancers_bed
        promoters_bed
        peaks
        bams

        gtf
        index

        genomic_region_size
        train_split
        random_seed
        top_quantile
        n_samples

        n_states
        n_bins
        pseudocount

    main:
        CAT_BACKGROUND(
            background.map{ it[1] }.collect().map { [[id: 'background'], it] }
        )

        SPLIT_DATASETS(
            CAT_BACKGROUND.out.file_out,
            enhancers_bed,
            promoters_bed,
            gtf,
            genomic_region_size,
            train_split,
            random_seed,
            top_quantile,
            n_samples
        )

        ch_bam_bai = background.map{ [it[2], it[3]] }.flatten().collect()

        LEARN_BACKGROUND_MODEL (
            SPLIT_DATASETS.out.trainBackground.map{[[id: "train_background", model: "BackgroundModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers.map{[[id: "train_enhancers", model: "EnhancerModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters.map{[[id: "train_promoters", model: "PromoterModel"], it]},
            ch_bam_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        CONSTRUCT_MODEL (
            LEARN_BACKGROUND_MODEL.out.rdata,
            LEARN_ENHANCER_MODEL.out.rdata,
            LEARN_PROMOTER_MODEL.out.rdata
        )

        MERGE_ASSAYS (
            peaks.map{meta, file_ -> [[id: meta["state"], state: meta["state"]], file_]}.groupTuple()
        )

        RENAME_BAMS (
            bams.map{meta, original -> [meta, original, meta["assay"]]}
        )

        ch_bams = RENAME_BAMS.out
                    .map{ meta, original, directory_ -> [[id: meta["state"], state: meta["state"]], directory_]}
                    .groupTuple()

        ch_combined = MERGE_ASSAYS.out.file_out.join(ch_bams)

        APPLY_MODEL (
            ch_combined,
            CONSTRUCT_MODEL.out.model,
            index,
            n_bins,
            pseudocount
        )

        ch_enhancers = APPLY_MODEL.out.enhancers
        ch_promoters = APPLY_MODEL.out.promoters
}
