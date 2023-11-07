include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"

workflow EHMM {
    take:
        background_bed
        background_bam
        background_bai

        enhancers_bed
        promoters_bed
        gtf

        genomic_region_size
        train_split
        random_seed
        top_quantile
        n_samples

        n_states
        n_bins
        pseudocount

    main:
        SPLIT_DATASETS(
            background_bed,
            enhancers_bed,
            promoters_bed,
            gtf,
            genomic_region_size,
            train_split,
            random_seed,
            top_quantile,
            n_samples
        )

        LEARN_BACKGROUND_MODEL (
            SPLIT_DATASETS.out.trainBackground.map{[[id: "train_background", model: "background"], it]},
            background_bam,
            background_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers.map{[[id: "train_enhancers", model: "enhancer"], it]},
            background_bam,
            background_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters.map{[[id: "train_promoters", model: "promoter"], it]},
            background_bam,
            background_bai,
            n_states,
            n_bins,
            pseudocount,
        )

        CONSTRUCT_MODEL (
            LEARN_BACKGROUND_MODEL.out.rdata,
            LEARN_ENHANCER_MODEL.out.rdata,
            LEARN_PROMOTER_MODEL.out.rdata
        )
}
