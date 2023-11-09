include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"

workflow EHMM {
    take:
        background

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
            background.map{ it[1] }.collect(),
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
            background.map{ it[2] }.collect(),
            background.map{ it[3] }.collect(),
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers.map{[[id: "train_enhancers", model: "enhancer"], it]},
            background.map{ it[2] }.collect(),
            background.map{ it[3] }.collect(),
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters.map{[[id: "train_promoters", model: "promoter"], it]},
            background.map{ it[2] }.collect(),
            background.map{ it[3] }.collect(),
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
