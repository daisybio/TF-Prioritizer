include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"
include { CAT_CAT as CAT_BACKGROUND } from "../../modules/nf-core/cat/cat/main"

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

        LEARN_BACKGROUND_MODEL (
            SPLIT_DATASETS.out.trainBackground.map{[[id: "train_background", model: "BackgroundModel"], it]},
            background.map{ it[2] }.collect().map{ [[id: "background_bam"], it]},
            background.map{ it[3] }.collect().map{ [[id: "background_bai"], it]},
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers.map{[[id: "train_enhancers", model: "EnhancerModel"], it]},
            background.map{ it[2] }.collect().map{ [[id: "background_bam"], it]},
            background.map{ it[3] }.collect().map{ [[id: "background_bai"], it]},
            n_states,
            n_bins,
            pseudocount,
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters.map{[[id: "train_promoters", model: "PromoterModel"], it]},
            background.map{ it[2] }.collect().map{ [[id: "background_bam"], it]},
            background.map{ it[3] }.collect().map{ [[id: "background_bai"], it]},
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
