include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"

workflow EHMM {
    take:
        background_bed
        enhancers_bed
        promoters_bed
        gtf

        genomic_region_size
        train_split
        random_seed
        top_quantile
        n_samples

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
}