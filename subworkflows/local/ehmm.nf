include { SCRAPE_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { DOWNLOAD_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { UPLIFT } from "../../modules/local/uplift"
include { EPD_NEW } from "../../modules/local/ehmm/epd_new"
include { CONCATENATE as CONCATENATE_EH_ATLAS } from "../../modules/local/concatenate"
include { CONCATENATE as CONCATENATE_CHIP_ATLAS } from "../../modules/local/concatenate"
include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"
include { LEARN_MODEL as LEARN_BACKGROUND_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_ENHANCER_MODEL } from "../../modules/local/ehmm/learn_model"
include { LEARN_MODEL as LEARN_PROMOTER_MODEL } from "../../modules/local/ehmm/learn_model"
include { REMOVE_CHR_PREFIX as REMOVE_EH_CHR_PREFIX } from "../../modules/local/remove_chr_prefix"
include { REMOVE_CHR_PREFIX as REMOVE_EPDNEW_CHR_PREFIX } from "../../modules/local/remove_chr_prefix"
include { CONSTRUCT_MODEL } from "../../modules/local/ehmm/construct_model"

workflow EHMM {
    take:
        ch_chipatlas // meta, bed, bam, bai
        ch_chromosome_lengths // chromosome-lengths-file
        genome
        tissues

    main:
        SCRAPE_EH_ATLAS (tissues, genome)

        ch_beds = DOWNLOAD_EH_ATLAS (
            SCRAPE_EH_ATLAS.out.splitCsv(header: true, sep: "\t").map { [it + ['id': it.tf], it.link] },
        )

        ch_beds = REMOVE_EH_CHR_PREFIX (
            ch_beds
        )

        ch_beds = UPLIFT (
            ch_beds,
            genome
        )

        CONCATENATE_EH_ATLAS (
            ch_beds.map { it[1] }.collect()
        )

        CONCATENATE_CHIP_ATLAS (
            ch_chipatlas.map { it[1] }.collect()
        )

        epd_new = EPD_NEW(genome)

        epd_new = REMOVE_EPDNEW_CHR_PREFIX (
            epd_new
        )

        SPLIT_DATASETS (
            CONCATENATE_CHIP_ATLAS.out,
            CONCATENATE_EH_ATLAS.out,
            epd_new.map { it[1] },
            Channel.value(params.tepic_gtf),
            Channel.value(params.ehmm_genomic_region_size),
            Channel.value(params.ehmm_train_split),
            Channel.value(params.ehmm_random_seed),
            Channel.value(params.ehmm_top_quantile),
            Channel.value(params.ehmm_n_samples)
        )

        bams = ch_chipatlas.map { it[2] }.collect()
        bais = ch_chipatlas.map { it[3] }.collect()

        LEARN_BACKGROUND_MODEL (
            SPLIT_DATASETS.out.trainBackground,
            bams,
            bais,
            Channel.value(params.ehmm_n_states),
            Channel.value(params.ehmm_n_bins),
            Channel.value(params.ehmm_pseudocount),
        )

        LEARN_ENHANCER_MODEL (
            SPLIT_DATASETS.out.trainEnhancers,
            bams,
            bais,
            Channel.value(params.ehmm_n_states),
            Channel.value(params.ehmm_n_bins),
            Channel.value(params.ehmm_pseudocount),
        )

        LEARN_PROMOTER_MODEL (
            SPLIT_DATASETS.out.trainPromoters,
            bams,
            bais,
            Channel.value(params.ehmm_n_states),
            Channel.value(params.ehmm_n_bins),
            Channel.value(params.ehmm_pseudocount),
        )

        CONSTRUCT_MODEL (
            LEARN_BACKGROUND_MODEL.out.rdata,
            LEARN_ENHANCER_MODEL.out.rdata,
            LEARN_PROMOTER_MODEL.out.rdata
        )
}