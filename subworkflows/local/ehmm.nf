include { SCRAPE_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { DOWNLOAD_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { UPLIFT } from "../../modules/local/uplift"
include { EPD_NEW } from "../../modules/local/ehmm/epd_new"
include { CONCATENATE as CONCATENATE_EH_ATLAS } from "../../modules/local/concatenate"
include { CONCATENATE as CONCATENATE_CHIP_ATLAS } from "../../modules/local/concatenate"
include { SPLIT_DATASETS } from "../../modules/local/ehmm/split_datasets"

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

        EPD_NEW(genome)

        SPLIT_DATASETS (
            CONCATENATE_CHIP_ATLAS.out,
            CONCATENATE_EH_ATLAS.out,
            EPD_NEW.out,
            Channel.value(params.tepic_gtf),
            Channel.value(params.ehmm_genomic_region_size),
            Channel.value(params.ehmm_train_split),
            Channel.value(params.ehmm_random_seed),
            Channel.value(params.ehmm_top_quantile),
            Channel.value(params.ehmm_sample_number)
        )
}