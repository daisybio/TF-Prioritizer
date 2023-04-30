include { SCRAPE_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { DOWNLOAD_EH_ATLAS } from "../../modules/local/ehmm/eh_atlas"
include { UPLIFT } from "../../modules/local/uplift"
include { EPD_NEW } from "../../modules/local/ehmm/epd_new"

workflow EHMM {
    take:
        ch_chipatlas // tf, tissue, chipatlas-file
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

        EPD_NEW(genome)
}