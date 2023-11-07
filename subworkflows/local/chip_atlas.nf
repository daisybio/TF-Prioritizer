include { DOWNLOAD } from "../../modules/local/chip_atlas/download"
include { FILTER } from "../../modules/local/chip_atlas/filter"

include { GAWK as CLEAN_ALL } from "../../modules/nf-core/gawk"

workflow CHIP_ATLAS {
    main:
        DOWNLOAD()
        FILTER(
            DOWNLOAD.out.map{it -> [[id: "chip_atlas"], it]},
            params.genome,
            params.chip_atlas_tissues
        )

        ch_entries = FILTER.out.map{ it[1] }.splitCsv(header: true, sep: "\t")

        ch_all = ch_entries
            .filter{ it['antigen'] == "" && it['cell_type'] == "" }
            .map{ file(it['file_url']) }
            .first()
            .map(it -> [[id: "chip_atlas_all"], it])

        CLEAN_ALL(ch_all, [])

    emit:
        all = CLEAN_ALL.out.output
}