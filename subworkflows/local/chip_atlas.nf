include { DOWNLOAD } from "../../modules/local/chip_atlas/download"
include { FILTER } from "../../modules/local/chip_atlas/filter"
include { GAWK as CLEAN_ALL } from "../../modules/nf-core/gawk"
include { BEDTOOLS_BEDTOBAM } from "../../modules/local/bedtools/bedtobam/main"
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { GAWK as REMOVE_CHR } from "../../modules/nf-core/gawk"

workflow CHIP_ATLAS {
    take:
        fai

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
        REMOVE_CHR(CLEAN_ALL.out.output, [])

        BEDTOOLS_BEDTOBAM(REMOVE_CHR.out.output, fai)
        SAMTOOLS_SORT(BEDTOOLS_BEDTOBAM.out.bam)
        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    emit:
        all = REMOVE_CHR.out.output
        all_bam = SAMTOOLS_SORT.out.bam
        all_bai = SAMTOOLS_INDEX.out.bai
}
