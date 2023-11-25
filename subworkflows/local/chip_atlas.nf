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
        antigens

    main:
        DOWNLOAD()

        // Filter by genome, tissues and keep only entry with lowest threshold
        FILTER(
            DOWNLOAD.out.map{it -> [[id: "chip_atlas"], it]},
            params.genome,
            params.chip_atlas_tissues
        )

        ch_entries = FILTER.out.map{ it[1] }.splitCsv(header: true, sep: "\t")
        
        ch_ehmm_prep = ch_entries
            .map{
                [antigen_key: (it['antigen_class'] + "_" + it['antigen']).strip("_")] + it
            }
            .filter{ antigens.contains(it['antigen_key'])}
            .map{ [[id: it['antigen_key'].replace("Histone_", "")], file(it['file_url'])] }

        CLEAN_ALL(ch_ehmm_prep, [])
        REMOVE_CHR(CLEAN_ALL.out.output, [])

        BEDTOOLS_BEDTOBAM(REMOVE_CHR.out.output, fai)
        SAMTOOLS_SORT(BEDTOOLS_BEDTOBAM.out.bam)
        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

        ch_ehmm = REMOVE_CHR.out.output
                    .join(SAMTOOLS_SORT.out.bam)
                    .join(SAMTOOLS_INDEX.out.bai)

    emit:
        ehmm = ch_ehmm
}
