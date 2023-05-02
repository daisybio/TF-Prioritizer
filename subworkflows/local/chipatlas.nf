include { FETCH_LINKS } from '../../modules/local/chip_atlas'
include { FETCH_BED } from '../../modules/local/chip_atlas'
include { BED_TO_BAM } from '../../modules/local/bed_to_bam'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index'
include { REMOVE_CHR_PREFIX } from '../../modules/local/remove_chr_prefix'

workflow CHIP_ATLAS {
    take:
        tissue_types
        tfs
        genome
        threshold
        ch_chromosome_lengths
    
    main:
        ch_links = FETCH_LINKS(tissue_types, tfs, genome, threshold)
            .splitCsv(header: true, sep: '\t')
            .map { [it + ["id": it.tf + '_' + it.tissue], it.url ]}

        FETCH_BED(ch_links)

        REMOVE_CHR_PREFIX(FETCH_BED.out)

        BED_TO_BAM(REMOVE_CHR_PREFIX.out, ch_chromosome_lengths)

        SAMTOOLS_SORT(BED_TO_BAM.out)

        SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

        ch_collected = REMOVE_CHR_PREFIX.out
            .combine(SAMTOOLS_SORT.out.bam, by: 0)
            .combine(SAMTOOLS_INDEX.out.bai, by: 0)

    emit:
        ch_collected
}