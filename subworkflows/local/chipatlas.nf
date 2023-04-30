include { FETCH_LINKS } from '../../modules/local/chip_atlas'
include { FETCH_BED } from '../../modules/local/chip_atlas'
include { BED_TO_BAM } from '../../modules/local/bed_to_bam'

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
            .map { [it, it.url ]}

        FETCH_BED(ch_links)

        BED_TO_BAM(FETCH_BED.out, ch_chromosome_lengths)
}