include { BED_TO_GFF } from "../../modules/local/rose/bed_to_gff"
include { REFORMAT_GFF } from "../../modules/local/rose/reformat_gff"
include { RUN_ROSE } from "../../modules/local/rose/run_rose"
include { GAWK as ROSE_OUTPUT_TO_BED } from "../../modules/nf-core/gawk/main"

workflow ROSE {
    take:
        input_bed
        ucsc_file
    
    main:
        BED_TO_GFF(input_bed)

        REFORMAT_GFF(BED_TO_GFF.out)

        RUN_ROSE(REFORMAT_GFF.out, ucsc_file)

        ROSE_OUTPUT_TO_BED(RUN_ROSE.out, [])

    emit:
        stitched_enhancers = ROSE_OUTPUT_TO_BED.out.output
}