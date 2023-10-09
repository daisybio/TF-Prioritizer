include { INTEGRATE_DATA } from             '../../modules/local/dynamite/integrate'
include { PREPARE_FOR_CLASSIFICATION } from '../../modules/local/dynamite/prepare_for_classification'
include { DYNAMITE as RUN_DYNAMITE } from                   '../../modules/local/dynamite/dynamite'
include { GAWK as FILTER } from '../../modules/nf-core/gawk/main'

workflow DYNAMITE {
    take:
        ch_data
        ofolds
        ifolds
        alpha
        randomize

    main:
        INTEGRATE_DATA (ch_data)
        PREPARE_FOR_CLASSIFICATION (INTEGRATE_DATA.out)

        RUN_DYNAMITE (
            PREPARE_FOR_CLASSIFICATION.out,
            ofolds, ifolds, alpha, randomize
        )

        
        FILTER (
            RUN_DYNAMITE.out,
            []
        )

    emit:
        FILTER.out.output
}