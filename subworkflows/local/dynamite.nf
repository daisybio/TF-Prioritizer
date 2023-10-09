include { INTEGRATE_DATA } from             '../../modules/local/dynamite/integrate'
include { PREPARE_FOR_CLASSIFICATION } from '../../modules/local/dynamite/prepare_for_classification'
include { DYNAMITE as RUN_DYNAMITE } from                   '../../modules/local/dynamite/dynamite'
//include { DYNAMITE_FILTER } from '../../modules/local/dynamite'

workflow DYNAMITE {
    take:
        ch_data
        ofolds
        ifolds
        alpha
        randomize
        min_regression

    main:
        INTEGRATE_DATA (ch_data)
        PREPARE_FOR_CLASSIFICATION (INTEGRATE_DATA.out)

        RUN_DYNAMITE (
            PREPARE_FOR_CLASSIFICATION.out,
            ofolds, ifolds, alpha, randomize
        )

        /*
        ch_dynamite = DYNAMITE_FILTER (
            RUN_DYNAMITE.out,
            min_regression
        )

    emit:
        ch_dynamite*/
}