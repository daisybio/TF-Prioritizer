include { INTEGRATE_DATA } from '../../modules/local/dynamite'
include { PREPARE_FOR_CLASSIFICATION } from '../../modules/local/dynamite'
include { DYNAMITE as RUN_DYNAMITE } from '../../modules/local/dynamite'
include { DYNAMITE_FILTER } from '../../modules/local/dynamite'

workflow DYNAMITE {
    take:
        ch_affinity_ratios // group1, group2, hm, affinityRatios
        ch_diff_expression // group1, group2, deseq2-file
        ofolds // value
        ifolds // value
        alpha // value
        randomize // value
        min_regression // value

    main:
        ch_integration = ch_affinity_ratios
            .combine(ch_diff_expression, by: [0, 1]) // group1, group2, hm, affinityRatios, diffExpression

        INTEGRATE_DATA (ch_integration)
        PREPARE_FOR_CLASSIFICATION (INTEGRATE_DATA.out)

        RUN_DYNAMITE (
            PREPARE_FOR_CLASSIFICATION.out,
            ofolds, ifolds, alpha, randomize
        )

        ch_dynamite = DYNAMITE_FILTER (
            RUN_DYNAMITE.out,
            min_regression
        )

    emit:
        ch_dynamite
}