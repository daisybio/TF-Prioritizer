include { COLLECT_TF_DATA } from '../../modules/local/collect'
include { TOP_TARGET_GENES } from '../../modules/local/distributionAnalysis/top_target_genes'
include { BIOPHYSICAL_MODELS } from '../../modules/local/biophysical_models'
include { HEATMAPS } from '../../modules/local/distributionAnalysis/heatmaps'
include { TF_SEQUENCE } from '../../modules/local/tf_sequence'
include { COLLECT_HEATMAPS } from '../../modules/local/collect'
include { REPORT as CREATE_REPORT } from '../../modules/local/report'
include { COLLECT } from '../../modules/local/collect'

workflow PREPARE_TFGROUPS {
    take:
        ch_tfgroups // tf-group
        top_target_genes // value
        ch_affinity_sums
        ch_counts_per_sample
        ch_rnaseq_samplesheet
        pwm_file
        ch_ensg_map

    main:
        TOP_TARGET_GENES (
            ch_tfgroups,
            top_target_genes,
            ch_affinity_sums,
            ch_counts_per_sample
        ) // group1, group2, hm, topTargetGenes

        ch_heatmaps = COLLECT_HEATMAPS (
            HEATMAPS (
                TOP_TARGET_GENES.out,
                ch_counts_per_sample,
                ch_rnaseq_samplesheet,
                ch_ensg_map
            )   .transpose() // group1, group2, hm, heatmap
                .map { it[3] }
                .collect()
        )   .flatten()
            .map { [it.name.replaceAll(/_heatmaps$/, ''), it] }

        BIOPHYSICAL_MODELS (
            ch_tfgroups,
            pwm_file
        )

        ch_bio_logos = BIOPHYSICAL_MODELS.out.plots.flatten().map { [it.name.replaceAll(/.png$/, ''), it] }
        ch_bio_models = BIOPHYSICAL_MODELS.out.pwms.flatten().map { [it.name.replaceAll(/.pwm$/, ''), it] }

        ch_jaspar = TF_SEQUENCE (ch_tfgroups)
            .flatten()
            .map { [it.name.replaceAll(/_jaspar$/, ''), it] }

        ch_biophysical = ch_bio_logos
            .combine(ch_bio_models, by: 0) // tf-group, logo, model

        ch_logos = ch_biophysical
            .combine(ch_jaspar, by: 0) // tf-group, biophysical_logo, model, jaspar_logos

        ch_plots = ch_logos
            .combine(ch_heatmaps, by: 0) // tf-group, biophysical_logo, model, jaspar_logos, heatmaps

        COLLECT_TF_DATA(ch_plots)

    emit:
        COLLECT_TF_DATA.out
}

workflow REPORT {
    take:
        ch_expression
        ch_ranks
        ch_tfgroups
        ch_groups_tf_map
        ch_tf_ensg_map

    main:
        COLLECT (
            ch_expression,
            ch_ranks,
            ch_tfgroups,
            ch_groups_tf_map,
            ch_tf_ensg_map
        )

        CREATE_REPORT ( COLLECT.out )
}