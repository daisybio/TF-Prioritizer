include { TFTG_SCORE } from '../../modules/local/distributionAnalysis/tftg_scores'
include { RANKING as CREATE_RANKING } from '../../modules/local/distributionAnalysis/ranking'
include { COLLECT_TFS } from '../../modules/local/distributionAnalysis/collect_tfs'
include { SYMBOL_ENSG } from '../../modules/local/symbol_ensg'
include { COLLECT_RANKS } from '../../modules/local/collect'


workflow RANKING {
    take:
        ch_dynamite
        ch_diff_expression // group1, group2, deseq2-file
        ch_affinity_sums // group1, group2, hm, affinitySum-file
        ch_map

    main:
        ch_tftg = ch_dynamite
            .combine(ch_diff_expression, by: [0, 1]) // group1, group2, hm, coefficients, diffExpression
            .combine(ch_affinity_sums, by: [0, 1, 2]) // group1, group2, hm, coefficients, diffExpression, affinitySums

        TFTG_SCORE (ch_tftg) // group1, group2, hm, tftgScore
        CREATE_RANKING (TFTG_SCORE.out) // group1, group2, hm, ranking

        COLLECT_TFS (
            CREATE_RANKING.out.map { it[3] } .collect()
        ) // tfs

        ch_tf_ensg = COLLECT_TFS.out.tfs
            .splitCsv(header: false)
            .map { it[0] }
            .combine(ch_map, by: 0) // tf, ensg

        tf_ensg_map = ch_tf_ensg
            .collectFile(name: "tf_map.tsv", storeDir: "${params.outdir}/results", newLine: true) { it.join('\t')}

        SYMBOL_ENSG(tf_ensg_map)

        ensgs = ch_tf_ensg
            .map { it[1] }
            .collect()

        ch_ranks = COLLECT_RANKS (
            CREATE_RANKING.out.map { it[3] } .collect()
        )

    emit:
        tfs = COLLECT_TFS.out.tfs
        groups = COLLECT_TFS.out.groups
        group_tfs_map = COLLECT_TFS.out.groups_tfs_map
        tf_ensg_map = SYMBOL_ENSG.out
        ensgs = ensgs
        ranks = ch_ranks

}