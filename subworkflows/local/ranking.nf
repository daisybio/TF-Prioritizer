include { TFTG_SCORE } from '../../modules/local/ranking/TFTG'
include { RANKING as CREATE_RANKING } from '../../modules/local/ranking/ranking'
// include { COLLECT_TFS } from '../../modules/local/distributionAnalysis/collect_tfs'
// include { SYMBOL_ENSG } from '../../modules/local/symbol_ensg'
// include { COLLECT_RANKS } from '../../modules/local/collect'


workflow RANKING {
    take:
        ch_affinitySum_deseq_dynamite

    main:
        TFTG_SCORE (ch_affinitySum_deseq_dynamite)
        CREATE_RANKING (TFTG_SCORE.out)
        /*
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
    */
}