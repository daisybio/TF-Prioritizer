include { TFTG_SCORE } from '../../modules/local/ranking/TFTG'
include { RANKING as CREATE_RANKING } from '../../modules/local/ranking/ranking'
include { COLLECT_TFS } from '../../modules/local/ranking/collect_tfs'


workflow RANKING {
    take:
        ch_affinitySum_deseq_dynamite

    main:
        TFTG_SCORE (ch_affinitySum_deseq_dynamite)
        CREATE_RANKING (TFTG_SCORE.out)

        ch_rankings = CREATE_RANKING.out.map { meta, ranking -> ranking }.collect()

        COLLECT_TFS (ch_rankings)

    emit:
        tfs = COLLECT_TFS.out.tfs
        groups = COLLECT_TFS.out.tfgroup_tfs
        ranks = CREATE_RANKING.out
}