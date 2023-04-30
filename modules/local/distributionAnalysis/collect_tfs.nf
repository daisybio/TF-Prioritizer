process COLLECT_TFS {
    container 'tfprio-python'

    input:
        path(rankings)

    output:
        path('tf_groups.txt'), emit: groups
        path('tf_groups.json'), emit: groups_tfs_map
        path('tfs.txt'), emit: tfs

    script:
        // Concatenate all files while removing headers
        """
        awk '\$1!="sum" {print \$1}' *_ranked.tsv | sort -u > tf_groups.txt
        tf_group_splitting.py -i tf_groups.txt -o tfs.txt -j tf_groups.json
        """
}