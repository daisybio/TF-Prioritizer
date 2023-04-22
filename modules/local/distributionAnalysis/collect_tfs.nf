process COLLECT_TFS {
    container 'ubuntu:22.04'

    input:
        path(rankings)

    output:
        path('tfs.txt')

    script:
        // Concatenate all files while removing headers
        """
        awk '\$1!="sum" {print \$1}' *_ranked.tsv | sort -u > tfs.txt
        """
}