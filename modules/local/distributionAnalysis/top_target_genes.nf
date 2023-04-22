process TOP_TARGET_GENES {
    conda "pandas scipy"
    container "tfprio-python"

    input:
        path(tfs)
        tuple val(group1), val(group2), val(hm), path(affinities)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_ttg.tsv")

    script:
        """
        top_target_genes.py -t $tfs -a $affinities -o ${group1}:${group2}_${hm}_ttg.tsv
        """
}