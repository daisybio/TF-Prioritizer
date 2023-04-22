process RANKING {
    conda "pandas scipy"
    container "tfprio-python"

    input:
        tuple val(group1), val(group2), val(hm), path(tftg_scores)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}.tsv")
    script:
        """
        ranking.py --input $tftg_scores --output ${group1}:${group2}_${hm}.tsv
        """
}