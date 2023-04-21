process MANN_WHITNEY_U {
    conda "pandas matplotlib numpy scipy"
    container "tfprio-python"

    input:
        tuple val(group1), val(group2), val(hm), path(tftg_scores)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}.tsv")
    script:
        """
        mann_whitney_u.py --input $tftg_scores --output ${group1}:${group2}_${hm}.tsv
        """
}