process TFTG_SCORE {
    conda "pandas"
    container "tfprio-python"

    input:
        tuple val(group1), val(group2), val(hm), path(regression_coefficients), path(diff_expression), path(affinity_sums)

    output:
        tuple val(group1), val(group2), val(hm), path("${group1}:${group2}_${hm}_tftg.tsv")

    script:
        """
        tftg_scores.py -l $diff_expression -c $regression_coefficients -a $affinity_sums -o ${group1}:${group2}_${hm}_tftg.tsv
        """
}