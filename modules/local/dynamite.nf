process INTEGRATE_DATA {
    container 'python'

    input:
        tuple val(group1), val(group2), val(hm), path(affinity_ratios), path(diff_expression)

    output:
        path("${group1}:${group2}_${hm}_integrated.tsv")

    script:
        """
        python3 $projectDir/assets/tepic/DYNAMITE/Scripts/integrateData.py \\
            $affinity_ratios \\
            $diff_expression \\
            "${group1}:${group2}_${hm}_integrated.tsv" \\
            --geneIDs 0 --expressionC 1
        """
}