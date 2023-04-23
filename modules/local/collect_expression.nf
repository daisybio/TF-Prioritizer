process COLLECT_EXPRESSION {
    conda "pandas"
    container "tfprio-python"

    input:
        path(tpm)
        path(deseq)
        path(counts)
        val(ensgs)

    output:
        path("*.json")

    script:
        """
        collect_expression.py -t $tpm -c $counts -d $deseq -e ${ensgs.join(" ")}
        """
}