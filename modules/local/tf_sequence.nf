process TF_SEQUENCE {
    conda "python logomaker matplotlib"
    container "tfprio-python"

    input:
        path(tfs)

    output:
        path("*.svg")
    
    script:
        """
        jaspar.py -t $tfs
        """
}