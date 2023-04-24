process TF_SEQUENCE {
    conda "python logomaker matplotlib"
    container "tfprio-python"

    input:
        path(tfs)

    output:
        path("*_jaspar")
    
    script:
        """
        jaspar.py -t $tfs
        """
}