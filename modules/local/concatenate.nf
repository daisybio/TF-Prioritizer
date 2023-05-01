process CONCATENATE {
    container 'ubuntu'

    input:
        path(files)

    output:
        path("*concatenated.txt")

    script:
    def prefix = task.ext.prefix + "_" ?: ""
        """
        cat ${files} > ${prefix}concatenated.txt
        """
}