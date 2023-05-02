process CONCATENATE {
    container 'ubuntu'

    input:
        tuple val(meta), path(files)

    output:
        tuple val(meta), path("${meta.id}_concatenated.${meta.ext}")

    script:
        """
        cat ${files} > ${meta.id}_concatenated.${meta.ext}
        """
}