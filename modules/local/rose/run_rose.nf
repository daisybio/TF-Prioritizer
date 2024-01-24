process RUN_ROSE {
    tag "${meta.state}"

    label "process_single"
    container "registry.hub.docker.com/library/python:3.11"

    input:
    tuple val(meta), path(gff)
    path ucsc_file

    output:
    tuple val(meta), path("${gff.baseName}_STITCHED.gff")

    script:
    """
    rose.py \
    -g ${ucsc_file} \
    -i ${gff} \
    -o ${gff.baseName}_STITCHED.gff \
    -s 12500 \
    -t 2500
    """

    stub:
    """
    touch "${gff.baseName}_STITCHED.gff"
    """
}