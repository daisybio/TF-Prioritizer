process BED_TO_GFF {
    tag "${meta.state}"

    label "process_single"
    container "registry.hub.docker.com/library/python:3.11"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${bed.baseName}.gff")

    script:
    """
    bed_to_gff.py --bed "${bed}" --gff "${bed.baseName}.gff"
    """

    stub:
    """
    touch "${bed.baseName}.gff"
    """

}