process REFORMAT_GFF {
    tag "${meta.state}"
    
    label "process_single"
    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    tuple val(meta), path(input_gff)

    output:
    tuple val(meta), path("${input_gff.baseName}_rose.gff")

    script:
    """
    reformat_gff.py \
        --input ${input_gff} \
        --output ${input_gff.baseName}_rose.gff
    """

    stub:
    """
    touch "${input_gff.baseName}_rose.gff"
    """
}