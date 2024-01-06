process ROSE_OUTPUT_TO_BED {
    tag "${meta.state}"

    label "process_single"
    container "quay.io/biocontainers/pandas:1.4.3"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${gff.baseName}.bed")
    
    script:
    """
    rose_output_to_bed.py --gff ${gff} --bed ${gff.baseName}.bed
    """

    stub:
    """
    touch ${gff.baseName}.bed
    """

}