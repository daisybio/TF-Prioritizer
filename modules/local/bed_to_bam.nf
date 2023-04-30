process BED_TO_BAM {
    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
       tuple val(meta), path(bed)
       path(chrom_sizes)

    output:
        tuple val(meta), path("${meta.id}.bam")

    script:
    """
    bedtools bedtobam -i ${bed} -g ${chrom_sizes} > ${meta.id}.bam
    """
}