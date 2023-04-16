process SEQUENCE_TO_BED {
    conda 'bioconda::bedtools'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(hm), val(group), path(binding_sequences)
    val affinity_cutoff

    output:
    tuple val(hm), val(group), path("*.bed")

    script:
    """
    binding_sequences.sh $binding_sequences $affinity_cutoff
    """
}