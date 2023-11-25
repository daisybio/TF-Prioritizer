process COPY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
        tuple val(meta), path(original), val(copy_name)
    
    output:
        tuple val(meta), path(original), path(copy_name)
    
    script:
    if (! (original.name == copy_name))
        """
        ln -s ${original} ${copy_name}
        """
    else
        """
        echo "No need to copy ${original} to ${copy_name}"
        """
}