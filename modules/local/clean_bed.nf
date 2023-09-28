process CLEAN_BED {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
        tuple val(meta), path(file)

    output:
        tuple val(meta), path("${meta.id}_cleaned.bed")

    script:
        """
        cut -f1-3 ${file} | awk -v sample=${meta.id} '{print \$0"\t"sample}' > ${meta.id}_cleaned.bed
        """
}