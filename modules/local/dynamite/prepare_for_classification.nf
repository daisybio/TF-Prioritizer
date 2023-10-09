process PREPARE_FOR_CLASSIFICATION {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-optparse:1.6.0' :
        'biocontainers/r-optparse:1.6.0' }"

    input:
        tuple val(meta), path(integrated_data)
    
    output:
        tuple val(meta), path("${meta.id}_classification.tsv")
    
    script:
        """
        Rscript $projectDir/lib/tepic/DYNAMITE/Scripts/prepareForClassification.R \\
            $integrated_data \\
            "${meta.id}_classification.tsv"
        """
}