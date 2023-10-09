process INTEGRATE_DATA {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(affinity_ratios), path(diff_expression)

    output:
        tuple val(meta), path("${meta.id}_integrated.tsv")

    script:
        """
        python3 $projectDir/lib/tepic/DYNAMITE/Scripts/integrateData.py \\
            $affinity_ratios \\
            $diff_expression \\
            "${meta.id}_integrated.tsv" \\
            --expressionC 2
        """
}