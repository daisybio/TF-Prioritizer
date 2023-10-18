process DOWNLOAD {
    tag "chip_atlas"
    label 'process_single'

    conda "conda-forge::curl=8.4.0 conda-forge::unzip=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/curl:7.80.0' :
        'biocontainers/curl:7.80.0' }"


    output:
        path("chip_atlas_file_list.csv")
    
    script:
    """
        curl -o overview.zip http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv
        unzip overview.zip
    """
}