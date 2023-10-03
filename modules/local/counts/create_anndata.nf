process CREATE_ANNDATA {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(counts), path(design)
        val(index_col)

    output:
        tuple val(meta), path("*.h5ad"), emit: anndata
        tuple val(meta), path("genes.txt"), emit: genes

    script:
    """
        create_anndata.py --counts ${counts} --metadata ${design} --output ${meta.id}.h5ad --index_col ${index_col}
    """
}