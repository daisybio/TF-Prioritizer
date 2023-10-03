process FILTER_ANNDATA {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(anndata)
        val(min_count)
        val(min_tpm)

    output:
        tuple val(meta), path("*.filtered.h5ad")

    script:
    """
        #!/usr/bin/env python3
        import anndata as ad

        adata = ad.read_h5ad("${anndata}")

        # Filter genes
        adata = adata[:, adata.X.sum(axis=0) >= ${min_count}]
        adata = adata[:, adata.layers["TPM"].sum(axis=0) >= ${min_tpm}]

        adata.write_h5ad("${meta.id}.filtered.h5ad")
    """
}