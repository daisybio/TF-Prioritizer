process CREATE_ANNDATA {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(counts), path(design)

    output:
        tuple val(meta), path("*.h5ad")

    script:
    """
        #!/usr/bin/env python3

        import anndata as ad
        import pandas as pd

        counts = pd.read_csv("$counts", index_col=0, sep="\\t")
        design = pd.read_csv("$design", index_col=0)

        # Create anndata object
        adata = ad.AnnData(X=counts.T.values, obs=design)
        adata.var.index = counts.index

        # Add TPMs as a layer
        adata.layers["TPM"] = counts.T.div(counts.T.sum(axis=1), axis=0) * 1e6
        adata.layers["raw_counts"] = counts.T.values

        # Save anndata object
        adata.write("counts.h5ad")
    """
}