process USE_SYMBOL_MAP {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(anndata)
        tuple val(meta2), path(symbol_map)

    output:
        tuple val(meta), path("*.symbols.h5ad")

    script:
    """
        #!/usr/bin/env python3
        import anndata as ad
        import json

        adata = ad.read_h5ad("${anndata}")

        with open("${symbol_map}") as f:
            symbol_map = json.load(f)
        
        adata.var.index = [symbol_map.get(gene, gene) for gene in adata.var.index]

        adata.write_h5ad("${meta.id}.symbols.h5ad")
    """
}