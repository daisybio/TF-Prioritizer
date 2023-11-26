process USE_SYMBOL_MAP {
    tag "$meta.id"
    label "process_low"

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(dataframe)
        tuple val(meta2), path(symbol_map)

    output:
        tuple val(meta), path("*.symbols.tsv")

    script:
    """
        #!/usr/bin/env python3
        import pandas as pd
        import json

        df = pd.read_csv("${dataframe}", sep="\\t", index_col=0)
        df_mapping = pd.read_csv("${symbol_map}", sep="\\t", index_col=0)
        
        df = df.rename(index=df_mapping["gene_name"])

        # Merge rows with the same gene symbol
        df = df.groupby(df.index).mean().round(0).astype(int)

        df.to_csv("${meta.id}.symbols.tsv", sep="\\t")
    """
}