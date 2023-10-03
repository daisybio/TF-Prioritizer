process NORMALIZE {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::pydeseq2==0.4.0--pyhdfd78af_0"
    container "biocontainers/pydeseq2:0.4.0--pyhdfd78af_0"

    input:
        tuple val(meta), path(anndata)

    output:
        tuple val(meta), path("*.normalized.h5ad")

    script:
    """
        #!/usr/bin/env python3
        import anndata as ad
        from pydeseq2.dds import DeseqDataSet
        import pandas as pd
        import json

        adata = ad.read_h5ad("${anndata}")

        design_factors = ['state', 'batch'] if 'batch' in adata.obs.columns and adata.obs['batch'].nunique() > 1 else ['state']

        dataset = DeseqDataSet(adata=adata, design_factors=design_factors, n_cpus=${task.cpus})
        dataset.deseq2()

        adata.layers['normalized'] = dataset.layers['normed_counts']

        adata.write_h5ad("${meta.id}.normalized.h5ad")
    """
}