process DESEQ2 {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::pydeseq2==0.4.0--pyhdfd78af_0"
    container "biocontainers/pydeseq2:0.4.0--pyhdfd78af_0"

    input:
        tuple val(meta), path(anndata)

    output:
        tuple val(meta), path("*.tsv")

    script:
    """
        #!/usr/bin/env python3
        import anndata as ad
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        import json

        adata = ad.read_h5ad("${anndata}")

        design_factors = ['state', 'batch'] if 'batch' in adata.obs.columns and adata.obs['batch'].nunique() > 1 else ['state']

        dataset = DeseqDataSet(adata=adata, design_factors=design_factors, n_cpus=${task.cpus})
        dataset.deseq2()

        states = adata.obs['state'].unique()

        for state1 in states:
            for state2 in states:
                if state1 >= state2:
                    continue

                stats = DeseqStats(dds=dataset, n_cpus=${task.cpus}, contrast=['state', state1, state2])
                stats.summary()

                stats.results_df.to_csv(f"{state1}:{state2}.tsv", sep="\\t")
    """
}