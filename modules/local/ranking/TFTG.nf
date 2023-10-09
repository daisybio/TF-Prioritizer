process TFTG_SCORE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"

    input:
        tuple val(meta), path(affinity_sum), path(differential_expression), path(dynamite)
  
    output:
        tuple val(meta), path("${meta.id}.tftg.tsv")
  
    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import argparse

        log2fc = pd.read_csv("${differential_expression.name}", sep='\\t', index_col=0)
        affinities = pd.read_csv("${affinity_sum.name}", sep='\\t', index_col=0)
        coefficients = pd.read_csv("${dynamite.name}", sep='\\t', index_col=0)

        # Restructure the affinities df so that its row names match the log2fc df index
        genes = list(set(log2fc.index).intersection(set(affinities.index)))
        genes.sort()

        if len(genes) == 0:
            print("No genes found")
            exit(1)

        affinities = affinities.loc[genes]
        log2fc = log2fc.loc[genes]

        # Restructure the affinities df so that its column names match the coefficients df index
        tfs = list(set(affinities.columns).intersection(set(coefficients.index)))
        tfs.sort()

        if len(tfs) == 0:
            print("No TFs found")
            exit(1)

        affinities = affinities[tfs]
        coefficients = coefficients.loc[tfs]

        # Calculate the TF-TG scores

        ## Multiply the log2FC by the affinities
        result = affinities.mul(abs(log2fc["log2FoldChange"]), axis=0)

        ## Multiply the result by the coefficients
        result = result.mul(abs(coefficients["value"]), axis=1)

        # Make sure results are not empty
        if result.empty:
            print("No results found")
            exit(1)

        # Save the result
        result.to_csv("${meta.id}.tftg.tsv", sep='\\t')
        """
}