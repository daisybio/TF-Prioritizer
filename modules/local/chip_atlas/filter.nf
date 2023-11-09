process FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6==fccb0c41a243c639e11dd1be7b74f563e624fcca-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0':
        'biocontainers/mulled-v2-2076f4a3fb468a04063c9e6b7747a630abb457f6:fccb0c41a243c639e11dd1be7b74f563e624fcca-0' }"
    
    input:
        tuple val(meta), path(sheet)
        val(genome)
        val(tissues)

    output:
        tuple val(meta), path("data.tsv")

    script:
    """
        #!/usr/bin/env python

        import pandas as pd
        import argparse

        cell_types = ["${tissues.join("', '")}"]

        df = pd.read_csv("${sheet}")

        # Filter by genome assembly
        df = df[((df["genome_assembly"] == "mm10") & (df["cell_type_class"].isin(cell_types)) & (df["cell_type"].isna()))]

        df = df[["antigen", "antigen_class", "cell_type_class", "threshold", "file_url"]]
        df["threshold"] = df["threshold"].fillna(0)
        df = df.fillna("")

        # Group by antigen, cell type class, cell type
        # Keep entry with lowest threshold
        df = df.sort_values("threshold")
        df = df.groupby(["antigen", "antigen_class", "cell_type_class"], as_index=False).first().reset_index()

        # Remove threshold column
        df.drop("threshold", axis=1, inplace=True)
        df.drop("index", axis=1, inplace=True)

        df.to_csv("data.tsv", sep="\\t", index=False)
    """
}