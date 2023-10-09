process RANKING {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::mulled-v2-cd5249a47f81a81b2e7785172c240f12497f55b4==c5c6cff7c28d3260400f938602ee600b1acf0323-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cd5249a47f81a81b2e7785172c240f12497f55b4:c5c6cff7c28d3260400f938602ee600b1acf0323-0':
        'biocontainers/mulled-v2-cd5249a47f81a81b2e7785172c240f12497f55b4:c5c6cff7c28d3260400f938602ee600b1acf0323-0' }"

    input:
        tuple val(meta), path(tftg)
    
    output:
        tuple val(meta), path("${meta.id}.ranking.tsv")
    
    script:
        """
        #!/usr/bin/env python3

        import pandas as pd
        import argparse
        import statistics as st
        import scipy.stats as stats

        df = pd.read_csv("${tftg.name}", sep='\\t', header=0, index_col=0).T

        # Drop all columns which contain exclusively NA values
        df = df.dropna(axis=1, how='all')

        # Save whole content of the dataframe in a single, flattened list
        background = df.values.flatten().tolist()
        background_median = st.median(background)

        def mann_whitney_u(background, foreground):
            _, p = stats.mannwhitneyu(background, foreground)
            return p

        # Transform df to have the following columns: sum, mean, q95, q99, median, p-value
        df['sum'] = df.sum(axis=1)
        df['mean'] = df.mean(axis=1)
        df['q95'] = df.quantile(0.95, axis=1)
        df['q99'] = df.quantile(0.99, axis=1)
        df['median'] = df.median(axis=1)
        df['p-value'] = df.apply(lambda x: mann_whitney_u(background, x), axis=1)

        df = df[['sum', 'mean', 'q95', 'q99', 'median', 'p-value']]
        df = df[(df['median'] > background_median) & (df['p-value'] < 0.01)]

        df.sort_values(by=['median'], ascending=False, inplace=True)

        length = len(df.index)
        df['rank'] = range(1, length + 1)
        df['dcg'] = 1 - (df['rank'] - 1) / length

        df.to_csv("${meta.id}.ranking.tsv", sep='\\t')
        """
}