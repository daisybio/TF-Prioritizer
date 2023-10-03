process FETCH_SYMBOL_MAP {
    tag "$meta.id"
    label "process_low"

    conda "bioconda::mygene==3.2.2--pyh5e36f6f_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mygene:3.0.0--py27_0' :
        'biocontainers/mygene:3.0.0--py_2' }"

    input:
        tuple val(meta), path(genes)
        val(tax_id)

    output:
        tuple val(meta), path("*.json")

    script:
    """
        #!/usr/bin/env python3
        import mygene as mg
        import json

        with open("${genes}") as f:
            genes = f.read().splitlines()

        genes_filtered = [gene for gene in genes if gene.startswith("ENS")]

        mg = mg.MyGeneInfo()
        res = mg.querymany(genes_filtered, scopes='ensembl.gene', fields='symbol', species="${tax_id}")

        formatted_res = {x['query']: x['symbol'] for x in res if 'symbol' in x}

        with open("genes.json", "w+") as f:
            json.dump(formatted_res, f, indent=4)
    """
}