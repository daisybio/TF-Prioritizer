process FETCH_LINKS {
    container 'tfprio-python'

    input:
        val(tissue_types)
        val(tfs)
        val(genome)
        val(threshold)

    output:
        path('links.tsv')

    script:
        """
        if [ ! -f chip_atlas_file_list.csv ]; then
            curl -o overview.zip http://togodb.biosciencedbc.jp/togodb/release/chip_atlas_file_list.csv
            unzip overview.zip
        fi

        chipatlas_links.py -l chip_atlas_file_list.csv \\
            -t ${tissue_types.join(' ')} \\
            -f ${tfs.join(' ')} \\
            -g ${genome} \\
            -s ${threshold} \\
            -o links.tsv
        """
}

process FETCH_BED {
    container 'tfprio-python'

    input:
        tuple val(meta), val(url)

    output:
        tuple val(meta), path("${meta.id}.bed")

    script:
        """
        #!/usr/bin/env python

        import requests

        r = requests.get('${url}')
        with open('${meta.id}.bed', 'w') as f:
            f.write(r.text)
        """
}