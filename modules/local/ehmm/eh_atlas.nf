process SCRAPE_EH_ATLAS {
    container 'tfprio-python'

    input:
        val(tissues)
        val(genome)

    output:
        path('eh_atlas.tsv')

    script:
        """
        scrape_eh_atlas.py -t ${tissues.join(' ')} -g ${genome} -o eh_atlas.tsv
        """
}

process DOWNLOAD_EH_ATLAS {
    container 'tfprio-python'

    input:
        tuple val(meta), val(link)

    output:
        tuple val(meta), path("${meta.id}.bed")

    script:
        """
        #!/usr/bin/env python

        import requests

        r = requests.get('${link}')
        with open('${meta.id}.bed', 'w+') as f:
            f.write(r.text)
        """
}