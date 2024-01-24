process FETCH_LINKS {
    tag "fetch_eh_links"
    label "process_single"

    container 'registry.hub.docker.com/bigdatainbiomedicine/inspect-python'

    input:
        val(genome)

    output:
        path('eh_atlas.tsv')

    script:
    """
    #!/usr/bin/env python3

    import bs4
    import requests
    import re
    import os
    import pandas as pd

    URL = "http://www.enhanceratlas.org/downloadv2.php"
    BED_BASE_URL = "http://www.enhanceratlas.org/"

    # Fetch the webpage from the URL
    page = requests.get(URL)

    # Create a BeautifulSoup object
    soup = bs4.BeautifulSoup(page.text, 'html.parser')

    link_regex = f"data/download/enhancer/$genome/.+\\.bed"

    # Find all the links on the page that match the regex
    link_objects = soup.find_all('a', href=re.compile(link_regex))
    links = [link.get('href') for link in link_objects]

    tissue_link = {os.path.splitext(os.path.basename(link))[0]: BED_BASE_URL + link for link in links}

    df = pd.DataFrame.from_dict(tissue_link, orient='index', columns=['url'])

    df.index.name = 'tissue'

    df.to_csv("eh_atlas.tsv", sep='\\t')
    """
}
