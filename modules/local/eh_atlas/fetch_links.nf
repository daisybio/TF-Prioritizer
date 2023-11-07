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
    #!/usr/bin/env python

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

    tf_link = {os.path.splitext(os.path.basename(link))[0]: BED_BASE_URL + link for link in links}

    df = pd.DataFrame.from_dict(tf_link, orient='index', columns=['link'])

    df.index.name = 'tf'

    df.to_csv("eh_atlas.tsv", sep='\\t')
    """
}
