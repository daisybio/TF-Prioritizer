#!/usr/bin/env python

import argparse
import bs4
import requests
import re
import os
import pandas as pd

URL = "http://www.enhanceratlas.org/downloadv2.php"
BED_BASE_URL = "http://www.enhanceratlas.org/"

parser = argparse.ArgumentParser(description='Scrape EH Atlas')

parser.add_argument('-g', '--genome', dest='genome', help='Genome', required=True)
parser.add_argument('-t', '--tissues', dest='tissues', help='Tissues', required=True, nargs='+')
parser.add_argument('-o', '--output', dest='output', help='Output file', required=True)

args = parser.parse_args()

genomeToKeyMap = {
    "hg38": "hs",
    "GRCh38": "hs",
    "hg19": "hs",
    "GRCh37": "hs",
    "mm10": "mm",
    "GRCm38": "mm",
    "mm9": "mm",
    "GRCm37": "mm",
    "zv10": "dr",
    "danRer10": "dr",
    "GRCz10": "dr",
    "zv9": "dr",
    "GRCz9": "dr",
    "dm6": "dm",
    "BDGP6": "dm",
    "ce10": "ce",
    "WBcel235": "ce",
    "rn6": "rn",
    "Rnor_6.0": "rn",
    "Galgal4": "gg",
    "galGal4": "gg",
    "susScr3": "ss",
    "Sscrofa10.2": "ss",
    "sacCer3": "sc",
    "R64-1-1": "sc",
}

enhancerVersionMap = {
    "hs": "hg19",
    "mm": "mm9",
    "dr": "danRer10",
    "dm": "dm3",
    "ce": "ce10",
    "rn": "rn5",
    "sc": "sacCer3",
    "gg": "galGal4",
    "ss": "susScr3"
}


# Get the key for the genome
if args.genome not in genomeToKeyMap:
    raise Exception(f"Genome {args.genome} not found in genomeToKeyMap")

genomeKey = genomeToKeyMap[args.genome]

# Fetch the webpage from the URL
page = requests.get(URL)

# Create a BeautifulSoup object
soup = bs4.BeautifulSoup(page.text, 'html.parser')

link_regex = f"data/download/enhancer/{genomeKey}/.+\.bed"

# Find all the links on the page that match the regex
link_objects = soup.find_all('a', href=re.compile(link_regex))
links = [link.get('href') for link in link_objects]

tf_link = {os.path.splitext(os.path.basename(link))[0]: BED_BASE_URL + link for link in links}

df = pd.DataFrame.from_dict(tf_link, orient='index', columns=['link'])
df['genome'] = enhancerVersionMap[genomeKey]

df.index.name = 'tf'

df.to_csv(args.output, sep='\t')