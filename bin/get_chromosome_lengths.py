#!/usr/bin/env python

import json
import argparse
import requests
import pandas as pd

parser = argparse.ArgumentParser(description='Fetch chromosome lengths from biomaRt')

parser.add_argument('-o', '--output', dest='output', type=str, required=True)
parser.add_argument('-s', '--species', dest='species', default='hsapiens_gene_ensembl')

args = parser.parse_args()

species = args.species.split('_')[0]

res = requests.get(f"https://rest.ensembl.org/info/assembly/{species}?content-type=application/json")

if not res.ok:
    res.raise_for_status()

top_level = res.json()["top_level_region"]

dictionary = {
    entry["name"]: entry["length"] 
        for entry in top_level 
        if entry['coord_system'] == 'chromosome'
}

df = pd.DataFrame.from_dict(dictionary, orient='index')

df.sort_index(inplace=True)

df.to_csv(args.output, sep='\t', header=False)