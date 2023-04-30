#!/usr/bin/env python

import argparse
import requests
import os

DATA_URL = 'https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_redundant_pfms_jaspar.txt'
LOGO_URL = 'https://jaspar.genereg.net/static/logos/all/svg/'

parser = argparse.ArgumentParser(description='Fetch jaspar logos')

parser.add_argument('-t', '--tfs', dest='tfs', type=argparse.FileType(), help='TFs file', required=True)

args = parser.parse_args()

# Read list of TFs
tfs = args.tfs.read().splitlines()

content = requests.get(DATA_URL).text

entries = [entry[1:].split() for entry in content.splitlines() if entry.startswith('>')]

for entry in entries:
    matrix_id = entry[0]
    tf = entry[1].upper().replace('::', '..')

    if tf not in tfs:
        continue

    file_name = f"{tf}_jaspar/{matrix_id}.svg"
    os.makedirs(os.path.dirname(file_name), exist_ok=True)

    current_url = LOGO_URL + matrix_id + '.svg'

    with open(file_name, 'wb') as f:
        f.write(requests.get(current_url).content)