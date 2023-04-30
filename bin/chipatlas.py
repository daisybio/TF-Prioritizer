#!/usr/bin/env python

import pandas as pd
import argparse
import os
import requests

parser = argparse.ArgumentParser(description='Fetch binding data of the transcription factors')

parser.add_argument('-l', '--list', dest='list', type=str, help='List file', required=True)
parser.add_argument('-g', '--genome', dest='genome', type=str, help='Genome', required=True)
parser.add_argument('-t', '--tissues', dest='tissues', type=str, help='Tissues to process', required=True, nargs='+')
parser.add_argument('-f', '--tfs', dest='tfs', type=str, help='TFs to process', required=True, nargs='+')

args = parser.parse_args()

df = pd.read_csv(args.list, index_col=0)

tissues = [tissue.lower() for tissue in args.tissues]
tfs = [tf.upper() for tf in args.tfs]

df['cell_type_class'] = df['cell_type_class'].str.lower()
df['antigen'] = df['antigen'].str.upper()

df = df[(df['genome_assembly'] == args.genome) & (df['cell_type_class'].isin(tissues)) & (df['antigen'].isin(tfs))]

# Group by antigen, keep only the row with minimum threshold
df = df.sort_values(by=['antigen', 'threshold']).groupby('antigen').first()

df = df[['cell_type_class', 'file_url']]

for index, row in df.iterrows():
    bed_file = f"{index}_{row['cell_type_class']}.bed"

    if os.path.exists(bed_file):
        print(f"{bed_file} already exists, skipping...")
        continue

    print(f"Downloading {bed_file} from {row['file_url']}...")

    r = requests.get(row['file_url'], allow_redirects=True)
    
    with open(bed_file, 'wb') as f:
        f.write(r.content)