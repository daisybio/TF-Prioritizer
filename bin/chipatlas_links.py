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
parser.add_argument('-s', '--threshold', dest='threshold', type=int, help='Threshold', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

df = pd.read_csv(args.list, index_col=0)

tissues = [tissue.lower() for tissue in args.tissues]
tfs = [tf.upper() for tf in args.tfs]

df['cell_type_class'] = df['cell_type_class'].str.lower()
df['antigen'] = df['antigen'].str.upper()

df = df[(df['genome_assembly'] == args.genome) & 
        (df['cell_type_class'].isin(tissues)) & 
        (df['antigen'].isin(tfs)) & 
        (df['threshold'] == args.threshold)]

# Group by antigen, keep only the row with minimum threshold
df.index = df['antigen']

df = df[['cell_type_class', 'file_url']]

# Rename columns
df.columns = ['tissue', 'url']
# Set index name
df.index.name = 'tf'

df.to_csv(args.output, sep='\t', index=True, header=True)