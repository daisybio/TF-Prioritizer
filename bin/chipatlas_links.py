#!/usr/bin/env python

import pandas as pd
import argparse
import os
import requests

parser = argparse.ArgumentParser(description='Fetch binding data of the transcription factors')

parser.add_argument('-l', '--list', dest='list', type=str, help='List file', required=True)
parser.add_argument('-g', '--genome', dest='genome', type=str, help='Genome', required=True)
parser.add_argument('-t', '--tissues', dest='tissues', type=str, help='Tissues to process', required=True, nargs='+')
parser.add_argument('-a', '--antigenes', dest='antigenes', type=str, help='Antigenes to process', required=True, nargs='+')
parser.add_argument('-s', '--threshold', dest='threshold', type=int, help='Threshold', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

df = pd.read_csv(args.list, index_col=0)
for col in df.columns:
        if col == 'threshold':
                continue
        df[col] = df[col].fillna('')
        df[col] = df[col].astype(str)

tissues = [tissue.lower() for tissue in args.tissues]
antigenes = args.antigenes

df['cell_type_class'] = df['cell_type_class'].str.lower()
df['key'] = df['antigen_class'] + '_' + df['antigen']
df['key'] = df['key'].str.rstrip('_')

df = df[(df['genome_assembly'] == args.genome) & 
        (df['cell_type_class'].isin(tissues)) & 
        (df['key'].isin(antigenes)) & 
        (df['cell_type'] == '') &
        (df['threshold'] == args.threshold)]

# Replace antigen with ALL if it is empty
df['antigen'] = df['antigen'].apply(lambda x: 'ALL' if x == '' else x)

df['id'] = df['antigen'] + "_" + str(args.threshold) + "_" + df['cell_type_class'] + "_" + df['antigen_class']

df = df[['file_url', 'id']]

df.to_csv(args.output, sep='\t', index=False, header=True)