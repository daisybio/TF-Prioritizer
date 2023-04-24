#!/usr/bin/env python

import pandas as pd
import argparse
import json
import os

parser = argparse.ArgumentParser(description='Collect ranks of the transcription factors')

parser.add_argument('-r', '--ranks', dest='ranks', type=str, help='Rank files', required=True, nargs='+')
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

result = {}

for file in args.ranks:
    cleaned_name = os.path.splitext(os.path.basename(file))[0]
    pairing, hm, _ = cleaned_name.split('_')

    df = pd.read_csv(file, sep='\t', index_col=0)
    
    if not hm in result:
        result[hm] = {}
    
    result[hm][pairing] = df.index.tolist()

with open(args.output, 'w+') as f:
    json.dump(result, f, indent=4)
