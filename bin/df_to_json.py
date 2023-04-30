#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Convert a CSV to JSON')
parser.add_argument('-i', '--input', dest='input',
                    help='Input CSV', required=True)
parser.add_argument('-o', '--output', dest='output', type=argparse.FileType('w'),
                    help='Output JSON', required=True)

args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', index_col=0,
                 header=None)

args.output.write(df[1].to_json(indent=4))
