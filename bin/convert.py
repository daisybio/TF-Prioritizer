#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file")
parser.add_argument("--map", help="Map file")
parser.add_argument("--output", help="output file")

args = parser.parse_args()

df_counts = pd.read_csv(args.input, sep="\t", index_col=0)
df_map = pd.read_csv(args.map, sep="\t", header=None, index_col=0)

original_size = len(df_counts.index)

# Create intersection of genes
ensg_dict = df_map.to_dict()[1]

df_counts.index = df_counts.index.map(lambda x: ensg_dict[x] if x in ensg_dict else x)

# Group by ENSG, sum counts
df_counts = df_counts.groupby(df_counts.index).sum()

# Crash if there are duplicates in ENSG
assert len(df_counts.index) == len(set(df_counts.index))

df_counts.to_csv(args.output, sep="\t", index=True)