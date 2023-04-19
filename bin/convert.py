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

df_counts["ENSG"] = df_counts.index.map(lambda x: ensg_dict[x] if x in ensg_dict else x)

# Put ENSG column first
cols = df_counts.columns.tolist()
cols = cols[-1:] + cols[:-1]
df_counts = df_counts[cols]

df_counts.to_csv(args.output, sep="\t", index=False)