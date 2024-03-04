#!/usr/bin/env python3

import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file")
parser.add_argument("--mapping", help="mapping file", default=None)
parser.add_argument("--output", help="output file")

args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", index_col=0)

# Remove additional metrics produced by STARE and keep only gene-TF affinities
df = df.drop(columns=["NumPeaks", "AvgPeakDistance", "AvgPeakSize"])

# Some column title are structured like "MEIS1(MA0498.1)" in this case we
# want to remove the annotation in parenthesis and keep only the TF name
# while calculating the mean of the expression values.

# Remove the annotation in parenthesis by using a regular expression
df = df.rename(columns=lambda x: re.sub(r"\(.*\)", "", x))

# Remove gene versions from row names
df = df.rename(index=lambda x: re.sub(r"\.\d+$", "", x))

if args.mapping is not None:
    # Map the gene ids to the gene names
    mapping_df = pd.read_csv(args.mapping, sep="\t", index_col=0)
    df = df.rename(index=mapping_df["gene_name"])

# Use max affinity for each TF
df = df.groupby(df.columns, axis=1).max()

# Use max for each gene
df = df.groupby(df.index).mean()

df.to_csv(args.output, sep="\t", index=True)