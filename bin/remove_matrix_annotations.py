#!/usr/bin/env python3

import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file")
parser.add_argument("--output", help="output file")

args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", index_col=0)

# Some column title are structured like "MEIS1(MA0498.1)" in this case we
# want to remove the annotation in parenthesis and keep only the TF name
# while calculating the mean of the expression values.

# Remove the annotation in parenthesis by using a regular expression
df = df.rename(columns=lambda x: re.sub(r"\(.*\)", "", x))

# Calculate the mean of the expression values for each TF
df = df.groupby(df.columns, axis=1).mean()

df.to_csv(args.output, sep="\t", index=True)