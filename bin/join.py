#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file", nargs="+")
parser.add_argument("--output", help="output file")

args = parser.parse_args()

dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in args.input]

df = pd.concat(dfs, axis=1, join="outer", sort=True)

df.to_csv(args.output, sep="\t", index=True)