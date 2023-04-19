#!/usr/bin/env python3

import pandas as pd
import mygene
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file")
parser.add_argument("--output", help="output file")
parser.add_argument("--taxonomy", help="species taxonomy as found in https://docs.mygene.info/en/latest/doc/data.html#species")

args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", header=None)

symbols = df[0].tolist()
mg = mygene.MyGeneInfo()
res = mg.querymany(symbols, scopes="symbol", fields="ensembl.gene", species=args.taxonomy, as_dataframe=True)

mapping = res["ensembl.gene"]
notna = mapping.dropna()

found_percentage = len(notna) / len(mapping)

if found_percentage < 0.5:
    print("WARNING: Only {0:.2f}% of genes were found".format(found_percentage * 100))
else:
    print("Found {0:.2f}% of genes".format(found_percentage * 100))

notna.to_csv(args.output, sep="\t", header=False)