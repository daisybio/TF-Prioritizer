#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Calculate the TF-TG scores')
parser.add_argument('-l', '--log2fc', dest='log2fc', type=str, help='Log2FC file', required=True)
parser.add_argument('-a', '--affinities', dest='affinities', type=str, help='Affinities file', required=True)
parser.add_argument('-c', '--coefficients', dest='coefficients', type=str, help='Coefficients file', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

log2fc = pd.read_csv(args.log2fc, sep=' ', index_col=0)
affinities = pd.read_csv(args.affinities, sep='\t', index_col=0)
coefficients = pd.read_csv(args.coefficients, sep='\t', index_col=0)

# Restructure the affinities df so that its row names match the log2fc df index
genes = list(set(log2fc.index).intersection(set(affinities.index)))
genes.sort()
affinities = affinities.loc[genes]
log2fc = log2fc.loc[genes]

# Restructure the affinities df so that its column names match the coefficients df index
tfs = list(set(affinities.columns).intersection(set(coefficients.index)))
tfs.sort()
affinities = affinities[tfs]
coefficients = coefficients.loc[tfs]

# Calculate the TF-TG scores

## Multiply the log2FC by the affinities
result = affinities.mul(abs(log2fc["log2FoldChange"]), axis=0)

## Multiply the result by the coefficients
result = result.mul(abs(coefficients["value"]), axis=1)

# Save the result
result.to_csv(args.output, sep='\t')