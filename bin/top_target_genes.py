#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Find the top target genes')
parser.add_argument('-t', '--tfs', dest='tfs', type=str, help='TFs file', required=True)
parser.add_argument('-a', '--affinities', dest='affinities', type=str, help='Affinities file', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

# Read list of TFs
tfs = pd.read_csv(args.tfs, sep='\t', index_col=0).index.tolist()

# Read affinities
affinities = pd.read_csv(args.affinities, sep='\t', index_col=0, usecols=tfs)

# Find the top target genes by sorting the genes based on the affinites
top_target_genes = affinities.apply(lambda x: x.sort_values(ascending=False).index.tolist(), axis=0)

# Save the result
top_target_genes.to_csv(args.output, sep='\t', index=False)
