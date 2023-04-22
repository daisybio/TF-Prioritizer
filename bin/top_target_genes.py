#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Find the top target genes')
parser.add_argument('-t', '--tfs', dest='tfs', type=str, help='TFs file', required=True)
parser.add_argument('-a', '--affinities', dest='affinities', type=str, help='Affinities file', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)
parser.add_argument('-k', '--top', dest='top', type=int, help='Number of top target genes', default=30)
parser.add_argument('-c', '--counts', dest='counts', type=str, help='Counts file', default=None)

args = parser.parse_args()

# Read list of TFs
tfs = pd.read_csv(args.tfs, sep='\t', index_col=0).index.tolist()

# Read affinities
affinities = pd.read_csv(args.affinities, sep='\t', index_col=0, usecols=tfs)

if args.counts is not None:
    # Read counts
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)

    # Build intersection of genes in affinities and counts
    intersection = affinities.index.intersection(counts.index)

    # Keep only affinities for genes which are present in the intersection
    affinities = affinities.loc[intersection]

# Find the top target genes by sorting the genes based on the affinites
top_target_genes = affinities.apply(lambda x: x.sort_values(ascending=False).index.tolist(), axis=0).head(args.top)

# Save the result
top_target_genes.to_csv(args.output, sep='\t', index=False)
