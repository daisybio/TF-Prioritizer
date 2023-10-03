#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse

parser = argparse.ArgumentParser(description="Create anndata object from count matrix and metadata")
parser.add_argument("--counts", type=str, help="Path to count matrix")
parser.add_argument("--metadata", type=str, help="Path to metadata")
parser.add_argument("--output", type=str, help="Path to output anndata object")
parser.add_argument("--index_col", type=str, help="Column to use as index in counts")
args = parser.parse_args()

metadata = pd.read_csv(args.metadata, index_col=0)

# Use only cols which are in the index of metadata
counts = pd.read_csv(args.counts, index_col=args.index_col, usecols=[*metadata.index, args.index_col], sep="\t")

# Create anndata object, rows are genes, cols are samples, values are raw counts
adata = ad.AnnData(X=counts.T, obs=metadata)

# Add TPMs as a layer
adata.layers["TPM"] = counts.T.div(counts.T.sum(axis=1), axis=0) * 1e6

# Save anndata object
adata.write(args.output)

adata.var.index.to_series().to_csv("genes.txt", index=False, header=False)