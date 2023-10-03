#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse
import os

parser = argparse.ArgumentParser(description="Create anndata object from count matrix and metadata")
parser.add_argument("--counts", type=str, help="Path to counts")
parser.add_argument("--metadata", type=str, help="Path to metadata")
parser.add_argument("--output", type=str, help="Path to output anndata object")
parser.add_argument("--index_col", type=str, help="Column to use as index in counts")
args = parser.parse_args()

metadata = pd.read_csv(args.metadata, index_col=0)

if not os.path.exists(args.counts):
    raise Exception("Counts input does not exist")

if os.path.isfile(args.counts):
    # Use only cols which are in the index of metadata
    counts = pd.read_csv(args.counts, index_col=args.index_col, usecols=[*metadata.index, args.index_col], sep="\t")
else:
    counts = pd.DataFrame(columns=metadata.index)
    for sub_name in os.listdir(args.counts):
        sub_path = os.path.join(args.counts, sub_name)

        if os.path.isfile(sub_path):
            with open(sub_path, "r") as f:
                counts.index = f.read().splitlines()
        else:
            for count_file in os.listdir(sub_path):
                sample_name = count_file.split(".")[0]

                if sample_name not in metadata.index:
                    print(f"Sample {sample_name} not in metadata, skipping")
                    continue

                sample_path = os.path.join(sub_path, count_file)
                with open(sample_path, "r") as f:
                    counts[sample_name] = f.read().splitlines()

                # Set dtype to int
                counts[sample_name] = counts[sample_name].astype(int)

# Create anndata object, rows are genes, cols are samples, values are raw counts
adata = ad.AnnData(X=counts.T, obs=metadata)

# Add TPMs as a layer
adata.layers["TPM"] = counts.T.div(counts.T.sum(axis=1), axis=0) * 1e6
adata.layers["raw_counts"] = counts.T

# Save anndata object
adata.write(args.output)

adata.var.index.to_series().to_csv("genes.txt", index=False, header=False)