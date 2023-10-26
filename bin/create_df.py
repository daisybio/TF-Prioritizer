#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import argparse
import os

parser = argparse.ArgumentParser(description="Create dataframe from count matrix and metadata")
parser.add_argument("--counts", type=str, help="Path to counts")
parser.add_argument("--metadata", type=str, help="Path to metadata")
parser.add_argument("--output", type=str, help="Path to output file")
args = parser.parse_args()

metadata = pd.read_csv(args.metadata, index_col=0, header=0)

if not os.path.exists(args.counts):
    raise Exception("Counts input does not exist")

if os.path.isfile(args.counts):
    # Use only cols which are in the index of metadata
    counts = pd.read_csv(args.counts, index_col=0, sep="\t", header=None)

    # If counts has no columns, add index name
    if len(counts.columns) == 0:
        counts.index.name = "gene_id"
    else:
        # Set first row as column names
        counts.columns = counts.iloc[0]
        # Remove first row
        counts = counts.iloc[1:]
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

for index, row in metadata.iterrows():
    if row["file"]:
        sample_df = pd.read_csv(row["file"], header=None)
        sample_counts = sample_df[0].to_list()
        counts[index] = sample_counts

samples = metadata.index.to_list()

if not all([sample in counts.columns for sample in samples]):
    raise Exception("Not all samples are in the counts matrix")

counts.to_csv(args.output, sep="\t")
counts.index.to_series().to_csv("genes.txt", index=False, header=False)