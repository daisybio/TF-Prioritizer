#!/usr/bin/env python3

import pandas as pd
import mygene
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", help="input file")
parser.add_argument("--gtf", help="gtf file")
parser.add_argument("--output", help="output file")

args = parser.parse_args()

# Proccess gtf file 
gtf_mapping = pd.read_csv('/nfs/data/COM2POSE/reference_data/gencode.vM25.annotation.gtf', sep='\t', header=None, comment='#')
gtf_mapping = gtf_mapping[gtf_mapping[2] == 'gene']
## Split the 9th column by semicolon into distinct columns
gtf_mapping = gtf_mapping[8].str.split(';', expand=True)[[0, 2]]
gtf_mapping.columns = ['gene_id', 'gene_name']
gtf_mapping["gene_id"] = gtf_mapping["gene_id"].str.lstrip('gene_id ').str.strip('"')
gtf_mapping["gene_name"] = gtf_mapping["gene_name"].str.lstrip('gene_name ').str.strip('"')
gtf_mapping.index = gtf_mapping["gene_name"]
gtf_mapping.drop("gene_name", axis=1, inplace=True)
gtf_mapping["gene_id"] = gtf_mapping["gene_id"].str.split('.', expand=True)[0]

## Group by index, select first gene id
gtf_mapping = gtf_mapping.groupby(level=0).first()

# Read input file

symbols = pd.read_csv(args.input, sep="\t", header=None)[0].tolist()

# Remove all symbols which are already in the gtf file
symbols = [s for s in symbols if s not in gtf_mapping.index]

mg = mygene.MyGeneInfo()
res = mg.querymany(symbols, scopes="symbol", fields="ensembl.gene", as_dataframe=True)


notna = res[["ensembl.gene"]].dropna()
notna.columns = ["gene_id"]

notna = notna[notna["gene_id"].str.startswith("ENSMUSG")]

# Concatenate the two mappings into a single series
mapping = pd.concat([gtf_mapping["gene_id"], notna["gene_id"]])

# Sort by index
mapping = mapping.sort_index()

# Remove duplicate indices
mapping = mapping[~mapping.index.duplicated(keep='first')]

mapping.to_csv(args.output, sep="\t", header=False)