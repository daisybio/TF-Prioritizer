#!/usr/bin/env python

import pandas as pd
import argparse
import json

parser = argparse.ArgumentParser(description='Find expression data of the transcription factors')

parser.add_argument('-t', '--tpm', dest='tpm', type=str, help='TPM file', required=True)
parser.add_argument('-c', '--counts', dest='counts', type=str, help='Counts file', default=None)
parser.add_argument('-d', '--deseq', dest='deseq', type=str, help='DESeq file', default=None)
parser.add_argument('-e', '--ensgs', dest='ensgs', type=str, help='ENSG IDs to process', required=True, nargs='+')

args = parser.parse_args()

# Read inputs
tpm = pd.read_csv(args.tpm, sep='\t', index_col=0)
counts = pd.read_csv(args.counts, sep='\t', index_col=0)
deseq = pd.read_csv(args.deseq, sep='\t', index_col=0)

# Build multiindex for each dataframe
tpm.columns = pd.MultiIndex.from_product([['TPM'], tpm.columns])
counts.columns = pd.MultiIndex.from_product([['Counts'], counts.columns])
deseq.columns = pd.MultiIndex.from_product([['DESeq'], deseq.columns])

# Combine all dataframes while making sure the index is the same
combined = pd.concat([tpm, counts, deseq], axis=1, join='outer')

# Build intersection of existing and wanted ENSG IDs
ensgs = list(set(combined.index).intersection(args.ensgs))
ensgs.sort()

# Filter out the ENSG IDs we want
combined = combined.loc[ensgs]

# Write output, json pretty-printed, separate file for each ENSG ID
# First level of the json is the ENSG ID, second level is the multiindex, third level is the sample

for ensg in ensgs:
    ensg_dict = {
        'TPM': combined.loc[ensg, 'TPM'].to_dict(),
        'Counts': combined.loc[ensg, 'Counts'].to_dict(),
        'DESeq': combined.loc[ensg, 'DESeq'].to_dict()
    }
    with open(ensg + '.json', 'w') as f:
        json.dump(ensg_dict, f, indent=4)