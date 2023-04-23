#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Combine multiple deseq files')

parser.add_argument('-d', '--deseq', dest='deseq', type=str, help='DESeq files', required=True, nargs='+')
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

combined = pd.DataFrame()

for deseq_file in args.deseq:
    pairing = deseq_file.split('/')[-1].split('.')[0]
    deseq = pd.read_csv(deseq_file, sep=' ', index_col=0)
    
    log2FoldChange = deseq['log2FoldChange']
    
    log2FoldChange.name = pairing

    if not combined.empty and not combined.index.equals(log2FoldChange.index):
        raise Exception('DESeq files do not have the same genes')

    combined = pd.concat([combined, log2FoldChange], axis=1)

combined.to_csv(args.output, sep='\t')