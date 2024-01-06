#!/usr/bin/env python3

import pandas as pd
import argparse


parser = argparse.ArgumentParser(prog="GFF to ROSE-GFF",
                                    description="Takes a standard GFF as input and reformats it into a ROSE input GFF")

parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)

args = parser.parse_args()
path_input = args.input
path_output = args.output


gff = pd.read_csv(path_input,
                  sep = "\t",
                  names = ["seqname", "source", "feature1", "start", "end", "score", "strand", "dot", "feature2"],
                  index_col=False).drop(columns=["score", "feature1", "dot", "feature2"])

if not all(gff['seqname'].str.startswith('chr')):
    gff['seqname'] = ["chr" + str(chrom) for chrom in gff['seqname'].tolist()]
gff['source'] = ['enhancer_'+ str(num) for num in range(gff.shape[0])]
gff['id2'] = gff['source']
gff['empty1'] = ''
gff['empty2'] = ''
gff['empty3'] = ''
gff = gff[['seqname', 'source', 'empty1', 'start', 'end', 'empty2', 'strand', 'empty3', 'id2']]
gff


gff.to_csv(path_output, sep = "\t", header = False, index=False)