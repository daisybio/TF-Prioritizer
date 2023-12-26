#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import pandas as pd


# Creates a cellmarkfiletable which is needed as input for ChromHMM
parser = argparse.ArgumentParser(description = "Script to remove full paths of input file to fit into nextflow workflow")
parser.add_argument("--input", help = "Input directory", required = True, type = str)
parser.add_argument("--output", help = "path for output file", required = True, type = str)

args = parser.parse_args()

input = args.input
output = args.output

table = pd.read_csv(input, sep = "\t", names=["state", "assay", "bam", "control"])

table["bam"] = [os.path.basename(path) for path in table["bam"]]
table["control"] = [os.path.basename(path) for path in table["control"]]
table.to_csv(output, header=False, sep="\t", index=False)
