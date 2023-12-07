#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description="Process ChromHMM output into bed file of predicted enhancers")

parser.add_argument("-e", "--emissions", type=str, required=True, help="Path to emission file")
parser.add_argument("-b", "--bed", type=str, required=True, help="Path to bed file")
parser.add_argument("-t", "--threshold", type=float, required=False, default=0.9, help="Threshold for state emissions")
parser.add_argument("-m", "--markers", nargs='+', required=False, default=["H3K27ac", "H3K4me3"], help="ChIP-Seq markers that indicate an enhancer")
parser.add_argument("-o", "--output", type=str, required=True, help="Path to output bed with enhancer positions")

args = parser.parse_args()

path_emissions = args.emissions
path_bed = args.bed
threshold = args.threshold
markers = args.markers
output = args.output


# Read emissions file for the provided markers
emissions = pd.read_csv(path_emissions, sep = "\t")[["State (Emission order)"] + markers].rename(columns={"State (Emission order)": "State"})


# Read input bed file and remove unecessary columns
bed = pd.read_csv(path_bed,
                  sep="\t",
                  skiprows=1,
                  names=["chr", "start", "end", "state", "score", "strand", "start_1", "end_1", "rgb"]
                 ).drop(columns=["strand", "score", "start_1", "end_1", "rgb"])


# Keep state if any of the markers is enriched > threshold for this state
states = emissions[np.any([emissions[marker] >= threshold for marker in markers], axis=0)]["State"].tolist()


# Subset bed file for selected states
out_bed = bed[np.isin(bed["state"], states)].drop(columns=["state"])

# Write output
out_bed.to_csv(output, index=False, sep="\t", header=False)

