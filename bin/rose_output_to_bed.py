#!/usr/bin/env python3

import pandas as pd
import argparse

def convert_rose_gff_to_bed(gff_file_path, bed_file_path):
    # Define column names for the GFF file
    columns = ["chr", "name1", "null1", "start", "stop", "null2", "strand", "null3", "name2"]
    
    # Read the GFF file
    gff = pd.read_csv(gff_file_path, sep="\t", header=None, index_col=False, names=columns)

    # adjust for 1-based to 0-based coordinate
    gff["start"] = gff["start"] - 1
    
    # Select relevant columns for BED format and write to file (use only 'chr', 'start', 'stop' for minimal bed file)
    gff[["chr", "start", "stop"]].to_csv(bed_file_path, sep="\t", header=False, index=False)

def main():
    parser = argparse.ArgumentParser(description='Convert GFF file to BED format.')
    parser.add_argument('--gff', '-i', type=str, help='Input GFF file path')
    parser.add_argument('--bed', '-o', type=str, help='Output BED file path')
    args = parser.parse_args()

    convert_rose_gff_to_bed(args.gff, args.bed)

if __name__ == "__main__":
    main()
