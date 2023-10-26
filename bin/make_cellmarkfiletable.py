#!/usr/bin/env python
import os
import argparse
import pandas as pd

# Creates a cellmarkfiletable which is needed as input for ChromHMM

parser = argparse.ArgumentParser(description = "Script to convert directory with different stages and their ChIP-Seq BAM files to a cellmarkfiletable (ChromHMM)")
parser.add_argument("--input_dir", help = "Input directory", required = True, type = str)
parser.add_argument("--output", help = "path for output file", required = True, type = str)
parser.add_argument("--markers", help = "path to file with possible markers (default markers are always set)", required = True, type = str)
parser.add_argument("--control", help = "file name marker for the control files", default = "CONTROL", type = str, required = False)

args = parser.parse_args()

input_dir = args.input_dir
output = args.output
path_markers = args.markers
control = args.control


# Extend default markers manually if we find new markers in new data
markers = {"H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "MED1", "NFIB", "PolII", "STAT5", "CONTROL"}

# Read marker file
with open(path_markers, "r") as f:
    add_markers = f.read().splitlines()

# Remove lines in file that were empty
while "" in add_markers:
    add_markers.remove("")

markers = markers.union(set(add_markers))

# Get stages from directory names
stages = [os.path.join(input_dir, element) for element in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, element))]

output_list = []
for stage in stages:
    # Get control file for stage
    control_file = next((f for f in os.listdir(stage) if control in f and f.endswith(".bam")))
    if control_file is None:
        sys.exit(f"No control file found for {os.path.basename(stage)}")
    # Iterate over files in current stage
    for file in os.listdir(stage):
        # Skip bai files without message
        if file.endswith(".bai"):
            continue
        # Sort for bam files only
        if not file.endswith(".bam"):
            print(f"Skipping file {file}, not a bam file")
            continue
        marker = next((m for m in markers if m in file), None)
        
        # Print message if file couldn't be assigned to any marker
        if marker is None:
            print(f"Skipping file {file}, no supported marker found!")
            continue
    
        output_list.append((os.path.basename(stage), marker, file, control_file))     

output_file = ""
for stage, marker, file, control in output_list:
    output_file += f"{stage}\t{marker}\t{file}\t{control}\n"

with open(output, "w") as f:
    f.write(output_file)

