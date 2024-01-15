#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

# Define the command-line arguments
parser = argparse.ArgumentParser(description="Calculate statistics between two multiple files.")
parser.add_argument("-i", "--input", type=str, nargs='+', help="List of input file paths", required=True)
parser.add_argument("-o", "--output", type=str, help="Output file path", required=True)
parser.add_argument("-m", "--method", type=str, choices=["mean", "sum", "ratio"], default="mean", help="Calculation method (mean, sum, ratio)")
args = parser.parse_args()

# Check if input and output paths are provided
if not args.input or not args.output:
    parser.error("Input and output paths are required.")

# Read all input files into a list of dataframes
dfs = [pd.read_csv(file, sep='\t', index_col=0) for file in args.input]

if args.method == "sum":
    index_union = dfs[0].index
    for df in dfs[1:]:
        index_union = index_union.union(df.index)

    # Add NA values for missing rows
    dfs = [df.reindex(index_union) for df in dfs]
else:
    index_intersection = dfs[0].index
    for df in dfs[1:]:
        index_intersection = index_intersection.intersection(df.index)
    
    print(f"Number of rows in intersection: {len(index_intersection)}")
    # Keep row indices which are available in all dataframes
    dfs = [df.loc[index_intersection] for df in dfs]

# Check if all dataframes have the same dimensions
if not all(df.shape == dfs[0].shape for df in dfs):
    raise ValueError(f"The input files must have the same dimensions. Got: {[df.shape for df in dfs]}")

# Check if all dataframes have the same row names
if not all(df.index.equals(dfs[0].index) for df in dfs):
    raise ValueError("The input files must have the same row names.")

# Check if all dataframes have the same column names
if not all(df.columns.equals(dfs[0].columns) for df in dfs):
    raise ValueError("The input files must have the same column names.")

# Calculate the selected statistic
if args.method == "mean":
    result = sum(dfs) / len(dfs)
elif args.method == "sum":
    result = sum(dfs)
elif args.method == "ratio":
    if len(dfs) != 2:
        raise ValueError("The ratio method requires exactly two input files.")
    
    # Replace 0 values with minimal existing float value
    dfs[1] = dfs[1].replace(0, np.finfo(float).eps)

    result = dfs[0] / dfs[1]

    print(f"Number of rows before dropping NA or inf values: {len(result)}")

    # Drop rows with NA or inf values (requirement for DYNAMITE)
    result = result.replace([np.inf, -np.inf], np.nan).dropna()

    print(f"Number of rows after dropping NA or inf values: {len(result)}")

# Write the result to a file
result.to_csv(args.output, sep='\t', index=True, quoting=0)