#!/usr/bin/env python3

import pandas as pd
import argparse
import statistics as st
import scipy.stats as stats

parser = argparse.ArgumentParser(description='Mann Whitney U test')
parser.add_argument('-i', '--input', help='Input file', required=True)
parser.add_argument('-o', '--output', help='Output file', required=True)

args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=0, index_col=0).T

# Drop all columns which contain exclusively NA values
df = df.dropna(axis=1, how='all')

# Save whole content of the dataframe in a single, flattened list
background = df.values.flatten().tolist()
background_median = st.median(background)

def mann_whitney_u(background, foreground):
    _, p = stats.mannwhitneyu(background, foreground)
    return p

# Transform df to have the following columns: sum, mean, q95, q99, median, p-value
df['sum'] = df.sum(axis=1)
df['mean'] = df.mean(axis=1)
df['q95'] = df.quantile(0.95, axis=1)
df['q99'] = df.quantile(0.99, axis=1)
df['median'] = df.median(axis=1)
df['p-value'] = df.apply(lambda x: mann_whitney_u(background, x), axis=1)

df = df[['sum', 'mean', 'q95', 'q99', 'median', 'p-value']]
df = df[(df['median'] > background_median) & (df['p-value'] < 0.01)]

df.sort_values(by=['median'], ascending=False, inplace=True)

length = len(df.index)
df['rank'] = range(1, length + 1)
df['dcg'] = 1 - (df['rank'] - 1) / length

df.to_csv(args.output, sep='\t')
