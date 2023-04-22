#!/usr/bin/env python

import pandas as pd
import argparse
from typing import List
import logomaker
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Extract the biophysical models')

parser.add_argument('-t', '--tfs', dest='tfs', type=str, help='TFs file', required=True)
parser.add_argument('-p', '--pwms', dest='pwms', type=argparse.FileType(), help='PWMs file', required=True)

args = parser.parse_args()

# Read list of TFs
tfs = pd.read_csv(args.tfs, sep='\t', index_col=0).index.tolist()

blocks: List[str] = [block for block in args.pwms.read().split('>') if not block.startswith('#')]

for block in blocks:
    lines = block.splitlines()
    tf = lines[0].split('\t')[1]

    if tf not in tfs:
        continue
    
    # Create dataframe from block content
    df = pd.DataFrame([line.split() for line in lines[1:]], columns=['A', 'C', 'G', 'T']).apply(pd.to_numeric)

    df.index = range(1, len(df) + 1)
    df.to_csv(tf + '.pwm', sep='\t')

    # Create logo
    logo = logomaker.Logo(df, color_scheme='classic')

    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    logo.ax.set_ylabel('Binding energy', labelpad=-1)
    logo.ax.xaxis.set_ticks_position('none')
    logo.ax.xaxis.set_tick_params(pad=-1)

    plt.savefig(tf + '.png', bbox_inches='tight', dpi=300)
    plt.close()
