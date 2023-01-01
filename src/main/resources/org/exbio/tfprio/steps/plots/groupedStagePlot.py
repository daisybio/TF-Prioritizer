import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("input", type=argparse.FileType("r"))
parser.add_argument("output", type=str)

args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t", index_col=0)

if not df.empty:
    change_colnames = True
    sns.set(font_scale=2)
    if len(df) > 41:
        sns.set(font_scale=1)
        change_colnames = False
    if len(df) > 61:
        sns.set(font_scale=0.8)
        change_colnames = False
    sns.set_style("whitegrid")
    plot = sns.heatmap(df.transpose(), cmap="Paired", square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5,
                       linecolor='black', xticklabels=True)
    plot.set_xlabel('Transcription Factor')
    plt.tight_layout()
    plt.savefig(args.output)
