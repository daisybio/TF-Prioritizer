import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("input", type=argparse.FileType("r"))
parser.add_argument("output", type=str)
parser.add_argument("hm", type=str)
parser.add_argument("group1", type=str)
parser.add_argument("group2", type=str)

args = parser.parse_args()

sns.set_context("notebook")
color = "#A6CEE3"
sns.set_context("talk")
sns.set(font_scale=2)
sns.set_style("whitegrid")
plt.figure(figsize=(26, 20))

df = pd.read_csv(args.input, sep="\t")

if not df.empty:
    sns.set(font_scale=3)
    if len(df) > 61:
        sns.set(font_scale=2)
    if len(df) > 81:
        sns.set(font_scale=1)
    sns.set_style("whitegrid")
    ax = sns.barplot(x="TF", y="value", data=df, color=color)
    ax.set_title(args.hm + "\n" + args.group1 + ":" + args.group2)
    ax.set_ylabel('Normalized feature value\n(regression coefficient)')
    ax.set_xlabel('Transcription Factor')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(args.output)
    plt.clf()
plt.figure(figsize=(26, 20))
