import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as scp
import seaborn as sns
import statistics as sts

parser = argparse.ArgumentParser()
parser.add_argument('--backgroundFile', type=str, required=True)
parser.add_argument('--statsFile', type=str, required=True)
parser.add_argument('--inputDirectory', type=str, required=True)
parser.add_argument('--outputDirectory', type=str, required=True)
args = parser.parse_args()

sns.set_context("notebook")
color = "#A6CEE3"
sns.set_context("talk")
sns.set(font_scale=4)
sns.set_style("whitegrid")

plt.figure(figsize=(20, 17))

background = pd.read_table(
    args.backgroundFile,
    comment="#", usecols=['TfTgScore']).sort_values(['TfTgScore'], ascending=False)
background["label"] = "background"
background.dropna(subset=["TfTgScore"], inplace=True)

background_sum = sum(background["TfTgScore"])
background_length = len(background)
background_mean = background_sum / background_length
background_median = sts.median(background['TfTgScore'])
background_quantile95 = np.percentile(background["TfTgScore"], 95)
background_quantile_99 = np.percentile(background["TfTgScore"], 99)

df_interesting_stats = pd.DataFrame(
    columns=['label', 'sum', 'count', 'mean', 'median', '95_quantile', '99_quantile'], dtype=float)
df_interesting_stats['label'] = df_interesting_stats['label'].astype(str)

df_interesting_stats.loc[0] = ['background', background_sum, background_length, background_mean, background_median,
                               background_quantile95, background_quantile_99]


def generate(input_file: str, plot_file: str, tf_name: str):
    df = pd.read_table(
        input_file,
        comment="#", usecols=['TfTgScore', 'TF']).sort_values(['TfTgScore'], ascending=False)
    df.columns = ['TfTgScore', 'label']
    df.dropna(subset=["TfTgScore"], inplace=True)

    length = len(df)
    if length < 1:
        return

    df_sum = df['TfTgScore'].sum()
    mean = df["TfTgScore"].mean()
    quantile99 = np.percentile(df['TfTgScore'], 99)
    quantile95 = np.percentile(df['TfTgScore'], 95)
    median = sts.median(df['TfTgScore'])
    mann_whitney_u = scp.mannwhitneyu(background['TfTgScore'], df['TfTgScore'])

    if median > background_median and mann_whitney_u[1] < 0.01:
        tf_background = pd.concat([background, df], axis=0)

        ax = sns.boxplot(x="label", y="TfTgScore", data=tf_background, palette="husl")
        ax.set_yscale("log")
        ax.set(xlabel='', ylabel='TF-TG score')
        plt.savefig(plot_file)
        plt.clf()

        return [tf_name, df_sum, length, mean, median, quantile95, quantile99]
    return None


for tf_file in os.listdir(args.inputDirectory):
    name = tf_file.rstrip(".tsv")
    in_file = os.path.join(args.inputDirectory, tf_file)
    pl_file = os.path.join(args.outputDirectory, name + ".png")
    result = generate(in_file, pl_file, name)
    if result is not None:
        df_interesting_stats.loc[len(df_interesting_stats)] = result

df_interesting_stats.sort_values("median", ascending=False, inplace=True)
df_interesting_stats.reset_index(inplace=True, drop=True)
df_interesting_stats.to_csv(args.statsFile, sep='\t')
