import scipy.stats as scp
import statistics as sts
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set_context("notebook")
color = "#A6CEE3"
sns.set_context("talk")
sns.set(font_scale=4)
sns.set_style("whitegrid")

plt.figure(figsize=(20, 17))

df_interesting_stats = pd.DataFrame(
    columns=['label', 'sum_all_values', 'number_target_genes', 'mean', 'median', '95_quantile', '99_quantile'])

background = pd.read_table(
    '{BACKGROUNDFILE}',
    comment="#", usecols=['TF_TG_SCORE']).sort_values(['TF_TG_SCORE'], ascending=False)
background["label"] = "background"
background.dropna(subset=["TF_TG_SCORE"], inplace=True)

background_sum = sum(background["TF_TG_SCORE"])
background_length = len(background)
background_mean = background_sum / background_length
background_median = sts.median(background['TF_TG_SCORE'])
background_quantile = np.percentile(background["TF_TG_SCORE"], 95)
background_quantile_99 = np.percentile(background["TF_TG_SCORE"], 99)

row_counter = 1


def generate(inputFile: str, plotFile: str, tfName: str):
    global row_counter
    df = pd.read_table(
        inputFile,
        comment="#", usecols=['TF_TG_SCORE', 'TF']).sort_values(['TF_TG_SCORE'], ascending=False)
    df.columns = ['TF_TG_SCORE', 'label']
    df.dropna(subset=["TF_TG_SCORE"], inplace=True)
    tfSum = sum(df['TF_TG_SCORE'])
    length = len(df)
    mean = 0
    quantile = 0
    quantile95 = 0
    median = 0
    mannWhitneyU = scp.mannwhitneyu(
        background['TF_TG_SCORE'], df['TF_TG_SCORE'])

    if (length > 0):
        mean = tfSum / length
        quantile = np.percentile(df['TF_TG_SCORE'], 99)
        quantile95 = np.percentile(df['TF_TG_SCORE'], 95)
        median = sts.median(df['TF_TG_SCORE'])

    if (median > background_median and mannWhitneyU[1] < 0.01):
        tfBackground = pd.concat([background, df], axis=0)
        ax = sns.boxplot(x="label", y="TF_TG_SCORE",
                         data=tfBackground, palette="husl")
        ax.set_yscale("log")
        ax.set(xlabel='', ylabel='TF-TG score')
        plt.savefig(plotFile)
        df_interesting_stats.loc[row_counter] = [tfName, tfSum, length, mean, median,
                                                 quantile95, quantile]
        row_counter += 1

    plt.clf()


{CALLS}

df_interesting_stats.sort_values("median", ascending=False, inplace=True)

df_interesting_stats.loc[0] = ['background', background_sum, background_length, background_mean, background_median,
                               background_quantile, background_quantile_99]
df_interesting_stats.sort_index(inplace=True)
df_interesting_stats.to_csv(
    '{STATSFILE}',
    sep='\t')
