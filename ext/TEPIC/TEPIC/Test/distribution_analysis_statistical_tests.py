import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])

import io
from base64 import b64encode
import_or_install("plotly.express")
import plotly.express as px
import_or_install("dash")
import_or_install("dash_core_components")
import_or_install("dash_html_components")
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

import_or_install("pandas")
import_or_install("seaborn")
import_or_install("matplotlib.pyplot")
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_context("notebook")
color = "#A6CEE3"
sns.set_context("talk")
sns.set_style("whitegrid")

import_or_install("numpy")
import numpy as np
import_or_install("sts")
import statistics as sts
import scipy.stats as scp
plt.figure(figsize=(20, 17))


df_interesting_stats=pd.DataFrame(columns=['label','sum_all_values','number_target_genes','mean','median','95_quantile','99_quantile'])
row_counter=0

background=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/01_BACKGROUND_DISTRIBUTION/01_ALL/distribution.csv', comment="#", usecols=['TF_TG_SCORE']).sort_values(['TF_TG_SCORE'], ascending=False)
background["label"] = "background"

background_sum = sum(background["TF_TG_SCORE"])
background_length = len(background)
background_mean = background_sum/background_length
background_median=sts.median(background['TF_TG_SCORE'])
background_quantile=np.percentile(background["TF_TG_SCORE"],95)
background_quantile_99=np.percentile(background["TF_TG_SCORE"],99)
df_interesting_stats.loc[0]=['background',background_sum,background_length,background_mean,background_median,background_quantile,background_quantile_99]
row_counter=1


JUNB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/JUNB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
JUNB.columns=['TF_TG_SCORE','label']
JUNB_sum = sum(JUNB['TF_TG_SCORE'])
JUNB_length = len(JUNB)
JUNB_mean=0
JUNB_quantile=0
JUNB_quantile_95=0
JUNB_median=0
JUNB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],JUNB['TF_TG_SCORE'])
if(JUNB_length>0):
    JUNB_mean=JUNB_sum/JUNB_length
    JUNB_quantile=np.percentile(JUNB['TF_TG_SCORE'], 99)
    JUNB_quantile_95=np.percentile(JUNB['TF_TG_SCORE'], 95)
    JUNB_median=sts.median(JUNB['TF_TG_SCORE'])
if(JUNB_median > background_median and JUNB_mannwhitneyU['pvalue']<0.01):
    background_JUNB = pd.concat([background,JUNB],axis=0)
    ax_JUNB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_JUNB,palette="Set3")
    ax_JUNB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/JUNB.png')
    del background_JUNB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['JUNB',JUNB_sum,JUNB_length,JUNB_mean,JUNB_median,JUNB_quantile_95,JUNB_quantile]
    row_counter=row_counter+1
del JUNB
plt.figure(figsize=(20, 17))


NR1D2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR1D2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR1D2.columns=['TF_TG_SCORE','label']
NR1D2_sum = sum(NR1D2['TF_TG_SCORE'])
NR1D2_length = len(NR1D2)
NR1D2_mean=0
NR1D2_quantile=0
NR1D2_quantile_95=0
NR1D2_median=0
NR1D2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR1D2['TF_TG_SCORE'])
if(NR1D2_length>0):
    NR1D2_mean=NR1D2_sum/NR1D2_length
    NR1D2_quantile=np.percentile(NR1D2['TF_TG_SCORE'], 99)
    NR1D2_quantile_95=np.percentile(NR1D2['TF_TG_SCORE'], 95)
    NR1D2_median=sts.median(NR1D2['TF_TG_SCORE'])
if(NR1D2_median > background_median and NR1D2_mannwhitneyU['pvalue']<0.01):
    background_NR1D2 = pd.concat([background,NR1D2],axis=0)
    ax_NR1D2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR1D2,palette="Set3")
    ax_NR1D2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR1D2.png')
    del background_NR1D2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR1D2',NR1D2_sum,NR1D2_length,NR1D2_mean,NR1D2_median,NR1D2_quantile_95,NR1D2_quantile]
    row_counter=row_counter+1
del NR1D2
plt.figure(figsize=(20, 17))


USF2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/USF2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
USF2.columns=['TF_TG_SCORE','label']
USF2_sum = sum(USF2['TF_TG_SCORE'])
USF2_length = len(USF2)
USF2_mean=0
USF2_quantile=0
USF2_quantile_95=0
USF2_median=0
USF2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],USF2['TF_TG_SCORE'])
if(USF2_length>0):
    USF2_mean=USF2_sum/USF2_length
    USF2_quantile=np.percentile(USF2['TF_TG_SCORE'], 99)
    USF2_quantile_95=np.percentile(USF2['TF_TG_SCORE'], 95)
    USF2_median=sts.median(USF2['TF_TG_SCORE'])
if(USF2_median > background_median and USF2_mannwhitneyU['pvalue']<0.01):
    background_USF2 = pd.concat([background,USF2],axis=0)
    ax_USF2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_USF2,palette="Set3")
    ax_USF2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/USF2.png')
    del background_USF2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['USF2',USF2_sum,USF2_length,USF2_mean,USF2_median,USF2_quantile_95,USF2_quantile]
    row_counter=row_counter+1
del USF2
plt.figure(figsize=(20, 17))


MECP2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MECP2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MECP2.columns=['TF_TG_SCORE','label']
MECP2_sum = sum(MECP2['TF_TG_SCORE'])
MECP2_length = len(MECP2)
MECP2_mean=0
MECP2_quantile=0
MECP2_quantile_95=0
MECP2_median=0
MECP2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MECP2['TF_TG_SCORE'])
if(MECP2_length>0):
    MECP2_mean=MECP2_sum/MECP2_length
    MECP2_quantile=np.percentile(MECP2['TF_TG_SCORE'], 99)
    MECP2_quantile_95=np.percentile(MECP2['TF_TG_SCORE'], 95)
    MECP2_median=sts.median(MECP2['TF_TG_SCORE'])
if(MECP2_median > background_median and MECP2_mannwhitneyU['pvalue']<0.01):
    background_MECP2 = pd.concat([background,MECP2],axis=0)
    ax_MECP2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MECP2,palette="Set3")
    ax_MECP2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MECP2.png')
    del background_MECP2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MECP2',MECP2_sum,MECP2_length,MECP2_mean,MECP2_median,MECP2_quantile_95,MECP2_quantile]
    row_counter=row_counter+1
del MECP2
plt.figure(figsize=(20, 17))


SNAI2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SNAI2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SNAI2.columns=['TF_TG_SCORE','label']
SNAI2_sum = sum(SNAI2['TF_TG_SCORE'])
SNAI2_length = len(SNAI2)
SNAI2_mean=0
SNAI2_quantile=0
SNAI2_quantile_95=0
SNAI2_median=0
SNAI2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SNAI2['TF_TG_SCORE'])
if(SNAI2_length>0):
    SNAI2_mean=SNAI2_sum/SNAI2_length
    SNAI2_quantile=np.percentile(SNAI2['TF_TG_SCORE'], 99)
    SNAI2_quantile_95=np.percentile(SNAI2['TF_TG_SCORE'], 95)
    SNAI2_median=sts.median(SNAI2['TF_TG_SCORE'])
if(SNAI2_median > background_median and SNAI2_mannwhitneyU['pvalue']<0.01):
    background_SNAI2 = pd.concat([background,SNAI2],axis=0)
    ax_SNAI2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SNAI2,palette="Set3")
    ax_SNAI2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SNAI2.png')
    del background_SNAI2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SNAI2',SNAI2_sum,SNAI2_length,SNAI2_mean,SNAI2_median,SNAI2_quantile_95,SNAI2_quantile]
    row_counter=row_counter+1
del SNAI2
plt.figure(figsize=(20, 17))


HOXD8=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HOXD8_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HOXD8.columns=['TF_TG_SCORE','label']
HOXD8_sum = sum(HOXD8['TF_TG_SCORE'])
HOXD8_length = len(HOXD8)
HOXD8_mean=0
HOXD8_quantile=0
HOXD8_quantile_95=0
HOXD8_median=0
HOXD8_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HOXD8['TF_TG_SCORE'])
if(HOXD8_length>0):
    HOXD8_mean=HOXD8_sum/HOXD8_length
    HOXD8_quantile=np.percentile(HOXD8['TF_TG_SCORE'], 99)
    HOXD8_quantile_95=np.percentile(HOXD8['TF_TG_SCORE'], 95)
    HOXD8_median=sts.median(HOXD8['TF_TG_SCORE'])
if(HOXD8_median > background_median and HOXD8_mannwhitneyU['pvalue']<0.01):
    background_HOXD8 = pd.concat([background,HOXD8],axis=0)
    ax_HOXD8 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HOXD8,palette="Set3")
    ax_HOXD8.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HOXD8.png')
    del background_HOXD8
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HOXD8',HOXD8_sum,HOXD8_length,HOXD8_mean,HOXD8_median,HOXD8_quantile_95,HOXD8_quantile]
    row_counter=row_counter+1
del HOXD8
plt.figure(figsize=(20, 17))


HIC1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HIC1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HIC1.columns=['TF_TG_SCORE','label']
HIC1_sum = sum(HIC1['TF_TG_SCORE'])
HIC1_length = len(HIC1)
HIC1_mean=0
HIC1_quantile=0
HIC1_quantile_95=0
HIC1_median=0
HIC1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HIC1['TF_TG_SCORE'])
if(HIC1_length>0):
    HIC1_mean=HIC1_sum/HIC1_length
    HIC1_quantile=np.percentile(HIC1['TF_TG_SCORE'], 99)
    HIC1_quantile_95=np.percentile(HIC1['TF_TG_SCORE'], 95)
    HIC1_median=sts.median(HIC1['TF_TG_SCORE'])
if(HIC1_median > background_median and HIC1_mannwhitneyU['pvalue']<0.01):
    background_HIC1 = pd.concat([background,HIC1],axis=0)
    ax_HIC1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HIC1,palette="Set3")
    ax_HIC1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HIC1.png')
    del background_HIC1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HIC1',HIC1_sum,HIC1_length,HIC1_mean,HIC1_median,HIC1_quantile_95,HIC1_quantile]
    row_counter=row_counter+1
del HIC1
plt.figure(figsize=(20, 17))


HAND1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HAND1..TCF3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HAND1.columns=['TF_TG_SCORE','label']
HAND1_sum = sum(HAND1['TF_TG_SCORE'])
HAND1_length = len(HAND1)
HAND1_mean=0
HAND1_quantile=0
HAND1_quantile_95=0
HAND1_median=0
HAND1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HAND1['TF_TG_SCORE'])
if(HAND1_length>0):
    HAND1_mean=HAND1_sum/HAND1_length
    HAND1_quantile=np.percentile(HAND1['TF_TG_SCORE'], 99)
    HAND1_quantile_95=np.percentile(HAND1['TF_TG_SCORE'], 95)
    HAND1_median=sts.median(HAND1['TF_TG_SCORE'])
if(HAND1_median > background_median and HAND1_mannwhitneyU['pvalue']<0.01):
    background_HAND1 = pd.concat([background,HAND1],axis=0)
    ax_HAND1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HAND1,palette="Set3")
    ax_HAND1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HAND1..TCF3.png')
    del background_HAND1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HAND1..TCF3',HAND1_sum,HAND1_length,HAND1_mean,HAND1_median,HAND1_quantile_95,HAND1_quantile]
    row_counter=row_counter+1
del HAND1
plt.figure(figsize=(20, 17))


SOX9=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SOX9_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SOX9.columns=['TF_TG_SCORE','label']
SOX9_sum = sum(SOX9['TF_TG_SCORE'])
SOX9_length = len(SOX9)
SOX9_mean=0
SOX9_quantile=0
SOX9_quantile_95=0
SOX9_median=0
SOX9_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SOX9['TF_TG_SCORE'])
if(SOX9_length>0):
    SOX9_mean=SOX9_sum/SOX9_length
    SOX9_quantile=np.percentile(SOX9['TF_TG_SCORE'], 99)
    SOX9_quantile_95=np.percentile(SOX9['TF_TG_SCORE'], 95)
    SOX9_median=sts.median(SOX9['TF_TG_SCORE'])
if(SOX9_median > background_median and SOX9_mannwhitneyU['pvalue']<0.01):
    background_SOX9 = pd.concat([background,SOX9],axis=0)
    ax_SOX9 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SOX9,palette="Set3")
    ax_SOX9.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SOX9.png')
    del background_SOX9
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SOX9',SOX9_sum,SOX9_length,SOX9_mean,SOX9_median,SOX9_quantile_95,SOX9_quantile]
    row_counter=row_counter+1
del SOX9
plt.figure(figsize=(20, 17))


PBX3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PBX3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PBX3.columns=['TF_TG_SCORE','label']
PBX3_sum = sum(PBX3['TF_TG_SCORE'])
PBX3_length = len(PBX3)
PBX3_mean=0
PBX3_quantile=0
PBX3_quantile_95=0
PBX3_median=0
PBX3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PBX3['TF_TG_SCORE'])
if(PBX3_length>0):
    PBX3_mean=PBX3_sum/PBX3_length
    PBX3_quantile=np.percentile(PBX3['TF_TG_SCORE'], 99)
    PBX3_quantile_95=np.percentile(PBX3['TF_TG_SCORE'], 95)
    PBX3_median=sts.median(PBX3['TF_TG_SCORE'])
if(PBX3_median > background_median and PBX3_mannwhitneyU['pvalue']<0.01):
    background_PBX3 = pd.concat([background,PBX3],axis=0)
    ax_PBX3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PBX3,palette="Set3")
    ax_PBX3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PBX3.png')
    del background_PBX3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PBX3',PBX3_sum,PBX3_length,PBX3_mean,PBX3_median,PBX3_quantile_95,PBX3_quantile]
    row_counter=row_counter+1
del PBX3
plt.figure(figsize=(20, 17))


IRF7=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF7_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF7.columns=['TF_TG_SCORE','label']
IRF7_sum = sum(IRF7['TF_TG_SCORE'])
IRF7_length = len(IRF7)
IRF7_mean=0
IRF7_quantile=0
IRF7_quantile_95=0
IRF7_median=0
IRF7_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF7['TF_TG_SCORE'])
if(IRF7_length>0):
    IRF7_mean=IRF7_sum/IRF7_length
    IRF7_quantile=np.percentile(IRF7['TF_TG_SCORE'], 99)
    IRF7_quantile_95=np.percentile(IRF7['TF_TG_SCORE'], 95)
    IRF7_median=sts.median(IRF7['TF_TG_SCORE'])
if(IRF7_median > background_median and IRF7_mannwhitneyU['pvalue']<0.01):
    background_IRF7 = pd.concat([background,IRF7],axis=0)
    ax_IRF7 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF7,palette="Set3")
    ax_IRF7.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF7.png')
    del background_IRF7
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF7',IRF7_sum,IRF7_length,IRF7_mean,IRF7_median,IRF7_quantile_95,IRF7_quantile]
    row_counter=row_counter+1
del IRF7
plt.figure(figsize=(20, 17))


MECOM=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MECOM_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MECOM.columns=['TF_TG_SCORE','label']
MECOM_sum = sum(MECOM['TF_TG_SCORE'])
MECOM_length = len(MECOM)
MECOM_mean=0
MECOM_quantile=0
MECOM_quantile_95=0
MECOM_median=0
MECOM_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MECOM['TF_TG_SCORE'])
if(MECOM_length>0):
    MECOM_mean=MECOM_sum/MECOM_length
    MECOM_quantile=np.percentile(MECOM['TF_TG_SCORE'], 99)
    MECOM_quantile_95=np.percentile(MECOM['TF_TG_SCORE'], 95)
    MECOM_median=sts.median(MECOM['TF_TG_SCORE'])
if(MECOM_median > background_median and MECOM_mannwhitneyU['pvalue']<0.01):
    background_MECOM = pd.concat([background,MECOM],axis=0)
    ax_MECOM = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MECOM,palette="Set3")
    ax_MECOM.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MECOM.png')
    del background_MECOM
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MECOM',MECOM_sum,MECOM_length,MECOM_mean,MECOM_median,MECOM_quantile_95,MECOM_quantile]
    row_counter=row_counter+1
del MECOM
plt.figure(figsize=(20, 17))


NR3C1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR3C1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR3C1.columns=['TF_TG_SCORE','label']
NR3C1_sum = sum(NR3C1['TF_TG_SCORE'])
NR3C1_length = len(NR3C1)
NR3C1_mean=0
NR3C1_quantile=0
NR3C1_quantile_95=0
NR3C1_median=0
NR3C1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR3C1['TF_TG_SCORE'])
if(NR3C1_length>0):
    NR3C1_mean=NR3C1_sum/NR3C1_length
    NR3C1_quantile=np.percentile(NR3C1['TF_TG_SCORE'], 99)
    NR3C1_quantile_95=np.percentile(NR3C1['TF_TG_SCORE'], 95)
    NR3C1_median=sts.median(NR3C1['TF_TG_SCORE'])
if(NR3C1_median > background_median and NR3C1_mannwhitneyU['pvalue']<0.01):
    background_NR3C1 = pd.concat([background,NR3C1],axis=0)
    ax_NR3C1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR3C1,palette="Set3")
    ax_NR3C1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR3C1.png')
    del background_NR3C1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR3C1',NR3C1_sum,NR3C1_length,NR3C1_mean,NR3C1_median,NR3C1_quantile_95,NR3C1_quantile]
    row_counter=row_counter+1
del NR3C1
plt.figure(figsize=(20, 17))


SMAD3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SMAD3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SMAD3.columns=['TF_TG_SCORE','label']
SMAD3_sum = sum(SMAD3['TF_TG_SCORE'])
SMAD3_length = len(SMAD3)
SMAD3_mean=0
SMAD3_quantile=0
SMAD3_quantile_95=0
SMAD3_median=0
SMAD3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SMAD3['TF_TG_SCORE'])
if(SMAD3_length>0):
    SMAD3_mean=SMAD3_sum/SMAD3_length
    SMAD3_quantile=np.percentile(SMAD3['TF_TG_SCORE'], 99)
    SMAD3_quantile_95=np.percentile(SMAD3['TF_TG_SCORE'], 95)
    SMAD3_median=sts.median(SMAD3['TF_TG_SCORE'])
if(SMAD3_median > background_median and SMAD3_mannwhitneyU['pvalue']<0.01):
    background_SMAD3 = pd.concat([background,SMAD3],axis=0)
    ax_SMAD3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SMAD3,palette="Set3")
    ax_SMAD3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SMAD3.png')
    del background_SMAD3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SMAD3',SMAD3_sum,SMAD3_length,SMAD3_mean,SMAD3_median,SMAD3_quantile_95,SMAD3_quantile]
    row_counter=row_counter+1
del SMAD3
plt.figure(figsize=(20, 17))


ETS1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ETS1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ETS1.columns=['TF_TG_SCORE','label']
ETS1_sum = sum(ETS1['TF_TG_SCORE'])
ETS1_length = len(ETS1)
ETS1_mean=0
ETS1_quantile=0
ETS1_quantile_95=0
ETS1_median=0
ETS1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ETS1['TF_TG_SCORE'])
if(ETS1_length>0):
    ETS1_mean=ETS1_sum/ETS1_length
    ETS1_quantile=np.percentile(ETS1['TF_TG_SCORE'], 99)
    ETS1_quantile_95=np.percentile(ETS1['TF_TG_SCORE'], 95)
    ETS1_median=sts.median(ETS1['TF_TG_SCORE'])
if(ETS1_median > background_median and ETS1_mannwhitneyU['pvalue']<0.01):
    background_ETS1 = pd.concat([background,ETS1],axis=0)
    ax_ETS1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ETS1,palette="Set3")
    ax_ETS1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ETS1.png')
    del background_ETS1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ETS1',ETS1_sum,ETS1_length,ETS1_mean,ETS1_median,ETS1_quantile_95,ETS1_quantile]
    row_counter=row_counter+1
del ETS1
plt.figure(figsize=(20, 17))


CEBPD=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CEBPD_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CEBPD.columns=['TF_TG_SCORE','label']
CEBPD_sum = sum(CEBPD['TF_TG_SCORE'])
CEBPD_length = len(CEBPD)
CEBPD_mean=0
CEBPD_quantile=0
CEBPD_quantile_95=0
CEBPD_median=0
CEBPD_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CEBPD['TF_TG_SCORE'])
if(CEBPD_length>0):
    CEBPD_mean=CEBPD_sum/CEBPD_length
    CEBPD_quantile=np.percentile(CEBPD['TF_TG_SCORE'], 99)
    CEBPD_quantile_95=np.percentile(CEBPD['TF_TG_SCORE'], 95)
    CEBPD_median=sts.median(CEBPD['TF_TG_SCORE'])
if(CEBPD_median > background_median and CEBPD_mannwhitneyU['pvalue']<0.01):
    background_CEBPD = pd.concat([background,CEBPD],axis=0)
    ax_CEBPD = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CEBPD,palette="Set3")
    ax_CEBPD.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CEBPD.png')
    del background_CEBPD
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CEBPD',CEBPD_sum,CEBPD_length,CEBPD_mean,CEBPD_median,CEBPD_quantile_95,CEBPD_quantile]
    row_counter=row_counter+1
del CEBPD
plt.figure(figsize=(20, 17))


HOXA9=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HOXA9_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HOXA9.columns=['TF_TG_SCORE','label']
HOXA9_sum = sum(HOXA9['TF_TG_SCORE'])
HOXA9_length = len(HOXA9)
HOXA9_mean=0
HOXA9_quantile=0
HOXA9_quantile_95=0
HOXA9_median=0
HOXA9_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HOXA9['TF_TG_SCORE'])
if(HOXA9_length>0):
    HOXA9_mean=HOXA9_sum/HOXA9_length
    HOXA9_quantile=np.percentile(HOXA9['TF_TG_SCORE'], 99)
    HOXA9_quantile_95=np.percentile(HOXA9['TF_TG_SCORE'], 95)
    HOXA9_median=sts.median(HOXA9['TF_TG_SCORE'])
if(HOXA9_median > background_median and HOXA9_mannwhitneyU['pvalue']<0.01):
    background_HOXA9 = pd.concat([background,HOXA9],axis=0)
    ax_HOXA9 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HOXA9,palette="Set3")
    ax_HOXA9.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HOXA9.png')
    del background_HOXA9
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HOXA9',HOXA9_sum,HOXA9_length,HOXA9_mean,HOXA9_median,HOXA9_quantile_95,HOXA9_quantile]
    row_counter=row_counter+1
del HOXA9
plt.figure(figsize=(20, 17))


SOX4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SOX4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SOX4.columns=['TF_TG_SCORE','label']
SOX4_sum = sum(SOX4['TF_TG_SCORE'])
SOX4_length = len(SOX4)
SOX4_mean=0
SOX4_quantile=0
SOX4_quantile_95=0
SOX4_median=0
SOX4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SOX4['TF_TG_SCORE'])
if(SOX4_length>0):
    SOX4_mean=SOX4_sum/SOX4_length
    SOX4_quantile=np.percentile(SOX4['TF_TG_SCORE'], 99)
    SOX4_quantile_95=np.percentile(SOX4['TF_TG_SCORE'], 95)
    SOX4_median=sts.median(SOX4['TF_TG_SCORE'])
if(SOX4_median > background_median and SOX4_mannwhitneyU['pvalue']<0.01):
    background_SOX4 = pd.concat([background,SOX4],axis=0)
    ax_SOX4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SOX4,palette="Set3")
    ax_SOX4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SOX4.png')
    del background_SOX4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SOX4',SOX4_sum,SOX4_length,SOX4_mean,SOX4_median,SOX4_quantile_95,SOX4_quantile]
    row_counter=row_counter+1
del SOX4
plt.figure(figsize=(20, 17))


EGR2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/EGR2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
EGR2.columns=['TF_TG_SCORE','label']
EGR2_sum = sum(EGR2['TF_TG_SCORE'])
EGR2_length = len(EGR2)
EGR2_mean=0
EGR2_quantile=0
EGR2_quantile_95=0
EGR2_median=0
EGR2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],EGR2['TF_TG_SCORE'])
if(EGR2_length>0):
    EGR2_mean=EGR2_sum/EGR2_length
    EGR2_quantile=np.percentile(EGR2['TF_TG_SCORE'], 99)
    EGR2_quantile_95=np.percentile(EGR2['TF_TG_SCORE'], 95)
    EGR2_median=sts.median(EGR2['TF_TG_SCORE'])
if(EGR2_median > background_median and EGR2_mannwhitneyU['pvalue']<0.01):
    background_EGR2 = pd.concat([background,EGR2],axis=0)
    ax_EGR2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_EGR2,palette="Set3")
    ax_EGR2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/EGR2.png')
    del background_EGR2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['EGR2',EGR2_sum,EGR2_length,EGR2_mean,EGR2_median,EGR2_quantile_95,EGR2_quantile]
    row_counter=row_counter+1
del EGR2
plt.figure(figsize=(20, 17))


NFKB1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFKB1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFKB1.columns=['TF_TG_SCORE','label']
NFKB1_sum = sum(NFKB1['TF_TG_SCORE'])
NFKB1_length = len(NFKB1)
NFKB1_mean=0
NFKB1_quantile=0
NFKB1_quantile_95=0
NFKB1_median=0
NFKB1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFKB1['TF_TG_SCORE'])
if(NFKB1_length>0):
    NFKB1_mean=NFKB1_sum/NFKB1_length
    NFKB1_quantile=np.percentile(NFKB1['TF_TG_SCORE'], 99)
    NFKB1_quantile_95=np.percentile(NFKB1['TF_TG_SCORE'], 95)
    NFKB1_median=sts.median(NFKB1['TF_TG_SCORE'])
if(NFKB1_median > background_median and NFKB1_mannwhitneyU['pvalue']<0.01):
    background_NFKB1 = pd.concat([background,NFKB1],axis=0)
    ax_NFKB1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFKB1,palette="Set3")
    ax_NFKB1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFKB1.png')
    del background_NFKB1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFKB1',NFKB1_sum,NFKB1_length,NFKB1_mean,NFKB1_median,NFKB1_quantile_95,NFKB1_quantile]
    row_counter=row_counter+1
del NFKB1
plt.figure(figsize=(20, 17))


IRF9=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF9_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF9.columns=['TF_TG_SCORE','label']
IRF9_sum = sum(IRF9['TF_TG_SCORE'])
IRF9_length = len(IRF9)
IRF9_mean=0
IRF9_quantile=0
IRF9_quantile_95=0
IRF9_median=0
IRF9_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF9['TF_TG_SCORE'])
if(IRF9_length>0):
    IRF9_mean=IRF9_sum/IRF9_length
    IRF9_quantile=np.percentile(IRF9['TF_TG_SCORE'], 99)
    IRF9_quantile_95=np.percentile(IRF9['TF_TG_SCORE'], 95)
    IRF9_median=sts.median(IRF9['TF_TG_SCORE'])
if(IRF9_median > background_median and IRF9_mannwhitneyU['pvalue']<0.01):
    background_IRF9 = pd.concat([background,IRF9],axis=0)
    ax_IRF9 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF9,palette="Set3")
    ax_IRF9.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF9.png')
    del background_IRF9
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF9',IRF9_sum,IRF9_length,IRF9_mean,IRF9_median,IRF9_quantile_95,IRF9_quantile]
    row_counter=row_counter+1
del IRF9
plt.figure(figsize=(20, 17))


EP300=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/EP300_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
EP300.columns=['TF_TG_SCORE','label']
EP300_sum = sum(EP300['TF_TG_SCORE'])
EP300_length = len(EP300)
EP300_mean=0
EP300_quantile=0
EP300_quantile_95=0
EP300_median=0
EP300_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],EP300['TF_TG_SCORE'])
if(EP300_length>0):
    EP300_mean=EP300_sum/EP300_length
    EP300_quantile=np.percentile(EP300['TF_TG_SCORE'], 99)
    EP300_quantile_95=np.percentile(EP300['TF_TG_SCORE'], 95)
    EP300_median=sts.median(EP300['TF_TG_SCORE'])
if(EP300_median > background_median and EP300_mannwhitneyU['pvalue']<0.01):
    background_EP300 = pd.concat([background,EP300],axis=0)
    ax_EP300 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_EP300,palette="Set3")
    ax_EP300.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/EP300.png')
    del background_EP300
    plt.clf()
    df_interesting_stats.loc[row_counter]=['EP300',EP300_sum,EP300_length,EP300_mean,EP300_median,EP300_quantile_95,EP300_quantile]
    row_counter=row_counter+1
del EP300
plt.figure(figsize=(20, 17))


RUNX2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RUNX2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RUNX2.columns=['TF_TG_SCORE','label']
RUNX2_sum = sum(RUNX2['TF_TG_SCORE'])
RUNX2_length = len(RUNX2)
RUNX2_mean=0
RUNX2_quantile=0
RUNX2_quantile_95=0
RUNX2_median=0
RUNX2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RUNX2['TF_TG_SCORE'])
if(RUNX2_length>0):
    RUNX2_mean=RUNX2_sum/RUNX2_length
    RUNX2_quantile=np.percentile(RUNX2['TF_TG_SCORE'], 99)
    RUNX2_quantile_95=np.percentile(RUNX2['TF_TG_SCORE'], 95)
    RUNX2_median=sts.median(RUNX2['TF_TG_SCORE'])
if(RUNX2_median > background_median and RUNX2_mannwhitneyU['pvalue']<0.01):
    background_RUNX2 = pd.concat([background,RUNX2],axis=0)
    ax_RUNX2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RUNX2,palette="Set3")
    ax_RUNX2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RUNX2.png')
    del background_RUNX2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RUNX2',RUNX2_sum,RUNX2_length,RUNX2_mean,RUNX2_median,RUNX2_quantile_95,RUNX2_quantile]
    row_counter=row_counter+1
del RUNX2
plt.figure(figsize=(20, 17))


VDR=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/VDR_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
VDR.columns=['TF_TG_SCORE','label']
VDR_sum = sum(VDR['TF_TG_SCORE'])
VDR_length = len(VDR)
VDR_mean=0
VDR_quantile=0
VDR_quantile_95=0
VDR_median=0
VDR_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],VDR['TF_TG_SCORE'])
if(VDR_length>0):
    VDR_mean=VDR_sum/VDR_length
    VDR_quantile=np.percentile(VDR['TF_TG_SCORE'], 99)
    VDR_quantile_95=np.percentile(VDR['TF_TG_SCORE'], 95)
    VDR_median=sts.median(VDR['TF_TG_SCORE'])
if(VDR_median > background_median and VDR_mannwhitneyU['pvalue']<0.01):
    background_VDR = pd.concat([background,VDR],axis=0)
    ax_VDR = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_VDR,palette="Set3")
    ax_VDR.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/VDR.png')
    del background_VDR
    plt.clf()
    df_interesting_stats.loc[row_counter]=['VDR',VDR_sum,VDR_length,VDR_mean,VDR_median,VDR_quantile_95,VDR_quantile]
    row_counter=row_counter+1
del VDR
plt.figure(figsize=(20, 17))


KLF4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/KLF4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
KLF4.columns=['TF_TG_SCORE','label']
KLF4_sum = sum(KLF4['TF_TG_SCORE'])
KLF4_length = len(KLF4)
KLF4_mean=0
KLF4_quantile=0
KLF4_quantile_95=0
KLF4_median=0
KLF4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],KLF4['TF_TG_SCORE'])
if(KLF4_length>0):
    KLF4_mean=KLF4_sum/KLF4_length
    KLF4_quantile=np.percentile(KLF4['TF_TG_SCORE'], 99)
    KLF4_quantile_95=np.percentile(KLF4['TF_TG_SCORE'], 95)
    KLF4_median=sts.median(KLF4['TF_TG_SCORE'])
if(KLF4_median > background_median and KLF4_mannwhitneyU['pvalue']<0.01):
    background_KLF4 = pd.concat([background,KLF4],axis=0)
    ax_KLF4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_KLF4,palette="Set3")
    ax_KLF4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/KLF4.png')
    del background_KLF4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['KLF4',KLF4_sum,KLF4_length,KLF4_mean,KLF4_median,KLF4_quantile_95,KLF4_quantile]
    row_counter=row_counter+1
del KLF4
plt.figure(figsize=(20, 17))


BDP1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BDP1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BDP1.columns=['TF_TG_SCORE','label']
BDP1_sum = sum(BDP1['TF_TG_SCORE'])
BDP1_length = len(BDP1)
BDP1_mean=0
BDP1_quantile=0
BDP1_quantile_95=0
BDP1_median=0
BDP1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BDP1['TF_TG_SCORE'])
if(BDP1_length>0):
    BDP1_mean=BDP1_sum/BDP1_length
    BDP1_quantile=np.percentile(BDP1['TF_TG_SCORE'], 99)
    BDP1_quantile_95=np.percentile(BDP1['TF_TG_SCORE'], 95)
    BDP1_median=sts.median(BDP1['TF_TG_SCORE'])
if(BDP1_median > background_median and BDP1_mannwhitneyU['pvalue']<0.01):
    background_BDP1 = pd.concat([background,BDP1],axis=0)
    ax_BDP1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BDP1,palette="Set3")
    ax_BDP1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BDP1.png')
    del background_BDP1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BDP1',BDP1_sum,BDP1_length,BDP1_mean,BDP1_median,BDP1_quantile_95,BDP1_quantile]
    row_counter=row_counter+1
del BDP1
plt.figure(figsize=(20, 17))


CLOCK=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CLOCK_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CLOCK.columns=['TF_TG_SCORE','label']
CLOCK_sum = sum(CLOCK['TF_TG_SCORE'])
CLOCK_length = len(CLOCK)
CLOCK_mean=0
CLOCK_quantile=0
CLOCK_quantile_95=0
CLOCK_median=0
CLOCK_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CLOCK['TF_TG_SCORE'])
if(CLOCK_length>0):
    CLOCK_mean=CLOCK_sum/CLOCK_length
    CLOCK_quantile=np.percentile(CLOCK['TF_TG_SCORE'], 99)
    CLOCK_quantile_95=np.percentile(CLOCK['TF_TG_SCORE'], 95)
    CLOCK_median=sts.median(CLOCK['TF_TG_SCORE'])
if(CLOCK_median > background_median and CLOCK_mannwhitneyU['pvalue']<0.01):
    background_CLOCK = pd.concat([background,CLOCK],axis=0)
    ax_CLOCK = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CLOCK,palette="Set3")
    ax_CLOCK.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CLOCK.png')
    del background_CLOCK
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CLOCK',CLOCK_sum,CLOCK_length,CLOCK_mean,CLOCK_median,CLOCK_quantile_95,CLOCK_quantile]
    row_counter=row_counter+1
del CLOCK
plt.figure(figsize=(20, 17))


SP3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SP3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SP3.columns=['TF_TG_SCORE','label']
SP3_sum = sum(SP3['TF_TG_SCORE'])
SP3_length = len(SP3)
SP3_mean=0
SP3_quantile=0
SP3_quantile_95=0
SP3_median=0
SP3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SP3['TF_TG_SCORE'])
if(SP3_length>0):
    SP3_mean=SP3_sum/SP3_length
    SP3_quantile=np.percentile(SP3['TF_TG_SCORE'], 99)
    SP3_quantile_95=np.percentile(SP3['TF_TG_SCORE'], 95)
    SP3_median=sts.median(SP3['TF_TG_SCORE'])
if(SP3_median > background_median and SP3_mannwhitneyU['pvalue']<0.01):
    background_SP3 = pd.concat([background,SP3],axis=0)
    ax_SP3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SP3,palette="Set3")
    ax_SP3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SP3.png')
    del background_SP3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SP3',SP3_sum,SP3_length,SP3_mean,SP3_median,SP3_quantile_95,SP3_quantile]
    row_counter=row_counter+1
del SP3
plt.figure(figsize=(20, 17))


MAZ=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MAZ_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MAZ.columns=['TF_TG_SCORE','label']
MAZ_sum = sum(MAZ['TF_TG_SCORE'])
MAZ_length = len(MAZ)
MAZ_mean=0
MAZ_quantile=0
MAZ_quantile_95=0
MAZ_median=0
MAZ_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MAZ['TF_TG_SCORE'])
if(MAZ_length>0):
    MAZ_mean=MAZ_sum/MAZ_length
    MAZ_quantile=np.percentile(MAZ['TF_TG_SCORE'], 99)
    MAZ_quantile_95=np.percentile(MAZ['TF_TG_SCORE'], 95)
    MAZ_median=sts.median(MAZ['TF_TG_SCORE'])
if(MAZ_median > background_median and MAZ_mannwhitneyU['pvalue']<0.01):
    background_MAZ = pd.concat([background,MAZ],axis=0)
    ax_MAZ = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MAZ,palette="Set3")
    ax_MAZ.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MAZ.png')
    del background_MAZ
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MAZ',MAZ_sum,MAZ_length,MAZ_mean,MAZ_median,MAZ_quantile_95,MAZ_quantile]
    row_counter=row_counter+1
del MAZ
plt.figure(figsize=(20, 17))


NFIC=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFIC_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFIC.columns=['TF_TG_SCORE','label']
NFIC_sum = sum(NFIC['TF_TG_SCORE'])
NFIC_length = len(NFIC)
NFIC_mean=0
NFIC_quantile=0
NFIC_quantile_95=0
NFIC_median=0
NFIC_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFIC['TF_TG_SCORE'])
if(NFIC_length>0):
    NFIC_mean=NFIC_sum/NFIC_length
    NFIC_quantile=np.percentile(NFIC['TF_TG_SCORE'], 99)
    NFIC_quantile_95=np.percentile(NFIC['TF_TG_SCORE'], 95)
    NFIC_median=sts.median(NFIC['TF_TG_SCORE'])
if(NFIC_median > background_median and NFIC_mannwhitneyU['pvalue']<0.01):
    background_NFIC = pd.concat([background,NFIC],axis=0)
    ax_NFIC = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFIC,palette="Set3")
    ax_NFIC.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFIC.png')
    del background_NFIC
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFIC',NFIC_sum,NFIC_length,NFIC_mean,NFIC_median,NFIC_quantile_95,NFIC_quantile]
    row_counter=row_counter+1
del NFIC
plt.figure(figsize=(20, 17))


OVOL1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/OVOL1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
OVOL1.columns=['TF_TG_SCORE','label']
OVOL1_sum = sum(OVOL1['TF_TG_SCORE'])
OVOL1_length = len(OVOL1)
OVOL1_mean=0
OVOL1_quantile=0
OVOL1_quantile_95=0
OVOL1_median=0
OVOL1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],OVOL1['TF_TG_SCORE'])
if(OVOL1_length>0):
    OVOL1_mean=OVOL1_sum/OVOL1_length
    OVOL1_quantile=np.percentile(OVOL1['TF_TG_SCORE'], 99)
    OVOL1_quantile_95=np.percentile(OVOL1['TF_TG_SCORE'], 95)
    OVOL1_median=sts.median(OVOL1['TF_TG_SCORE'])
if(OVOL1_median > background_median and OVOL1_mannwhitneyU['pvalue']<0.01):
    background_OVOL1 = pd.concat([background,OVOL1],axis=0)
    ax_OVOL1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_OVOL1,palette="Set3")
    ax_OVOL1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/OVOL1.png')
    del background_OVOL1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['OVOL1',OVOL1_sum,OVOL1_length,OVOL1_mean,OVOL1_median,OVOL1_quantile_95,OVOL1_quantile]
    row_counter=row_counter+1
del OVOL1
plt.figure(figsize=(20, 17))


NR2C2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR2C2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR2C2.columns=['TF_TG_SCORE','label']
NR2C2_sum = sum(NR2C2['TF_TG_SCORE'])
NR2C2_length = len(NR2C2)
NR2C2_mean=0
NR2C2_quantile=0
NR2C2_quantile_95=0
NR2C2_median=0
NR2C2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR2C2['TF_TG_SCORE'])
if(NR2C2_length>0):
    NR2C2_mean=NR2C2_sum/NR2C2_length
    NR2C2_quantile=np.percentile(NR2C2['TF_TG_SCORE'], 99)
    NR2C2_quantile_95=np.percentile(NR2C2['TF_TG_SCORE'], 95)
    NR2C2_median=sts.median(NR2C2['TF_TG_SCORE'])
if(NR2C2_median > background_median and NR2C2_mannwhitneyU['pvalue']<0.01):
    background_NR2C2 = pd.concat([background,NR2C2],axis=0)
    ax_NR2C2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR2C2,palette="Set3")
    ax_NR2C2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR2C2.png')
    del background_NR2C2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR2C2',NR2C2_sum,NR2C2_length,NR2C2_mean,NR2C2_median,NR2C2_quantile_95,NR2C2_quantile]
    row_counter=row_counter+1
del NR2C2
plt.figure(figsize=(20, 17))


MLXIP=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MLXIP_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MLXIP.columns=['TF_TG_SCORE','label']
MLXIP_sum = sum(MLXIP['TF_TG_SCORE'])
MLXIP_length = len(MLXIP)
MLXIP_mean=0
MLXIP_quantile=0
MLXIP_quantile_95=0
MLXIP_median=0
MLXIP_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MLXIP['TF_TG_SCORE'])
if(MLXIP_length>0):
    MLXIP_mean=MLXIP_sum/MLXIP_length
    MLXIP_quantile=np.percentile(MLXIP['TF_TG_SCORE'], 99)
    MLXIP_quantile_95=np.percentile(MLXIP['TF_TG_SCORE'], 95)
    MLXIP_median=sts.median(MLXIP['TF_TG_SCORE'])
if(MLXIP_median > background_median and MLXIP_mannwhitneyU['pvalue']<0.01):
    background_MLXIP = pd.concat([background,MLXIP],axis=0)
    ax_MLXIP = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MLXIP,palette="Set3")
    ax_MLXIP.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MLXIP.png')
    del background_MLXIP
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MLXIP',MLXIP_sum,MLXIP_length,MLXIP_mean,MLXIP_median,MLXIP_quantile_95,MLXIP_quantile]
    row_counter=row_counter+1
del MLXIP
plt.figure(figsize=(20, 17))


NR1H3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR1H3..RXRA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR1H3.columns=['TF_TG_SCORE','label']
NR1H3_sum = sum(NR1H3['TF_TG_SCORE'])
NR1H3_length = len(NR1H3)
NR1H3_mean=0
NR1H3_quantile=0
NR1H3_quantile_95=0
NR1H3_median=0
NR1H3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR1H3['TF_TG_SCORE'])
if(NR1H3_length>0):
    NR1H3_mean=NR1H3_sum/NR1H3_length
    NR1H3_quantile=np.percentile(NR1H3['TF_TG_SCORE'], 99)
    NR1H3_quantile_95=np.percentile(NR1H3['TF_TG_SCORE'], 95)
    NR1H3_median=sts.median(NR1H3['TF_TG_SCORE'])
if(NR1H3_median > background_median and NR1H3_mannwhitneyU['pvalue']<0.01):
    background_NR1H3 = pd.concat([background,NR1H3],axis=0)
    ax_NR1H3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR1H3,palette="Set3")
    ax_NR1H3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR1H3..RXRA.png')
    del background_NR1H3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR1H3..RXRA',NR1H3_sum,NR1H3_length,NR1H3_mean,NR1H3_median,NR1H3_quantile_95,NR1H3_quantile]
    row_counter=row_counter+1
del NR1H3
plt.figure(figsize=(20, 17))


RUNX3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RUNX3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RUNX3.columns=['TF_TG_SCORE','label']
RUNX3_sum = sum(RUNX3['TF_TG_SCORE'])
RUNX3_length = len(RUNX3)
RUNX3_mean=0
RUNX3_quantile=0
RUNX3_quantile_95=0
RUNX3_median=0
RUNX3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RUNX3['TF_TG_SCORE'])
if(RUNX3_length>0):
    RUNX3_mean=RUNX3_sum/RUNX3_length
    RUNX3_quantile=np.percentile(RUNX3['TF_TG_SCORE'], 99)
    RUNX3_quantile_95=np.percentile(RUNX3['TF_TG_SCORE'], 95)
    RUNX3_median=sts.median(RUNX3['TF_TG_SCORE'])
if(RUNX3_median > background_median and RUNX3_mannwhitneyU['pvalue']<0.01):
    background_RUNX3 = pd.concat([background,RUNX3],axis=0)
    ax_RUNX3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RUNX3,palette="Set3")
    ax_RUNX3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RUNX3.png')
    del background_RUNX3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RUNX3',RUNX3_sum,RUNX3_length,RUNX3_mean,RUNX3_median,RUNX3_quantile_95,RUNX3_quantile]
    row_counter=row_counter+1
del RUNX3
plt.figure(figsize=(20, 17))


TBP=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TBP_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TBP.columns=['TF_TG_SCORE','label']
TBP_sum = sum(TBP['TF_TG_SCORE'])
TBP_length = len(TBP)
TBP_mean=0
TBP_quantile=0
TBP_quantile_95=0
TBP_median=0
TBP_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TBP['TF_TG_SCORE'])
if(TBP_length>0):
    TBP_mean=TBP_sum/TBP_length
    TBP_quantile=np.percentile(TBP['TF_TG_SCORE'], 99)
    TBP_quantile_95=np.percentile(TBP['TF_TG_SCORE'], 95)
    TBP_median=sts.median(TBP['TF_TG_SCORE'])
if(TBP_median > background_median and TBP_mannwhitneyU['pvalue']<0.01):
    background_TBP = pd.concat([background,TBP],axis=0)
    ax_TBP = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TBP,palette="Set3")
    ax_TBP.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TBP.png')
    del background_TBP
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TBP',TBP_sum,TBP_length,TBP_mean,TBP_median,TBP_quantile_95,TBP_quantile]
    row_counter=row_counter+1
del TBP
plt.figure(figsize=(20, 17))


MYC=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MYC_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MYC.columns=['TF_TG_SCORE','label']
MYC_sum = sum(MYC['TF_TG_SCORE'])
MYC_length = len(MYC)
MYC_mean=0
MYC_quantile=0
MYC_quantile_95=0
MYC_median=0
MYC_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MYC['TF_TG_SCORE'])
if(MYC_length>0):
    MYC_mean=MYC_sum/MYC_length
    MYC_quantile=np.percentile(MYC['TF_TG_SCORE'], 99)
    MYC_quantile_95=np.percentile(MYC['TF_TG_SCORE'], 95)
    MYC_median=sts.median(MYC['TF_TG_SCORE'])
if(MYC_median > background_median and MYC_mannwhitneyU['pvalue']<0.01):
    background_MYC = pd.concat([background,MYC],axis=0)
    ax_MYC = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MYC,palette="Set3")
    ax_MYC.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MYC.png')
    del background_MYC
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MYC',MYC_sum,MYC_length,MYC_mean,MYC_median,MYC_quantile_95,MYC_quantile]
    row_counter=row_counter+1
del MYC
plt.figure(figsize=(20, 17))


BHLHA15=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BHLHA15_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BHLHA15.columns=['TF_TG_SCORE','label']
BHLHA15_sum = sum(BHLHA15['TF_TG_SCORE'])
BHLHA15_length = len(BHLHA15)
BHLHA15_mean=0
BHLHA15_quantile=0
BHLHA15_quantile_95=0
BHLHA15_median=0
BHLHA15_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BHLHA15['TF_TG_SCORE'])
if(BHLHA15_length>0):
    BHLHA15_mean=BHLHA15_sum/BHLHA15_length
    BHLHA15_quantile=np.percentile(BHLHA15['TF_TG_SCORE'], 99)
    BHLHA15_quantile_95=np.percentile(BHLHA15['TF_TG_SCORE'], 95)
    BHLHA15_median=sts.median(BHLHA15['TF_TG_SCORE'])
if(BHLHA15_median > background_median and BHLHA15_mannwhitneyU['pvalue']<0.01):
    background_BHLHA15 = pd.concat([background,BHLHA15],axis=0)
    ax_BHLHA15 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BHLHA15,palette="Set3")
    ax_BHLHA15.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BHLHA15.png')
    del background_BHLHA15
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BHLHA15',BHLHA15_sum,BHLHA15_length,BHLHA15_mean,BHLHA15_median,BHLHA15_quantile_95,BHLHA15_quantile]
    row_counter=row_counter+1
del BHLHA15
plt.figure(figsize=(20, 17))


E2F1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/E2F1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
E2F1.columns=['TF_TG_SCORE','label']
E2F1_sum = sum(E2F1['TF_TG_SCORE'])
E2F1_length = len(E2F1)
E2F1_mean=0
E2F1_quantile=0
E2F1_quantile_95=0
E2F1_median=0
E2F1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],E2F1['TF_TG_SCORE'])
if(E2F1_length>0):
    E2F1_mean=E2F1_sum/E2F1_length
    E2F1_quantile=np.percentile(E2F1['TF_TG_SCORE'], 99)
    E2F1_quantile_95=np.percentile(E2F1['TF_TG_SCORE'], 95)
    E2F1_median=sts.median(E2F1['TF_TG_SCORE'])
if(E2F1_median > background_median and E2F1_mannwhitneyU['pvalue']<0.01):
    background_E2F1 = pd.concat([background,E2F1],axis=0)
    ax_E2F1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_E2F1,palette="Set3")
    ax_E2F1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/E2F1.png')
    del background_E2F1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['E2F1',E2F1_sum,E2F1_length,E2F1_mean,E2F1_median,E2F1_quantile_95,E2F1_quantile]
    row_counter=row_counter+1
del E2F1
plt.figure(figsize=(20, 17))


HSF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HSF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HSF1.columns=['TF_TG_SCORE','label']
HSF1_sum = sum(HSF1['TF_TG_SCORE'])
HSF1_length = len(HSF1)
HSF1_mean=0
HSF1_quantile=0
HSF1_quantile_95=0
HSF1_median=0
HSF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HSF1['TF_TG_SCORE'])
if(HSF1_length>0):
    HSF1_mean=HSF1_sum/HSF1_length
    HSF1_quantile=np.percentile(HSF1['TF_TG_SCORE'], 99)
    HSF1_quantile_95=np.percentile(HSF1['TF_TG_SCORE'], 95)
    HSF1_median=sts.median(HSF1['TF_TG_SCORE'])
if(HSF1_median > background_median and HSF1_mannwhitneyU['pvalue']<0.01):
    background_HSF1 = pd.concat([background,HSF1],axis=0)
    ax_HSF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HSF1,palette="Set3")
    ax_HSF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HSF1.png')
    del background_HSF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HSF1',HSF1_sum,HSF1_length,HSF1_mean,HSF1_median,HSF1_quantile_95,HSF1_quantile]
    row_counter=row_counter+1
del HSF1
plt.figure(figsize=(20, 17))


FOXI1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXI1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXI1.columns=['TF_TG_SCORE','label']
FOXI1_sum = sum(FOXI1['TF_TG_SCORE'])
FOXI1_length = len(FOXI1)
FOXI1_mean=0
FOXI1_quantile=0
FOXI1_quantile_95=0
FOXI1_median=0
FOXI1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXI1['TF_TG_SCORE'])
if(FOXI1_length>0):
    FOXI1_mean=FOXI1_sum/FOXI1_length
    FOXI1_quantile=np.percentile(FOXI1['TF_TG_SCORE'], 99)
    FOXI1_quantile_95=np.percentile(FOXI1['TF_TG_SCORE'], 95)
    FOXI1_median=sts.median(FOXI1['TF_TG_SCORE'])
if(FOXI1_median > background_median and FOXI1_mannwhitneyU['pvalue']<0.01):
    background_FOXI1 = pd.concat([background,FOXI1],axis=0)
    ax_FOXI1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXI1,palette="Set3")
    ax_FOXI1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXI1.png')
    del background_FOXI1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXI1',FOXI1_sum,FOXI1_length,FOXI1_mean,FOXI1_median,FOXI1_quantile_95,FOXI1_quantile]
    row_counter=row_counter+1
del FOXI1
plt.figure(figsize=(20, 17))


SOX6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SOX6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SOX6.columns=['TF_TG_SCORE','label']
SOX6_sum = sum(SOX6['TF_TG_SCORE'])
SOX6_length = len(SOX6)
SOX6_mean=0
SOX6_quantile=0
SOX6_quantile_95=0
SOX6_median=0
SOX6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SOX6['TF_TG_SCORE'])
if(SOX6_length>0):
    SOX6_mean=SOX6_sum/SOX6_length
    SOX6_quantile=np.percentile(SOX6['TF_TG_SCORE'], 99)
    SOX6_quantile_95=np.percentile(SOX6['TF_TG_SCORE'], 95)
    SOX6_median=sts.median(SOX6['TF_TG_SCORE'])
if(SOX6_median > background_median and SOX6_mannwhitneyU['pvalue']<0.01):
    background_SOX6 = pd.concat([background,SOX6],axis=0)
    ax_SOX6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SOX6,palette="Set3")
    ax_SOX6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SOX6.png')
    del background_SOX6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SOX6',SOX6_sum,SOX6_length,SOX6_mean,SOX6_median,SOX6_quantile_95,SOX6_quantile]
    row_counter=row_counter+1
del SOX6
plt.figure(figsize=(20, 17))


HES1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HES1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HES1.columns=['TF_TG_SCORE','label']
HES1_sum = sum(HES1['TF_TG_SCORE'])
HES1_length = len(HES1)
HES1_mean=0
HES1_quantile=0
HES1_quantile_95=0
HES1_median=0
HES1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HES1['TF_TG_SCORE'])
if(HES1_length>0):
    HES1_mean=HES1_sum/HES1_length
    HES1_quantile=np.percentile(HES1['TF_TG_SCORE'], 99)
    HES1_quantile_95=np.percentile(HES1['TF_TG_SCORE'], 95)
    HES1_median=sts.median(HES1['TF_TG_SCORE'])
if(HES1_median > background_median and HES1_mannwhitneyU['pvalue']<0.01):
    background_HES1 = pd.concat([background,HES1],axis=0)
    ax_HES1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HES1,palette="Set3")
    ax_HES1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HES1.png')
    del background_HES1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HES1',HES1_sum,HES1_length,HES1_mean,HES1_median,HES1_quantile_95,HES1_quantile]
    row_counter=row_counter+1
del HES1
plt.figure(figsize=(20, 17))


RELB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RELB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RELB.columns=['TF_TG_SCORE','label']
RELB_sum = sum(RELB['TF_TG_SCORE'])
RELB_length = len(RELB)
RELB_mean=0
RELB_quantile=0
RELB_quantile_95=0
RELB_median=0
RELB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RELB['TF_TG_SCORE'])
if(RELB_length>0):
    RELB_mean=RELB_sum/RELB_length
    RELB_quantile=np.percentile(RELB['TF_TG_SCORE'], 99)
    RELB_quantile_95=np.percentile(RELB['TF_TG_SCORE'], 95)
    RELB_median=sts.median(RELB['TF_TG_SCORE'])
if(RELB_median > background_median and RELB_mannwhitneyU['pvalue']<0.01):
    background_RELB = pd.concat([background,RELB],axis=0)
    ax_RELB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RELB,palette="Set3")
    ax_RELB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RELB.png')
    del background_RELB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RELB',RELB_sum,RELB_length,RELB_mean,RELB_median,RELB_quantile_95,RELB_quantile]
    row_counter=row_counter+1
del RELB
plt.figure(figsize=(20, 17))


RARG=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RARG_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RARG.columns=['TF_TG_SCORE','label']
RARG_sum = sum(RARG['TF_TG_SCORE'])
RARG_length = len(RARG)
RARG_mean=0
RARG_quantile=0
RARG_quantile_95=0
RARG_median=0
RARG_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RARG['TF_TG_SCORE'])
if(RARG_length>0):
    RARG_mean=RARG_sum/RARG_length
    RARG_quantile=np.percentile(RARG['TF_TG_SCORE'], 99)
    RARG_quantile_95=np.percentile(RARG['TF_TG_SCORE'], 95)
    RARG_median=sts.median(RARG['TF_TG_SCORE'])
if(RARG_median > background_median and RARG_mannwhitneyU['pvalue']<0.01):
    background_RARG = pd.concat([background,RARG],axis=0)
    ax_RARG = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RARG,palette="Set3")
    ax_RARG.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RARG.png')
    del background_RARG
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RARG',RARG_sum,RARG_length,RARG_mean,RARG_median,RARG_quantile_95,RARG_quantile]
    row_counter=row_counter+1
del RARG
plt.figure(figsize=(20, 17))


CTCF=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CTCF_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CTCF.columns=['TF_TG_SCORE','label']
CTCF_sum = sum(CTCF['TF_TG_SCORE'])
CTCF_length = len(CTCF)
CTCF_mean=0
CTCF_quantile=0
CTCF_quantile_95=0
CTCF_median=0
CTCF_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CTCF['TF_TG_SCORE'])
if(CTCF_length>0):
    CTCF_mean=CTCF_sum/CTCF_length
    CTCF_quantile=np.percentile(CTCF['TF_TG_SCORE'], 99)
    CTCF_quantile_95=np.percentile(CTCF['TF_TG_SCORE'], 95)
    CTCF_median=sts.median(CTCF['TF_TG_SCORE'])
if(CTCF_median > background_median and CTCF_mannwhitneyU['pvalue']<0.01):
    background_CTCF = pd.concat([background,CTCF],axis=0)
    ax_CTCF = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CTCF,palette="Set3")
    ax_CTCF.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CTCF.png')
    del background_CTCF
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CTCF',CTCF_sum,CTCF_length,CTCF_mean,CTCF_median,CTCF_quantile_95,CTCF_quantile]
    row_counter=row_counter+1
del CTCF
plt.figure(figsize=(20, 17))


LHX6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/LHX6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
LHX6.columns=['TF_TG_SCORE','label']
LHX6_sum = sum(LHX6['TF_TG_SCORE'])
LHX6_length = len(LHX6)
LHX6_mean=0
LHX6_quantile=0
LHX6_quantile_95=0
LHX6_median=0
LHX6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],LHX6['TF_TG_SCORE'])
if(LHX6_length>0):
    LHX6_mean=LHX6_sum/LHX6_length
    LHX6_quantile=np.percentile(LHX6['TF_TG_SCORE'], 99)
    LHX6_quantile_95=np.percentile(LHX6['TF_TG_SCORE'], 95)
    LHX6_median=sts.median(LHX6['TF_TG_SCORE'])
if(LHX6_median > background_median and LHX6_mannwhitneyU['pvalue']<0.01):
    background_LHX6 = pd.concat([background,LHX6],axis=0)
    ax_LHX6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_LHX6,palette="Set3")
    ax_LHX6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/LHX6.png')
    del background_LHX6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['LHX6',LHX6_sum,LHX6_length,LHX6_mean,LHX6_median,LHX6_quantile_95,LHX6_quantile]
    row_counter=row_counter+1
del LHX6
plt.figure(figsize=(20, 17))


MAF=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MAF_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MAF.columns=['TF_TG_SCORE','label']
MAF_sum = sum(MAF['TF_TG_SCORE'])
MAF_length = len(MAF)
MAF_mean=0
MAF_quantile=0
MAF_quantile_95=0
MAF_median=0
MAF_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MAF['TF_TG_SCORE'])
if(MAF_length>0):
    MAF_mean=MAF_sum/MAF_length
    MAF_quantile=np.percentile(MAF['TF_TG_SCORE'], 99)
    MAF_quantile_95=np.percentile(MAF['TF_TG_SCORE'], 95)
    MAF_median=sts.median(MAF['TF_TG_SCORE'])
if(MAF_median > background_median and MAF_mannwhitneyU['pvalue']<0.01):
    background_MAF = pd.concat([background,MAF],axis=0)
    ax_MAF = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MAF,palette="Set3")
    ax_MAF.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MAF.png')
    del background_MAF
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MAF',MAF_sum,MAF_length,MAF_mean,MAF_median,MAF_quantile_95,MAF_quantile]
    row_counter=row_counter+1
del MAF
plt.figure(figsize=(20, 17))


STAT6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/STAT6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
STAT6.columns=['TF_TG_SCORE','label']
STAT6_sum = sum(STAT6['TF_TG_SCORE'])
STAT6_length = len(STAT6)
STAT6_mean=0
STAT6_quantile=0
STAT6_quantile_95=0
STAT6_median=0
STAT6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],STAT6['TF_TG_SCORE'])
if(STAT6_length>0):
    STAT6_mean=STAT6_sum/STAT6_length
    STAT6_quantile=np.percentile(STAT6['TF_TG_SCORE'], 99)
    STAT6_quantile_95=np.percentile(STAT6['TF_TG_SCORE'], 95)
    STAT6_median=sts.median(STAT6['TF_TG_SCORE'])
if(STAT6_median > background_median and STAT6_mannwhitneyU['pvalue']<0.01):
    background_STAT6 = pd.concat([background,STAT6],axis=0)
    ax_STAT6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_STAT6,palette="Set3")
    ax_STAT6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/STAT6.png')
    del background_STAT6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['STAT6',STAT6_sum,STAT6_length,STAT6_mean,STAT6_median,STAT6_quantile_95,STAT6_quantile]
    row_counter=row_counter+1
del STAT6
plt.figure(figsize=(20, 17))


MAFG=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MAFG_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MAFG.columns=['TF_TG_SCORE','label']
MAFG_sum = sum(MAFG['TF_TG_SCORE'])
MAFG_length = len(MAFG)
MAFG_mean=0
MAFG_quantile=0
MAFG_quantile_95=0
MAFG_median=0
MAFG_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MAFG['TF_TG_SCORE'])
if(MAFG_length>0):
    MAFG_mean=MAFG_sum/MAFG_length
    MAFG_quantile=np.percentile(MAFG['TF_TG_SCORE'], 99)
    MAFG_quantile_95=np.percentile(MAFG['TF_TG_SCORE'], 95)
    MAFG_median=sts.median(MAFG['TF_TG_SCORE'])
if(MAFG_median > background_median and MAFG_mannwhitneyU['pvalue']<0.01):
    background_MAFG = pd.concat([background,MAFG],axis=0)
    ax_MAFG = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MAFG,palette="Set3")
    ax_MAFG.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MAFG.png')
    del background_MAFG
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MAFG',MAFG_sum,MAFG_length,MAFG_mean,MAFG_median,MAFG_quantile_95,MAFG_quantile]
    row_counter=row_counter+1
del MAFG
plt.figure(figsize=(20, 17))


MYB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MYB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MYB.columns=['TF_TG_SCORE','label']
MYB_sum = sum(MYB['TF_TG_SCORE'])
MYB_length = len(MYB)
MYB_mean=0
MYB_quantile=0
MYB_quantile_95=0
MYB_median=0
MYB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MYB['TF_TG_SCORE'])
if(MYB_length>0):
    MYB_mean=MYB_sum/MYB_length
    MYB_quantile=np.percentile(MYB['TF_TG_SCORE'], 99)
    MYB_quantile_95=np.percentile(MYB['TF_TG_SCORE'], 95)
    MYB_median=sts.median(MYB['TF_TG_SCORE'])
if(MYB_median > background_median and MYB_mannwhitneyU['pvalue']<0.01):
    background_MYB = pd.concat([background,MYB],axis=0)
    ax_MYB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MYB,palette="Set3")
    ax_MYB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MYB.png')
    del background_MYB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MYB',MYB_sum,MYB_length,MYB_mean,MYB_median,MYB_quantile_95,MYB_quantile]
    row_counter=row_counter+1
del MYB
plt.figure(figsize=(20, 17))


NFYA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFYA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFYA.columns=['TF_TG_SCORE','label']
NFYA_sum = sum(NFYA['TF_TG_SCORE'])
NFYA_length = len(NFYA)
NFYA_mean=0
NFYA_quantile=0
NFYA_quantile_95=0
NFYA_median=0
NFYA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFYA['TF_TG_SCORE'])
if(NFYA_length>0):
    NFYA_mean=NFYA_sum/NFYA_length
    NFYA_quantile=np.percentile(NFYA['TF_TG_SCORE'], 99)
    NFYA_quantile_95=np.percentile(NFYA['TF_TG_SCORE'], 95)
    NFYA_median=sts.median(NFYA['TF_TG_SCORE'])
if(NFYA_median > background_median and NFYA_mannwhitneyU['pvalue']<0.01):
    background_NFYA = pd.concat([background,NFYA],axis=0)
    ax_NFYA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFYA,palette="Set3")
    ax_NFYA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFYA.png')
    del background_NFYA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFYA',NFYA_sum,NFYA_length,NFYA_mean,NFYA_median,NFYA_quantile_95,NFYA_quantile]
    row_counter=row_counter+1
del NFYA
plt.figure(figsize=(20, 17))


NR4A1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR4A1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR4A1.columns=['TF_TG_SCORE','label']
NR4A1_sum = sum(NR4A1['TF_TG_SCORE'])
NR4A1_length = len(NR4A1)
NR4A1_mean=0
NR4A1_quantile=0
NR4A1_quantile_95=0
NR4A1_median=0
NR4A1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR4A1['TF_TG_SCORE'])
if(NR4A1_length>0):
    NR4A1_mean=NR4A1_sum/NR4A1_length
    NR4A1_quantile=np.percentile(NR4A1['TF_TG_SCORE'], 99)
    NR4A1_quantile_95=np.percentile(NR4A1['TF_TG_SCORE'], 95)
    NR4A1_median=sts.median(NR4A1['TF_TG_SCORE'])
if(NR4A1_median > background_median and NR4A1_mannwhitneyU['pvalue']<0.01):
    background_NR4A1 = pd.concat([background,NR4A1],axis=0)
    ax_NR4A1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR4A1,palette="Set3")
    ax_NR4A1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR4A1.png')
    del background_NR4A1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR4A1',NR4A1_sum,NR4A1_length,NR4A1_mean,NR4A1_median,NR4A1_quantile_95,NR4A1_quantile]
    row_counter=row_counter+1
del NR4A1
plt.figure(figsize=(20, 17))


SMAD2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SMAD2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SMAD2.columns=['TF_TG_SCORE','label']
SMAD2_sum = sum(SMAD2['TF_TG_SCORE'])
SMAD2_length = len(SMAD2)
SMAD2_mean=0
SMAD2_quantile=0
SMAD2_quantile_95=0
SMAD2_median=0
SMAD2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SMAD2['TF_TG_SCORE'])
if(SMAD2_length>0):
    SMAD2_mean=SMAD2_sum/SMAD2_length
    SMAD2_quantile=np.percentile(SMAD2['TF_TG_SCORE'], 99)
    SMAD2_quantile_95=np.percentile(SMAD2['TF_TG_SCORE'], 95)
    SMAD2_median=sts.median(SMAD2['TF_TG_SCORE'])
if(SMAD2_median > background_median and SMAD2_mannwhitneyU['pvalue']<0.01):
    background_SMAD2 = pd.concat([background,SMAD2],axis=0)
    ax_SMAD2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SMAD2,palette="Set3")
    ax_SMAD2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SMAD2.png')
    del background_SMAD2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SMAD2',SMAD2_sum,SMAD2_length,SMAD2_mean,SMAD2_median,SMAD2_quantile_95,SMAD2_quantile]
    row_counter=row_counter+1
del SMAD2
plt.figure(figsize=(20, 17))


CEBPB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CEBPB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CEBPB.columns=['TF_TG_SCORE','label']
CEBPB_sum = sum(CEBPB['TF_TG_SCORE'])
CEBPB_length = len(CEBPB)
CEBPB_mean=0
CEBPB_quantile=0
CEBPB_quantile_95=0
CEBPB_median=0
CEBPB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CEBPB['TF_TG_SCORE'])
if(CEBPB_length>0):
    CEBPB_mean=CEBPB_sum/CEBPB_length
    CEBPB_quantile=np.percentile(CEBPB['TF_TG_SCORE'], 99)
    CEBPB_quantile_95=np.percentile(CEBPB['TF_TG_SCORE'], 95)
    CEBPB_median=sts.median(CEBPB['TF_TG_SCORE'])
if(CEBPB_median > background_median and CEBPB_mannwhitneyU['pvalue']<0.01):
    background_CEBPB = pd.concat([background,CEBPB],axis=0)
    ax_CEBPB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CEBPB,palette="Set3")
    ax_CEBPB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CEBPB.png')
    del background_CEBPB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CEBPB',CEBPB_sum,CEBPB_length,CEBPB_mean,CEBPB_median,CEBPB_quantile_95,CEBPB_quantile]
    row_counter=row_counter+1
del CEBPB
plt.figure(figsize=(20, 17))


CREB1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CREB1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CREB1.columns=['TF_TG_SCORE','label']
CREB1_sum = sum(CREB1['TF_TG_SCORE'])
CREB1_length = len(CREB1)
CREB1_mean=0
CREB1_quantile=0
CREB1_quantile_95=0
CREB1_median=0
CREB1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CREB1['TF_TG_SCORE'])
if(CREB1_length>0):
    CREB1_mean=CREB1_sum/CREB1_length
    CREB1_quantile=np.percentile(CREB1['TF_TG_SCORE'], 99)
    CREB1_quantile_95=np.percentile(CREB1['TF_TG_SCORE'], 95)
    CREB1_median=sts.median(CREB1['TF_TG_SCORE'])
if(CREB1_median > background_median and CREB1_mannwhitneyU['pvalue']<0.01):
    background_CREB1 = pd.concat([background,CREB1],axis=0)
    ax_CREB1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CREB1,palette="Set3")
    ax_CREB1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CREB1.png')
    del background_CREB1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CREB1',CREB1_sum,CREB1_length,CREB1_mean,CREB1_median,CREB1_quantile_95,CREB1_quantile]
    row_counter=row_counter+1
del CREB1
plt.figure(figsize=(20, 17))


MTF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MTF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MTF1.columns=['TF_TG_SCORE','label']
MTF1_sum = sum(MTF1['TF_TG_SCORE'])
MTF1_length = len(MTF1)
MTF1_mean=0
MTF1_quantile=0
MTF1_quantile_95=0
MTF1_median=0
MTF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MTF1['TF_TG_SCORE'])
if(MTF1_length>0):
    MTF1_mean=MTF1_sum/MTF1_length
    MTF1_quantile=np.percentile(MTF1['TF_TG_SCORE'], 99)
    MTF1_quantile_95=np.percentile(MTF1['TF_TG_SCORE'], 95)
    MTF1_median=sts.median(MTF1['TF_TG_SCORE'])
if(MTF1_median > background_median and MTF1_mannwhitneyU['pvalue']<0.01):
    background_MTF1 = pd.concat([background,MTF1],axis=0)
    ax_MTF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MTF1,palette="Set3")
    ax_MTF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MTF1.png')
    del background_MTF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MTF1',MTF1_sum,MTF1_length,MTF1_mean,MTF1_median,MTF1_quantile_95,MTF1_quantile]
    row_counter=row_counter+1
del MTF1
plt.figure(figsize=(20, 17))


MEF2A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MEF2A_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MEF2A.columns=['TF_TG_SCORE','label']
MEF2A_sum = sum(MEF2A['TF_TG_SCORE'])
MEF2A_length = len(MEF2A)
MEF2A_mean=0
MEF2A_quantile=0
MEF2A_quantile_95=0
MEF2A_median=0
MEF2A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MEF2A['TF_TG_SCORE'])
if(MEF2A_length>0):
    MEF2A_mean=MEF2A_sum/MEF2A_length
    MEF2A_quantile=np.percentile(MEF2A['TF_TG_SCORE'], 99)
    MEF2A_quantile_95=np.percentile(MEF2A['TF_TG_SCORE'], 95)
    MEF2A_median=sts.median(MEF2A['TF_TG_SCORE'])
if(MEF2A_median > background_median and MEF2A_mannwhitneyU['pvalue']<0.01):
    background_MEF2A = pd.concat([background,MEF2A],axis=0)
    ax_MEF2A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MEF2A,palette="Set3")
    ax_MEF2A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MEF2A.png')
    del background_MEF2A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MEF2A',MEF2A_sum,MEF2A_length,MEF2A_mean,MEF2A_median,MEF2A_quantile_95,MEF2A_quantile]
    row_counter=row_counter+1
del MEF2A
plt.figure(figsize=(20, 17))


FOXJ2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXJ2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXJ2.columns=['TF_TG_SCORE','label']
FOXJ2_sum = sum(FOXJ2['TF_TG_SCORE'])
FOXJ2_length = len(FOXJ2)
FOXJ2_mean=0
FOXJ2_quantile=0
FOXJ2_quantile_95=0
FOXJ2_median=0
FOXJ2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXJ2['TF_TG_SCORE'])
if(FOXJ2_length>0):
    FOXJ2_mean=FOXJ2_sum/FOXJ2_length
    FOXJ2_quantile=np.percentile(FOXJ2['TF_TG_SCORE'], 99)
    FOXJ2_quantile_95=np.percentile(FOXJ2['TF_TG_SCORE'], 95)
    FOXJ2_median=sts.median(FOXJ2['TF_TG_SCORE'])
if(FOXJ2_median > background_median and FOXJ2_mannwhitneyU['pvalue']<0.01):
    background_FOXJ2 = pd.concat([background,FOXJ2],axis=0)
    ax_FOXJ2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXJ2,palette="Set3")
    ax_FOXJ2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXJ2.png')
    del background_FOXJ2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXJ2',FOXJ2_sum,FOXJ2_length,FOXJ2_mean,FOXJ2_median,FOXJ2_quantile_95,FOXJ2_quantile]
    row_counter=row_counter+1
del FOXJ2
plt.figure(figsize=(20, 17))


ZBTB33=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ZBTB33_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ZBTB33.columns=['TF_TG_SCORE','label']
ZBTB33_sum = sum(ZBTB33['TF_TG_SCORE'])
ZBTB33_length = len(ZBTB33)
ZBTB33_mean=0
ZBTB33_quantile=0
ZBTB33_quantile_95=0
ZBTB33_median=0
ZBTB33_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ZBTB33['TF_TG_SCORE'])
if(ZBTB33_length>0):
    ZBTB33_mean=ZBTB33_sum/ZBTB33_length
    ZBTB33_quantile=np.percentile(ZBTB33['TF_TG_SCORE'], 99)
    ZBTB33_quantile_95=np.percentile(ZBTB33['TF_TG_SCORE'], 95)
    ZBTB33_median=sts.median(ZBTB33['TF_TG_SCORE'])
if(ZBTB33_median > background_median and ZBTB33_mannwhitneyU['pvalue']<0.01):
    background_ZBTB33 = pd.concat([background,ZBTB33],axis=0)
    ax_ZBTB33 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ZBTB33,palette="Set3")
    ax_ZBTB33.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ZBTB33.png')
    del background_ZBTB33
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ZBTB33',ZBTB33_sum,ZBTB33_length,ZBTB33_mean,ZBTB33_median,ZBTB33_quantile_95,ZBTB33_quantile]
    row_counter=row_counter+1
del ZBTB33
plt.figure(figsize=(20, 17))


FOXM1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXM1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXM1.columns=['TF_TG_SCORE','label']
FOXM1_sum = sum(FOXM1['TF_TG_SCORE'])
FOXM1_length = len(FOXM1)
FOXM1_mean=0
FOXM1_quantile=0
FOXM1_quantile_95=0
FOXM1_median=0
FOXM1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXM1['TF_TG_SCORE'])
if(FOXM1_length>0):
    FOXM1_mean=FOXM1_sum/FOXM1_length
    FOXM1_quantile=np.percentile(FOXM1['TF_TG_SCORE'], 99)
    FOXM1_quantile_95=np.percentile(FOXM1['TF_TG_SCORE'], 95)
    FOXM1_median=sts.median(FOXM1['TF_TG_SCORE'])
if(FOXM1_median > background_median and FOXM1_mannwhitneyU['pvalue']<0.01):
    background_FOXM1 = pd.concat([background,FOXM1],axis=0)
    ax_FOXM1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXM1,palette="Set3")
    ax_FOXM1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXM1.png')
    del background_FOXM1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXM1',FOXM1_sum,FOXM1_length,FOXM1_mean,FOXM1_median,FOXM1_quantile_95,FOXM1_quantile]
    row_counter=row_counter+1
del FOXM1
plt.figure(figsize=(20, 17))


ETV4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ETV4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ETV4.columns=['TF_TG_SCORE','label']
ETV4_sum = sum(ETV4['TF_TG_SCORE'])
ETV4_length = len(ETV4)
ETV4_mean=0
ETV4_quantile=0
ETV4_quantile_95=0
ETV4_median=0
ETV4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ETV4['TF_TG_SCORE'])
if(ETV4_length>0):
    ETV4_mean=ETV4_sum/ETV4_length
    ETV4_quantile=np.percentile(ETV4['TF_TG_SCORE'], 99)
    ETV4_quantile_95=np.percentile(ETV4['TF_TG_SCORE'], 95)
    ETV4_median=sts.median(ETV4['TF_TG_SCORE'])
if(ETV4_median > background_median and ETV4_mannwhitneyU['pvalue']<0.01):
    background_ETV4 = pd.concat([background,ETV4],axis=0)
    ax_ETV4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ETV4,palette="Set3")
    ax_ETV4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ETV4.png')
    del background_ETV4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ETV4',ETV4_sum,ETV4_length,ETV4_mean,ETV4_median,ETV4_quantile_95,ETV4_quantile]
    row_counter=row_counter+1
del ETV4
plt.figure(figsize=(20, 17))


NR2F6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR2F6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR2F6.columns=['TF_TG_SCORE','label']
NR2F6_sum = sum(NR2F6['TF_TG_SCORE'])
NR2F6_length = len(NR2F6)
NR2F6_mean=0
NR2F6_quantile=0
NR2F6_quantile_95=0
NR2F6_median=0
NR2F6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR2F6['TF_TG_SCORE'])
if(NR2F6_length>0):
    NR2F6_mean=NR2F6_sum/NR2F6_length
    NR2F6_quantile=np.percentile(NR2F6['TF_TG_SCORE'], 99)
    NR2F6_quantile_95=np.percentile(NR2F6['TF_TG_SCORE'], 95)
    NR2F6_median=sts.median(NR2F6['TF_TG_SCORE'])
if(NR2F6_median > background_median and NR2F6_mannwhitneyU['pvalue']<0.01):
    background_NR2F6 = pd.concat([background,NR2F6],axis=0)
    ax_NR2F6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR2F6,palette="Set3")
    ax_NR2F6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR2F6.png')
    del background_NR2F6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR2F6',NR2F6_sum,NR2F6_length,NR2F6_mean,NR2F6_median,NR2F6_quantile_95,NR2F6_quantile]
    row_counter=row_counter+1
del NR2F6
plt.figure(figsize=(20, 17))


FOXJ3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXJ3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXJ3.columns=['TF_TG_SCORE','label']
FOXJ3_sum = sum(FOXJ3['TF_TG_SCORE'])
FOXJ3_length = len(FOXJ3)
FOXJ3_mean=0
FOXJ3_quantile=0
FOXJ3_quantile_95=0
FOXJ3_median=0
FOXJ3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXJ3['TF_TG_SCORE'])
if(FOXJ3_length>0):
    FOXJ3_mean=FOXJ3_sum/FOXJ3_length
    FOXJ3_quantile=np.percentile(FOXJ3['TF_TG_SCORE'], 99)
    FOXJ3_quantile_95=np.percentile(FOXJ3['TF_TG_SCORE'], 95)
    FOXJ3_median=sts.median(FOXJ3['TF_TG_SCORE'])
if(FOXJ3_median > background_median and FOXJ3_mannwhitneyU['pvalue']<0.01):
    background_FOXJ3 = pd.concat([background,FOXJ3],axis=0)
    ax_FOXJ3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXJ3,palette="Set3")
    ax_FOXJ3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXJ3.png')
    del background_FOXJ3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXJ3',FOXJ3_sum,FOXJ3_length,FOXJ3_mean,FOXJ3_median,FOXJ3_quantile_95,FOXJ3_quantile]
    row_counter=row_counter+1
del FOXJ3
plt.figure(figsize=(20, 17))


GATA3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/GATA3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
GATA3.columns=['TF_TG_SCORE','label']
GATA3_sum = sum(GATA3['TF_TG_SCORE'])
GATA3_length = len(GATA3)
GATA3_mean=0
GATA3_quantile=0
GATA3_quantile_95=0
GATA3_median=0
GATA3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],GATA3['TF_TG_SCORE'])
if(GATA3_length>0):
    GATA3_mean=GATA3_sum/GATA3_length
    GATA3_quantile=np.percentile(GATA3['TF_TG_SCORE'], 99)
    GATA3_quantile_95=np.percentile(GATA3['TF_TG_SCORE'], 95)
    GATA3_median=sts.median(GATA3['TF_TG_SCORE'])
if(GATA3_median > background_median and GATA3_mannwhitneyU['pvalue']<0.01):
    background_GATA3 = pd.concat([background,GATA3],axis=0)
    ax_GATA3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_GATA3,palette="Set3")
    ax_GATA3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/GATA3.png')
    del background_GATA3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['GATA3',GATA3_sum,GATA3_length,GATA3_mean,GATA3_median,GATA3_quantile_95,GATA3_quantile]
    row_counter=row_counter+1
del GATA3
plt.figure(figsize=(20, 17))


RXRB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RXRB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RXRB.columns=['TF_TG_SCORE','label']
RXRB_sum = sum(RXRB['TF_TG_SCORE'])
RXRB_length = len(RXRB)
RXRB_mean=0
RXRB_quantile=0
RXRB_quantile_95=0
RXRB_median=0
RXRB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RXRB['TF_TG_SCORE'])
if(RXRB_length>0):
    RXRB_mean=RXRB_sum/RXRB_length
    RXRB_quantile=np.percentile(RXRB['TF_TG_SCORE'], 99)
    RXRB_quantile_95=np.percentile(RXRB['TF_TG_SCORE'], 95)
    RXRB_median=sts.median(RXRB['TF_TG_SCORE'])
if(RXRB_median > background_median and RXRB_mannwhitneyU['pvalue']<0.01):
    background_RXRB = pd.concat([background,RXRB],axis=0)
    ax_RXRB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RXRB,palette="Set3")
    ax_RXRB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RXRB.png')
    del background_RXRB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RXRB',RXRB_sum,RXRB_length,RXRB_mean,RXRB_median,RXRB_quantile_95,RXRB_quantile]
    row_counter=row_counter+1
del RXRB
plt.figure(figsize=(20, 17))


ESRRB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ESRRB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ESRRB.columns=['TF_TG_SCORE','label']
ESRRB_sum = sum(ESRRB['TF_TG_SCORE'])
ESRRB_length = len(ESRRB)
ESRRB_mean=0
ESRRB_quantile=0
ESRRB_quantile_95=0
ESRRB_median=0
ESRRB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ESRRB['TF_TG_SCORE'])
if(ESRRB_length>0):
    ESRRB_mean=ESRRB_sum/ESRRB_length
    ESRRB_quantile=np.percentile(ESRRB['TF_TG_SCORE'], 99)
    ESRRB_quantile_95=np.percentile(ESRRB['TF_TG_SCORE'], 95)
    ESRRB_median=sts.median(ESRRB['TF_TG_SCORE'])
if(ESRRB_median > background_median and ESRRB_mannwhitneyU['pvalue']<0.01):
    background_ESRRB = pd.concat([background,ESRRB],axis=0)
    ax_ESRRB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ESRRB,palette="Set3")
    ax_ESRRB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ESRRB.png')
    del background_ESRRB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ESRRB',ESRRB_sum,ESRRB_length,ESRRB_mean,ESRRB_median,ESRRB_quantile_95,ESRRB_quantile]
    row_counter=row_counter+1
del ESRRB
plt.figure(figsize=(20, 17))


SIRT6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SIRT6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SIRT6.columns=['TF_TG_SCORE','label']
SIRT6_sum = sum(SIRT6['TF_TG_SCORE'])
SIRT6_length = len(SIRT6)
SIRT6_mean=0
SIRT6_quantile=0
SIRT6_quantile_95=0
SIRT6_median=0
SIRT6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SIRT6['TF_TG_SCORE'])
if(SIRT6_length>0):
    SIRT6_mean=SIRT6_sum/SIRT6_length
    SIRT6_quantile=np.percentile(SIRT6['TF_TG_SCORE'], 99)
    SIRT6_quantile_95=np.percentile(SIRT6['TF_TG_SCORE'], 95)
    SIRT6_median=sts.median(SIRT6['TF_TG_SCORE'])
if(SIRT6_median > background_median and SIRT6_mannwhitneyU['pvalue']<0.01):
    background_SIRT6 = pd.concat([background,SIRT6],axis=0)
    ax_SIRT6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SIRT6,palette="Set3")
    ax_SIRT6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SIRT6.png')
    del background_SIRT6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SIRT6',SIRT6_sum,SIRT6_length,SIRT6_mean,SIRT6_median,SIRT6_quantile_95,SIRT6_quantile]
    row_counter=row_counter+1
del SIRT6
plt.figure(figsize=(20, 17))


TAL1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TAL1..GATA1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TAL1.columns=['TF_TG_SCORE','label']
TAL1_sum = sum(TAL1['TF_TG_SCORE'])
TAL1_length = len(TAL1)
TAL1_mean=0
TAL1_quantile=0
TAL1_quantile_95=0
TAL1_median=0
TAL1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TAL1['TF_TG_SCORE'])
if(TAL1_length>0):
    TAL1_mean=TAL1_sum/TAL1_length
    TAL1_quantile=np.percentile(TAL1['TF_TG_SCORE'], 99)
    TAL1_quantile_95=np.percentile(TAL1['TF_TG_SCORE'], 95)
    TAL1_median=sts.median(TAL1['TF_TG_SCORE'])
if(TAL1_median > background_median and TAL1_mannwhitneyU['pvalue']<0.01):
    background_TAL1 = pd.concat([background,TAL1],axis=0)
    ax_TAL1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TAL1,palette="Set3")
    ax_TAL1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TAL1..GATA1.png')
    del background_TAL1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TAL1..GATA1',TAL1_sum,TAL1_length,TAL1_mean,TAL1_median,TAL1_quantile_95,TAL1_quantile]
    row_counter=row_counter+1
del TAL1
plt.figure(figsize=(20, 17))


ZBTB7A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ZBTB7A_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ZBTB7A.columns=['TF_TG_SCORE','label']
ZBTB7A_sum = sum(ZBTB7A['TF_TG_SCORE'])
ZBTB7A_length = len(ZBTB7A)
ZBTB7A_mean=0
ZBTB7A_quantile=0
ZBTB7A_quantile_95=0
ZBTB7A_median=0
ZBTB7A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ZBTB7A['TF_TG_SCORE'])
if(ZBTB7A_length>0):
    ZBTB7A_mean=ZBTB7A_sum/ZBTB7A_length
    ZBTB7A_quantile=np.percentile(ZBTB7A['TF_TG_SCORE'], 99)
    ZBTB7A_quantile_95=np.percentile(ZBTB7A['TF_TG_SCORE'], 95)
    ZBTB7A_median=sts.median(ZBTB7A['TF_TG_SCORE'])
if(ZBTB7A_median > background_median and ZBTB7A_mannwhitneyU['pvalue']<0.01):
    background_ZBTB7A = pd.concat([background,ZBTB7A],axis=0)
    ax_ZBTB7A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ZBTB7A,palette="Set3")
    ax_ZBTB7A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ZBTB7A.png')
    del background_ZBTB7A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ZBTB7A',ZBTB7A_sum,ZBTB7A_length,ZBTB7A_mean,ZBTB7A_median,ZBTB7A_quantile_95,ZBTB7A_quantile]
    row_counter=row_counter+1
del ZBTB7A
plt.figure(figsize=(20, 17))


ATF2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ATF2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ATF2.columns=['TF_TG_SCORE','label']
ATF2_sum = sum(ATF2['TF_TG_SCORE'])
ATF2_length = len(ATF2)
ATF2_mean=0
ATF2_quantile=0
ATF2_quantile_95=0
ATF2_median=0
ATF2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ATF2['TF_TG_SCORE'])
if(ATF2_length>0):
    ATF2_mean=ATF2_sum/ATF2_length
    ATF2_quantile=np.percentile(ATF2['TF_TG_SCORE'], 99)
    ATF2_quantile_95=np.percentile(ATF2['TF_TG_SCORE'], 95)
    ATF2_median=sts.median(ATF2['TF_TG_SCORE'])
if(ATF2_median > background_median and ATF2_mannwhitneyU['pvalue']<0.01):
    background_ATF2 = pd.concat([background,ATF2],axis=0)
    ax_ATF2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ATF2,palette="Set3")
    ax_ATF2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ATF2.png')
    del background_ATF2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ATF2',ATF2_sum,ATF2_length,ATF2_mean,ATF2_median,ATF2_quantile_95,ATF2_quantile]
    row_counter=row_counter+1
del ATF2
plt.figure(figsize=(20, 17))


NFYB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFYB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFYB.columns=['TF_TG_SCORE','label']
NFYB_sum = sum(NFYB['TF_TG_SCORE'])
NFYB_length = len(NFYB)
NFYB_mean=0
NFYB_quantile=0
NFYB_quantile_95=0
NFYB_median=0
NFYB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFYB['TF_TG_SCORE'])
if(NFYB_length>0):
    NFYB_mean=NFYB_sum/NFYB_length
    NFYB_quantile=np.percentile(NFYB['TF_TG_SCORE'], 99)
    NFYB_quantile_95=np.percentile(NFYB['TF_TG_SCORE'], 95)
    NFYB_median=sts.median(NFYB['TF_TG_SCORE'])
if(NFYB_median > background_median and NFYB_mannwhitneyU['pvalue']<0.01):
    background_NFYB = pd.concat([background,NFYB],axis=0)
    ax_NFYB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFYB,palette="Set3")
    ax_NFYB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFYB.png')
    del background_NFYB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFYB',NFYB_sum,NFYB_length,NFYB_mean,NFYB_median,NFYB_quantile_95,NFYB_quantile]
    row_counter=row_counter+1
del NFYB
plt.figure(figsize=(20, 17))


FOXO1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXO1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXO1.columns=['TF_TG_SCORE','label']
FOXO1_sum = sum(FOXO1['TF_TG_SCORE'])
FOXO1_length = len(FOXO1)
FOXO1_mean=0
FOXO1_quantile=0
FOXO1_quantile_95=0
FOXO1_median=0
FOXO1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXO1['TF_TG_SCORE'])
if(FOXO1_length>0):
    FOXO1_mean=FOXO1_sum/FOXO1_length
    FOXO1_quantile=np.percentile(FOXO1['TF_TG_SCORE'], 99)
    FOXO1_quantile_95=np.percentile(FOXO1['TF_TG_SCORE'], 95)
    FOXO1_median=sts.median(FOXO1['TF_TG_SCORE'])
if(FOXO1_median > background_median and FOXO1_mannwhitneyU['pvalue']<0.01):
    background_FOXO1 = pd.concat([background,FOXO1],axis=0)
    ax_FOXO1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXO1,palette="Set3")
    ax_FOXO1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXO1.png')
    del background_FOXO1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXO1',FOXO1_sum,FOXO1_length,FOXO1_mean,FOXO1_median,FOXO1_quantile_95,FOXO1_quantile]
    row_counter=row_counter+1
del FOXO1
plt.figure(figsize=(20, 17))


FOXO4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXO4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXO4.columns=['TF_TG_SCORE','label']
FOXO4_sum = sum(FOXO4['TF_TG_SCORE'])
FOXO4_length = len(FOXO4)
FOXO4_mean=0
FOXO4_quantile=0
FOXO4_quantile_95=0
FOXO4_median=0
FOXO4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXO4['TF_TG_SCORE'])
if(FOXO4_length>0):
    FOXO4_mean=FOXO4_sum/FOXO4_length
    FOXO4_quantile=np.percentile(FOXO4['TF_TG_SCORE'], 99)
    FOXO4_quantile_95=np.percentile(FOXO4['TF_TG_SCORE'], 95)
    FOXO4_median=sts.median(FOXO4['TF_TG_SCORE'])
if(FOXO4_median > background_median and FOXO4_mannwhitneyU['pvalue']<0.01):
    background_FOXO4 = pd.concat([background,FOXO4],axis=0)
    ax_FOXO4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXO4,palette="Set3")
    ax_FOXO4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXO4.png')
    del background_FOXO4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXO4',FOXO4_sum,FOXO4_length,FOXO4_mean,FOXO4_median,FOXO4_quantile_95,FOXO4_quantile]
    row_counter=row_counter+1
del FOXO4
plt.figure(figsize=(20, 17))


HMGN3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HMGN3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HMGN3.columns=['TF_TG_SCORE','label']
HMGN3_sum = sum(HMGN3['TF_TG_SCORE'])
HMGN3_length = len(HMGN3)
HMGN3_mean=0
HMGN3_quantile=0
HMGN3_quantile_95=0
HMGN3_median=0
HMGN3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HMGN3['TF_TG_SCORE'])
if(HMGN3_length>0):
    HMGN3_mean=HMGN3_sum/HMGN3_length
    HMGN3_quantile=np.percentile(HMGN3['TF_TG_SCORE'], 99)
    HMGN3_quantile_95=np.percentile(HMGN3['TF_TG_SCORE'], 95)
    HMGN3_median=sts.median(HMGN3['TF_TG_SCORE'])
if(HMGN3_median > background_median and HMGN3_mannwhitneyU['pvalue']<0.01):
    background_HMGN3 = pd.concat([background,HMGN3],axis=0)
    ax_HMGN3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HMGN3,palette="Set3")
    ax_HMGN3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HMGN3.png')
    del background_HMGN3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HMGN3',HMGN3_sum,HMGN3_length,HMGN3_mean,HMGN3_median,HMGN3_quantile_95,HMGN3_quantile]
    row_counter=row_counter+1
del HMGN3
plt.figure(figsize=(20, 17))


IRF8=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF8_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF8.columns=['TF_TG_SCORE','label']
IRF8_sum = sum(IRF8['TF_TG_SCORE'])
IRF8_length = len(IRF8)
IRF8_mean=0
IRF8_quantile=0
IRF8_quantile_95=0
IRF8_median=0
IRF8_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF8['TF_TG_SCORE'])
if(IRF8_length>0):
    IRF8_mean=IRF8_sum/IRF8_length
    IRF8_quantile=np.percentile(IRF8['TF_TG_SCORE'], 99)
    IRF8_quantile_95=np.percentile(IRF8['TF_TG_SCORE'], 95)
    IRF8_median=sts.median(IRF8['TF_TG_SCORE'])
if(IRF8_median > background_median and IRF8_mannwhitneyU['pvalue']<0.01):
    background_IRF8 = pd.concat([background,IRF8],axis=0)
    ax_IRF8 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF8,palette="Set3")
    ax_IRF8.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF8.png')
    del background_IRF8
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF8',IRF8_sum,IRF8_length,IRF8_mean,IRF8_median,IRF8_quantile_95,IRF8_quantile]
    row_counter=row_counter+1
del IRF8
plt.figure(figsize=(20, 17))


SOX10=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SOX10_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SOX10.columns=['TF_TG_SCORE','label']
SOX10_sum = sum(SOX10['TF_TG_SCORE'])
SOX10_length = len(SOX10)
SOX10_mean=0
SOX10_quantile=0
SOX10_quantile_95=0
SOX10_median=0
SOX10_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SOX10['TF_TG_SCORE'])
if(SOX10_length>0):
    SOX10_mean=SOX10_sum/SOX10_length
    SOX10_quantile=np.percentile(SOX10['TF_TG_SCORE'], 99)
    SOX10_quantile_95=np.percentile(SOX10['TF_TG_SCORE'], 95)
    SOX10_median=sts.median(SOX10['TF_TG_SCORE'])
if(SOX10_median > background_median and SOX10_mannwhitneyU['pvalue']<0.01):
    background_SOX10 = pd.concat([background,SOX10],axis=0)
    ax_SOX10 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SOX10,palette="Set3")
    ax_SOX10.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SOX10.png')
    del background_SOX10
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SOX10',SOX10_sum,SOX10_length,SOX10_mean,SOX10_median,SOX10_quantile_95,SOX10_quantile]
    row_counter=row_counter+1
del SOX10
plt.figure(figsize=(20, 17))


RAD21=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RAD21_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RAD21.columns=['TF_TG_SCORE','label']
RAD21_sum = sum(RAD21['TF_TG_SCORE'])
RAD21_length = len(RAD21)
RAD21_mean=0
RAD21_quantile=0
RAD21_quantile_95=0
RAD21_median=0
RAD21_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RAD21['TF_TG_SCORE'])
if(RAD21_length>0):
    RAD21_mean=RAD21_sum/RAD21_length
    RAD21_quantile=np.percentile(RAD21['TF_TG_SCORE'], 99)
    RAD21_quantile_95=np.percentile(RAD21['TF_TG_SCORE'], 95)
    RAD21_median=sts.median(RAD21['TF_TG_SCORE'])
if(RAD21_median > background_median and RAD21_mannwhitneyU['pvalue']<0.01):
    background_RAD21 = pd.concat([background,RAD21],axis=0)
    ax_RAD21 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RAD21,palette="Set3")
    ax_RAD21.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RAD21.png')
    del background_RAD21
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RAD21',RAD21_sum,RAD21_length,RAD21_mean,RAD21_median,RAD21_quantile_95,RAD21_quantile]
    row_counter=row_counter+1
del RAD21
plt.figure(figsize=(20, 17))


IRF3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF3.columns=['TF_TG_SCORE','label']
IRF3_sum = sum(IRF3['TF_TG_SCORE'])
IRF3_length = len(IRF3)
IRF3_mean=0
IRF3_quantile=0
IRF3_quantile_95=0
IRF3_median=0
IRF3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF3['TF_TG_SCORE'])
if(IRF3_length>0):
    IRF3_mean=IRF3_sum/IRF3_length
    IRF3_quantile=np.percentile(IRF3['TF_TG_SCORE'], 99)
    IRF3_quantile_95=np.percentile(IRF3['TF_TG_SCORE'], 95)
    IRF3_median=sts.median(IRF3['TF_TG_SCORE'])
if(IRF3_median > background_median and IRF3_mannwhitneyU['pvalue']<0.01):
    background_IRF3 = pd.concat([background,IRF3],axis=0)
    ax_IRF3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF3,palette="Set3")
    ax_IRF3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF3.png')
    del background_IRF3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF3',IRF3_sum,IRF3_length,IRF3_mean,IRF3_median,IRF3_quantile_95,IRF3_quantile]
    row_counter=row_counter+1
del IRF3
plt.figure(figsize=(20, 17))


ELF2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ELF2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ELF2.columns=['TF_TG_SCORE','label']
ELF2_sum = sum(ELF2['TF_TG_SCORE'])
ELF2_length = len(ELF2)
ELF2_mean=0
ELF2_quantile=0
ELF2_quantile_95=0
ELF2_median=0
ELF2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ELF2['TF_TG_SCORE'])
if(ELF2_length>0):
    ELF2_mean=ELF2_sum/ELF2_length
    ELF2_quantile=np.percentile(ELF2['TF_TG_SCORE'], 99)
    ELF2_quantile_95=np.percentile(ELF2['TF_TG_SCORE'], 95)
    ELF2_median=sts.median(ELF2['TF_TG_SCORE'])
if(ELF2_median > background_median and ELF2_mannwhitneyU['pvalue']<0.01):
    background_ELF2 = pd.concat([background,ELF2],axis=0)
    ax_ELF2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ELF2,palette="Set3")
    ax_ELF2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ELF2.png')
    del background_ELF2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ELF2',ELF2_sum,ELF2_length,ELF2_mean,ELF2_median,ELF2_quantile_95,ELF2_quantile]
    row_counter=row_counter+1
del ELF2
plt.figure(figsize=(20, 17))


HINFP=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HINFP_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HINFP.columns=['TF_TG_SCORE','label']
HINFP_sum = sum(HINFP['TF_TG_SCORE'])
HINFP_length = len(HINFP)
HINFP_mean=0
HINFP_quantile=0
HINFP_quantile_95=0
HINFP_median=0
HINFP_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HINFP['TF_TG_SCORE'])
if(HINFP_length>0):
    HINFP_mean=HINFP_sum/HINFP_length
    HINFP_quantile=np.percentile(HINFP['TF_TG_SCORE'], 99)
    HINFP_quantile_95=np.percentile(HINFP['TF_TG_SCORE'], 95)
    HINFP_median=sts.median(HINFP['TF_TG_SCORE'])
if(HINFP_median > background_median and HINFP_mannwhitneyU['pvalue']<0.01):
    background_HINFP = pd.concat([background,HINFP],axis=0)
    ax_HINFP = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HINFP,palette="Set3")
    ax_HINFP.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HINFP.png')
    del background_HINFP
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HINFP',HINFP_sum,HINFP_length,HINFP_mean,HINFP_median,HINFP_quantile_95,HINFP_quantile]
    row_counter=row_counter+1
del HINFP
plt.figure(figsize=(20, 17))


PPARA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PPARA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PPARA.columns=['TF_TG_SCORE','label']
PPARA_sum = sum(PPARA['TF_TG_SCORE'])
PPARA_length = len(PPARA)
PPARA_mean=0
PPARA_quantile=0
PPARA_quantile_95=0
PPARA_median=0
PPARA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PPARA['TF_TG_SCORE'])
if(PPARA_length>0):
    PPARA_mean=PPARA_sum/PPARA_length
    PPARA_quantile=np.percentile(PPARA['TF_TG_SCORE'], 99)
    PPARA_quantile_95=np.percentile(PPARA['TF_TG_SCORE'], 95)
    PPARA_median=sts.median(PPARA['TF_TG_SCORE'])
if(PPARA_median > background_median and PPARA_mannwhitneyU['pvalue']<0.01):
    background_PPARA = pd.concat([background,PPARA],axis=0)
    ax_PPARA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PPARA,palette="Set3")
    ax_PPARA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PPARA.png')
    del background_PPARA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PPARA',PPARA_sum,PPARA_length,PPARA_mean,PPARA_median,PPARA_quantile_95,PPARA_quantile]
    row_counter=row_counter+1
del PPARA
plt.figure(figsize=(20, 17))


EPAS1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/EPAS1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
EPAS1.columns=['TF_TG_SCORE','label']
EPAS1_sum = sum(EPAS1['TF_TG_SCORE'])
EPAS1_length = len(EPAS1)
EPAS1_mean=0
EPAS1_quantile=0
EPAS1_quantile_95=0
EPAS1_median=0
EPAS1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],EPAS1['TF_TG_SCORE'])
if(EPAS1_length>0):
    EPAS1_mean=EPAS1_sum/EPAS1_length
    EPAS1_quantile=np.percentile(EPAS1['TF_TG_SCORE'], 99)
    EPAS1_quantile_95=np.percentile(EPAS1['TF_TG_SCORE'], 95)
    EPAS1_median=sts.median(EPAS1['TF_TG_SCORE'])
if(EPAS1_median > background_median and EPAS1_mannwhitneyU['pvalue']<0.01):
    background_EPAS1 = pd.concat([background,EPAS1],axis=0)
    ax_EPAS1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_EPAS1,palette="Set3")
    ax_EPAS1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/EPAS1.png')
    del background_EPAS1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['EPAS1',EPAS1_sum,EPAS1_length,EPAS1_mean,EPAS1_median,EPAS1_quantile_95,EPAS1_quantile]
    row_counter=row_counter+1
del EPAS1
plt.figure(figsize=(20, 17))


MBD2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MBD2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MBD2.columns=['TF_TG_SCORE','label']
MBD2_sum = sum(MBD2['TF_TG_SCORE'])
MBD2_length = len(MBD2)
MBD2_mean=0
MBD2_quantile=0
MBD2_quantile_95=0
MBD2_median=0
MBD2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MBD2['TF_TG_SCORE'])
if(MBD2_length>0):
    MBD2_mean=MBD2_sum/MBD2_length
    MBD2_quantile=np.percentile(MBD2['TF_TG_SCORE'], 99)
    MBD2_quantile_95=np.percentile(MBD2['TF_TG_SCORE'], 95)
    MBD2_median=sts.median(MBD2['TF_TG_SCORE'])
if(MBD2_median > background_median and MBD2_mannwhitneyU['pvalue']<0.01):
    background_MBD2 = pd.concat([background,MBD2],axis=0)
    ax_MBD2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MBD2,palette="Set3")
    ax_MBD2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MBD2.png')
    del background_MBD2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MBD2',MBD2_sum,MBD2_length,MBD2_mean,MBD2_median,MBD2_quantile_95,MBD2_quantile]
    row_counter=row_counter+1
del MBD2
plt.figure(figsize=(20, 17))


KLF6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/KLF6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
KLF6.columns=['TF_TG_SCORE','label']
KLF6_sum = sum(KLF6['TF_TG_SCORE'])
KLF6_length = len(KLF6)
KLF6_mean=0
KLF6_quantile=0
KLF6_quantile_95=0
KLF6_median=0
KLF6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],KLF6['TF_TG_SCORE'])
if(KLF6_length>0):
    KLF6_mean=KLF6_sum/KLF6_length
    KLF6_quantile=np.percentile(KLF6['TF_TG_SCORE'], 99)
    KLF6_quantile_95=np.percentile(KLF6['TF_TG_SCORE'], 95)
    KLF6_median=sts.median(KLF6['TF_TG_SCORE'])
if(KLF6_median > background_median and KLF6_mannwhitneyU['pvalue']<0.01):
    background_KLF6 = pd.concat([background,KLF6],axis=0)
    ax_KLF6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_KLF6,palette="Set3")
    ax_KLF6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/KLF6.png')
    del background_KLF6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['KLF6',KLF6_sum,KLF6_length,KLF6_mean,KLF6_median,KLF6_quantile_95,KLF6_quantile]
    row_counter=row_counter+1
del KLF6
plt.figure(figsize=(20, 17))


GRHL2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/GRHL2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
GRHL2.columns=['TF_TG_SCORE','label']
GRHL2_sum = sum(GRHL2['TF_TG_SCORE'])
GRHL2_length = len(GRHL2)
GRHL2_mean=0
GRHL2_quantile=0
GRHL2_quantile_95=0
GRHL2_median=0
GRHL2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],GRHL2['TF_TG_SCORE'])
if(GRHL2_length>0):
    GRHL2_mean=GRHL2_sum/GRHL2_length
    GRHL2_quantile=np.percentile(GRHL2['TF_TG_SCORE'], 99)
    GRHL2_quantile_95=np.percentile(GRHL2['TF_TG_SCORE'], 95)
    GRHL2_median=sts.median(GRHL2['TF_TG_SCORE'])
if(GRHL2_median > background_median and GRHL2_mannwhitneyU['pvalue']<0.01):
    background_GRHL2 = pd.concat([background,GRHL2],axis=0)
    ax_GRHL2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_GRHL2,palette="Set3")
    ax_GRHL2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/GRHL2.png')
    del background_GRHL2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['GRHL2',GRHL2_sum,GRHL2_length,GRHL2_mean,GRHL2_median,GRHL2_quantile_95,GRHL2_quantile]
    row_counter=row_counter+1
del GRHL2
plt.figure(figsize=(20, 17))


CEBPA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CEBPA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CEBPA.columns=['TF_TG_SCORE','label']
CEBPA_sum = sum(CEBPA['TF_TG_SCORE'])
CEBPA_length = len(CEBPA)
CEBPA_mean=0
CEBPA_quantile=0
CEBPA_quantile_95=0
CEBPA_median=0
CEBPA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CEBPA['TF_TG_SCORE'])
if(CEBPA_length>0):
    CEBPA_mean=CEBPA_sum/CEBPA_length
    CEBPA_quantile=np.percentile(CEBPA['TF_TG_SCORE'], 99)
    CEBPA_quantile_95=np.percentile(CEBPA['TF_TG_SCORE'], 95)
    CEBPA_median=sts.median(CEBPA['TF_TG_SCORE'])
if(CEBPA_median > background_median and CEBPA_mannwhitneyU['pvalue']<0.01):
    background_CEBPA = pd.concat([background,CEBPA],axis=0)
    ax_CEBPA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CEBPA,palette="Set3")
    ax_CEBPA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CEBPA.png')
    del background_CEBPA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CEBPA',CEBPA_sum,CEBPA_length,CEBPA_mean,CEBPA_median,CEBPA_quantile_95,CEBPA_quantile]
    row_counter=row_counter+1
del CEBPA
plt.figure(figsize=(20, 17))


NFIB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFIB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFIB.columns=['TF_TG_SCORE','label']
NFIB_sum = sum(NFIB['TF_TG_SCORE'])
NFIB_length = len(NFIB)
NFIB_mean=0
NFIB_quantile=0
NFIB_quantile_95=0
NFIB_median=0
NFIB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFIB['TF_TG_SCORE'])
if(NFIB_length>0):
    NFIB_mean=NFIB_sum/NFIB_length
    NFIB_quantile=np.percentile(NFIB['TF_TG_SCORE'], 99)
    NFIB_quantile_95=np.percentile(NFIB['TF_TG_SCORE'], 95)
    NFIB_median=sts.median(NFIB['TF_TG_SCORE'])
if(NFIB_median > background_median and NFIB_mannwhitneyU['pvalue']<0.01):
    background_NFIB = pd.concat([background,NFIB],axis=0)
    ax_NFIB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFIB,palette="Set3")
    ax_NFIB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFIB.png')
    del background_NFIB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFIB',NFIB_sum,NFIB_length,NFIB_mean,NFIB_median,NFIB_quantile_95,NFIB_quantile]
    row_counter=row_counter+1
del NFIB
plt.figure(figsize=(20, 17))


HDAC2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HDAC2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HDAC2.columns=['TF_TG_SCORE','label']
HDAC2_sum = sum(HDAC2['TF_TG_SCORE'])
HDAC2_length = len(HDAC2)
HDAC2_mean=0
HDAC2_quantile=0
HDAC2_quantile_95=0
HDAC2_median=0
HDAC2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HDAC2['TF_TG_SCORE'])
if(HDAC2_length>0):
    HDAC2_mean=HDAC2_sum/HDAC2_length
    HDAC2_quantile=np.percentile(HDAC2['TF_TG_SCORE'], 99)
    HDAC2_quantile_95=np.percentile(HDAC2['TF_TG_SCORE'], 95)
    HDAC2_median=sts.median(HDAC2['TF_TG_SCORE'])
if(HDAC2_median > background_median and HDAC2_mannwhitneyU['pvalue']<0.01):
    background_HDAC2 = pd.concat([background,HDAC2],axis=0)
    ax_HDAC2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HDAC2,palette="Set3")
    ax_HDAC2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HDAC2.png')
    del background_HDAC2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HDAC2',HDAC2_sum,HDAC2_length,HDAC2_mean,HDAC2_median,HDAC2_quantile_95,HDAC2_quantile]
    row_counter=row_counter+1
del HDAC2
plt.figure(figsize=(20, 17))


SMAD4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SMAD4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SMAD4.columns=['TF_TG_SCORE','label']
SMAD4_sum = sum(SMAD4['TF_TG_SCORE'])
SMAD4_length = len(SMAD4)
SMAD4_mean=0
SMAD4_quantile=0
SMAD4_quantile_95=0
SMAD4_median=0
SMAD4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SMAD4['TF_TG_SCORE'])
if(SMAD4_length>0):
    SMAD4_mean=SMAD4_sum/SMAD4_length
    SMAD4_quantile=np.percentile(SMAD4['TF_TG_SCORE'], 99)
    SMAD4_quantile_95=np.percentile(SMAD4['TF_TG_SCORE'], 95)
    SMAD4_median=sts.median(SMAD4['TF_TG_SCORE'])
if(SMAD4_median > background_median and SMAD4_mannwhitneyU['pvalue']<0.01):
    background_SMAD4 = pd.concat([background,SMAD4],axis=0)
    ax_SMAD4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SMAD4,palette="Set3")
    ax_SMAD4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SMAD4.png')
    del background_SMAD4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SMAD4',SMAD4_sum,SMAD4_length,SMAD4_mean,SMAD4_median,SMAD4_quantile_95,SMAD4_quantile]
    row_counter=row_counter+1
del SMAD4
plt.figure(figsize=(20, 17))


MAX=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MAX_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MAX.columns=['TF_TG_SCORE','label']
MAX_sum = sum(MAX['TF_TG_SCORE'])
MAX_length = len(MAX)
MAX_mean=0
MAX_quantile=0
MAX_quantile_95=0
MAX_median=0
MAX_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MAX['TF_TG_SCORE'])
if(MAX_length>0):
    MAX_mean=MAX_sum/MAX_length
    MAX_quantile=np.percentile(MAX['TF_TG_SCORE'], 99)
    MAX_quantile_95=np.percentile(MAX['TF_TG_SCORE'], 95)
    MAX_median=sts.median(MAX['TF_TG_SCORE'])
if(MAX_median > background_median and MAX_mannwhitneyU['pvalue']<0.01):
    background_MAX = pd.concat([background,MAX],axis=0)
    ax_MAX = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MAX,palette="Set3")
    ax_MAX.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MAX.png')
    del background_MAX
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MAX',MAX_sum,MAX_length,MAX_mean,MAX_median,MAX_quantile_95,MAX_quantile]
    row_counter=row_counter+1
del MAX
plt.figure(figsize=(20, 17))


KLF3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/KLF3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
KLF3.columns=['TF_TG_SCORE','label']
KLF3_sum = sum(KLF3['TF_TG_SCORE'])
KLF3_length = len(KLF3)
KLF3_mean=0
KLF3_quantile=0
KLF3_quantile_95=0
KLF3_median=0
KLF3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],KLF3['TF_TG_SCORE'])
if(KLF3_length>0):
    KLF3_mean=KLF3_sum/KLF3_length
    KLF3_quantile=np.percentile(KLF3['TF_TG_SCORE'], 99)
    KLF3_quantile_95=np.percentile(KLF3['TF_TG_SCORE'], 95)
    KLF3_median=sts.median(KLF3['TF_TG_SCORE'])
if(KLF3_median > background_median and KLF3_mannwhitneyU['pvalue']<0.01):
    background_KLF3 = pd.concat([background,KLF3],axis=0)
    ax_KLF3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_KLF3,palette="Set3")
    ax_KLF3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/KLF3.png')
    del background_KLF3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['KLF3',KLF3_sum,KLF3_length,KLF3_mean,KLF3_median,KLF3_quantile_95,KLF3_quantile]
    row_counter=row_counter+1
del KLF3
plt.figure(figsize=(20, 17))


BACH1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BACH1..MAFK_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BACH1.columns=['TF_TG_SCORE','label']
BACH1_sum = sum(BACH1['TF_TG_SCORE'])
BACH1_length = len(BACH1)
BACH1_mean=0
BACH1_quantile=0
BACH1_quantile_95=0
BACH1_median=0
BACH1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BACH1['TF_TG_SCORE'])
if(BACH1_length>0):
    BACH1_mean=BACH1_sum/BACH1_length
    BACH1_quantile=np.percentile(BACH1['TF_TG_SCORE'], 99)
    BACH1_quantile_95=np.percentile(BACH1['TF_TG_SCORE'], 95)
    BACH1_median=sts.median(BACH1['TF_TG_SCORE'])
if(BACH1_median > background_median and BACH1_mannwhitneyU['pvalue']<0.01):
    background_BACH1 = pd.concat([background,BACH1],axis=0)
    ax_BACH1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BACH1,palette="Set3")
    ax_BACH1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BACH1..MAFK.png')
    del background_BACH1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BACH1..MAFK',BACH1_sum,BACH1_length,BACH1_mean,BACH1_median,BACH1_quantile_95,BACH1_quantile]
    row_counter=row_counter+1
del BACH1
plt.figure(figsize=(20, 17))


RORA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RORA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RORA.columns=['TF_TG_SCORE','label']
RORA_sum = sum(RORA['TF_TG_SCORE'])
RORA_length = len(RORA)
RORA_mean=0
RORA_quantile=0
RORA_quantile_95=0
RORA_median=0
RORA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RORA['TF_TG_SCORE'])
if(RORA_length>0):
    RORA_mean=RORA_sum/RORA_length
    RORA_quantile=np.percentile(RORA['TF_TG_SCORE'], 99)
    RORA_quantile_95=np.percentile(RORA['TF_TG_SCORE'], 95)
    RORA_median=sts.median(RORA['TF_TG_SCORE'])
if(RORA_median > background_median and RORA_mannwhitneyU['pvalue']<0.01):
    background_RORA = pd.concat([background,RORA],axis=0)
    ax_RORA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RORA,palette="Set3")
    ax_RORA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RORA.png')
    del background_RORA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RORA',RORA_sum,RORA_length,RORA_mean,RORA_median,RORA_quantile_95,RORA_quantile]
    row_counter=row_counter+1
del RORA
plt.figure(figsize=(20, 17))


MAFB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MAFB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MAFB.columns=['TF_TG_SCORE','label']
MAFB_sum = sum(MAFB['TF_TG_SCORE'])
MAFB_length = len(MAFB)
MAFB_mean=0
MAFB_quantile=0
MAFB_quantile_95=0
MAFB_median=0
MAFB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MAFB['TF_TG_SCORE'])
if(MAFB_length>0):
    MAFB_mean=MAFB_sum/MAFB_length
    MAFB_quantile=np.percentile(MAFB['TF_TG_SCORE'], 99)
    MAFB_quantile_95=np.percentile(MAFB['TF_TG_SCORE'], 95)
    MAFB_median=sts.median(MAFB['TF_TG_SCORE'])
if(MAFB_median > background_median and MAFB_mannwhitneyU['pvalue']<0.01):
    background_MAFB = pd.concat([background,MAFB],axis=0)
    ax_MAFB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MAFB,palette="Set3")
    ax_MAFB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MAFB.png')
    del background_MAFB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MAFB',MAFB_sum,MAFB_length,MAFB_mean,MAFB_median,MAFB_quantile_95,MAFB_quantile]
    row_counter=row_counter+1
del MAFB
plt.figure(figsize=(20, 17))


THAP1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/THAP1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
THAP1.columns=['TF_TG_SCORE','label']
THAP1_sum = sum(THAP1['TF_TG_SCORE'])
THAP1_length = len(THAP1)
THAP1_mean=0
THAP1_quantile=0
THAP1_quantile_95=0
THAP1_median=0
THAP1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],THAP1['TF_TG_SCORE'])
if(THAP1_length>0):
    THAP1_mean=THAP1_sum/THAP1_length
    THAP1_quantile=np.percentile(THAP1['TF_TG_SCORE'], 99)
    THAP1_quantile_95=np.percentile(THAP1['TF_TG_SCORE'], 95)
    THAP1_median=sts.median(THAP1['TF_TG_SCORE'])
if(THAP1_median > background_median and THAP1_mannwhitneyU['pvalue']<0.01):
    background_THAP1 = pd.concat([background,THAP1],axis=0)
    ax_THAP1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_THAP1,palette="Set3")
    ax_THAP1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/THAP1.png')
    del background_THAP1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['THAP1',THAP1_sum,THAP1_length,THAP1_mean,THAP1_median,THAP1_quantile_95,THAP1_quantile]
    row_counter=row_counter+1
del THAP1
plt.figure(figsize=(20, 17))


GMEB1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/GMEB1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
GMEB1.columns=['TF_TG_SCORE','label']
GMEB1_sum = sum(GMEB1['TF_TG_SCORE'])
GMEB1_length = len(GMEB1)
GMEB1_mean=0
GMEB1_quantile=0
GMEB1_quantile_95=0
GMEB1_median=0
GMEB1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],GMEB1['TF_TG_SCORE'])
if(GMEB1_length>0):
    GMEB1_mean=GMEB1_sum/GMEB1_length
    GMEB1_quantile=np.percentile(GMEB1['TF_TG_SCORE'], 99)
    GMEB1_quantile_95=np.percentile(GMEB1['TF_TG_SCORE'], 95)
    GMEB1_median=sts.median(GMEB1['TF_TG_SCORE'])
if(GMEB1_median > background_median and GMEB1_mannwhitneyU['pvalue']<0.01):
    background_GMEB1 = pd.concat([background,GMEB1],axis=0)
    ax_GMEB1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_GMEB1,palette="Set3")
    ax_GMEB1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/GMEB1.png')
    del background_GMEB1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['GMEB1',GMEB1_sum,GMEB1_length,GMEB1_mean,GMEB1_median,GMEB1_quantile_95,GMEB1_quantile]
    row_counter=row_counter+1
del GMEB1
plt.figure(figsize=(20, 17))


TRIM28=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TRIM28_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TRIM28.columns=['TF_TG_SCORE','label']
TRIM28_sum = sum(TRIM28['TF_TG_SCORE'])
TRIM28_length = len(TRIM28)
TRIM28_mean=0
TRIM28_quantile=0
TRIM28_quantile_95=0
TRIM28_median=0
TRIM28_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TRIM28['TF_TG_SCORE'])
if(TRIM28_length>0):
    TRIM28_mean=TRIM28_sum/TRIM28_length
    TRIM28_quantile=np.percentile(TRIM28['TF_TG_SCORE'], 99)
    TRIM28_quantile_95=np.percentile(TRIM28['TF_TG_SCORE'], 95)
    TRIM28_median=sts.median(TRIM28['TF_TG_SCORE'])
if(TRIM28_median > background_median and TRIM28_mannwhitneyU['pvalue']<0.01):
    background_TRIM28 = pd.concat([background,TRIM28],axis=0)
    ax_TRIM28 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TRIM28,palette="Set3")
    ax_TRIM28.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TRIM28.png')
    del background_TRIM28
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TRIM28',TRIM28_sum,TRIM28_length,TRIM28_mean,TRIM28_median,TRIM28_quantile_95,TRIM28_quantile]
    row_counter=row_counter+1
del TRIM28
plt.figure(figsize=(20, 17))


ELK4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ELK4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ELK4.columns=['TF_TG_SCORE','label']
ELK4_sum = sum(ELK4['TF_TG_SCORE'])
ELK4_length = len(ELK4)
ELK4_mean=0
ELK4_quantile=0
ELK4_quantile_95=0
ELK4_median=0
ELK4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ELK4['TF_TG_SCORE'])
if(ELK4_length>0):
    ELK4_mean=ELK4_sum/ELK4_length
    ELK4_quantile=np.percentile(ELK4['TF_TG_SCORE'], 99)
    ELK4_quantile_95=np.percentile(ELK4['TF_TG_SCORE'], 95)
    ELK4_median=sts.median(ELK4['TF_TG_SCORE'])
if(ELK4_median > background_median and ELK4_mannwhitneyU['pvalue']<0.01):
    background_ELK4 = pd.concat([background,ELK4],axis=0)
    ax_ELK4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ELK4,palette="Set3")
    ax_ELK4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ELK4.png')
    del background_ELK4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ELK4',ELK4_sum,ELK4_length,ELK4_mean,ELK4_median,ELK4_quantile_95,ELK4_quantile]
    row_counter=row_counter+1
del ELK4
plt.figure(figsize=(20, 17))


TEAD1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TEAD1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TEAD1.columns=['TF_TG_SCORE','label']
TEAD1_sum = sum(TEAD1['TF_TG_SCORE'])
TEAD1_length = len(TEAD1)
TEAD1_mean=0
TEAD1_quantile=0
TEAD1_quantile_95=0
TEAD1_median=0
TEAD1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TEAD1['TF_TG_SCORE'])
if(TEAD1_length>0):
    TEAD1_mean=TEAD1_sum/TEAD1_length
    TEAD1_quantile=np.percentile(TEAD1['TF_TG_SCORE'], 99)
    TEAD1_quantile_95=np.percentile(TEAD1['TF_TG_SCORE'], 95)
    TEAD1_median=sts.median(TEAD1['TF_TG_SCORE'])
if(TEAD1_median > background_median and TEAD1_mannwhitneyU['pvalue']<0.01):
    background_TEAD1 = pd.concat([background,TEAD1],axis=0)
    ax_TEAD1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TEAD1,palette="Set3")
    ax_TEAD1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TEAD1.png')
    del background_TEAD1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TEAD1',TEAD1_sum,TEAD1_length,TEAD1_mean,TEAD1_median,TEAD1_quantile_95,TEAD1_quantile]
    row_counter=row_counter+1
del TEAD1
plt.figure(figsize=(20, 17))


E2F2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/E2F2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
E2F2.columns=['TF_TG_SCORE','label']
E2F2_sum = sum(E2F2['TF_TG_SCORE'])
E2F2_length = len(E2F2)
E2F2_mean=0
E2F2_quantile=0
E2F2_quantile_95=0
E2F2_median=0
E2F2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],E2F2['TF_TG_SCORE'])
if(E2F2_length>0):
    E2F2_mean=E2F2_sum/E2F2_length
    E2F2_quantile=np.percentile(E2F2['TF_TG_SCORE'], 99)
    E2F2_quantile_95=np.percentile(E2F2['TF_TG_SCORE'], 95)
    E2F2_median=sts.median(E2F2['TF_TG_SCORE'])
if(E2F2_median > background_median and E2F2_mannwhitneyU['pvalue']<0.01):
    background_E2F2 = pd.concat([background,E2F2],axis=0)
    ax_E2F2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_E2F2,palette="Set3")
    ax_E2F2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/E2F2.png')
    del background_E2F2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['E2F2',E2F2_sum,E2F2_length,E2F2_mean,E2F2_median,E2F2_quantile_95,E2F2_quantile]
    row_counter=row_counter+1
del E2F2
plt.figure(figsize=(20, 17))


HOXD9=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HOXD9_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HOXD9.columns=['TF_TG_SCORE','label']
HOXD9_sum = sum(HOXD9['TF_TG_SCORE'])
HOXD9_length = len(HOXD9)
HOXD9_mean=0
HOXD9_quantile=0
HOXD9_quantile_95=0
HOXD9_median=0
HOXD9_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HOXD9['TF_TG_SCORE'])
if(HOXD9_length>0):
    HOXD9_mean=HOXD9_sum/HOXD9_length
    HOXD9_quantile=np.percentile(HOXD9['TF_TG_SCORE'], 99)
    HOXD9_quantile_95=np.percentile(HOXD9['TF_TG_SCORE'], 95)
    HOXD9_median=sts.median(HOXD9['TF_TG_SCORE'])
if(HOXD9_median > background_median and HOXD9_mannwhitneyU['pvalue']<0.01):
    background_HOXD9 = pd.concat([background,HOXD9],axis=0)
    ax_HOXD9 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HOXD9,palette="Set3")
    ax_HOXD9.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HOXD9.png')
    del background_HOXD9
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HOXD9',HOXD9_sum,HOXD9_length,HOXD9_mean,HOXD9_median,HOXD9_quantile_95,HOXD9_quantile]
    row_counter=row_counter+1
del HOXD9
plt.figure(figsize=(20, 17))


NPAS2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NPAS2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NPAS2.columns=['TF_TG_SCORE','label']
NPAS2_sum = sum(NPAS2['TF_TG_SCORE'])
NPAS2_length = len(NPAS2)
NPAS2_mean=0
NPAS2_quantile=0
NPAS2_quantile_95=0
NPAS2_median=0
NPAS2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NPAS2['TF_TG_SCORE'])
if(NPAS2_length>0):
    NPAS2_mean=NPAS2_sum/NPAS2_length
    NPAS2_quantile=np.percentile(NPAS2['TF_TG_SCORE'], 99)
    NPAS2_quantile_95=np.percentile(NPAS2['TF_TG_SCORE'], 95)
    NPAS2_median=sts.median(NPAS2['TF_TG_SCORE'])
if(NPAS2_median > background_median and NPAS2_mannwhitneyU['pvalue']<0.01):
    background_NPAS2 = pd.concat([background,NPAS2],axis=0)
    ax_NPAS2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NPAS2,palette="Set3")
    ax_NPAS2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NPAS2.png')
    del background_NPAS2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NPAS2',NPAS2_sum,NPAS2_length,NPAS2_mean,NPAS2_median,NPAS2_quantile_95,NPAS2_quantile]
    row_counter=row_counter+1
del NPAS2
plt.figure(figsize=(20, 17))


ESRRA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ESRRA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ESRRA.columns=['TF_TG_SCORE','label']
ESRRA_sum = sum(ESRRA['TF_TG_SCORE'])
ESRRA_length = len(ESRRA)
ESRRA_mean=0
ESRRA_quantile=0
ESRRA_quantile_95=0
ESRRA_median=0
ESRRA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ESRRA['TF_TG_SCORE'])
if(ESRRA_length>0):
    ESRRA_mean=ESRRA_sum/ESRRA_length
    ESRRA_quantile=np.percentile(ESRRA['TF_TG_SCORE'], 99)
    ESRRA_quantile_95=np.percentile(ESRRA['TF_TG_SCORE'], 95)
    ESRRA_median=sts.median(ESRRA['TF_TG_SCORE'])
if(ESRRA_median > background_median and ESRRA_mannwhitneyU['pvalue']<0.01):
    background_ESRRA = pd.concat([background,ESRRA],axis=0)
    ax_ESRRA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ESRRA,palette="Set3")
    ax_ESRRA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ESRRA.png')
    del background_ESRRA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ESRRA',ESRRA_sum,ESRRA_length,ESRRA_mean,ESRRA_median,ESRRA_quantile_95,ESRRA_quantile]
    row_counter=row_counter+1
del ESRRA
plt.figure(figsize=(20, 17))


ZFX=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ZFX_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ZFX.columns=['TF_TG_SCORE','label']
ZFX_sum = sum(ZFX['TF_TG_SCORE'])
ZFX_length = len(ZFX)
ZFX_mean=0
ZFX_quantile=0
ZFX_quantile_95=0
ZFX_median=0
ZFX_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ZFX['TF_TG_SCORE'])
if(ZFX_length>0):
    ZFX_mean=ZFX_sum/ZFX_length
    ZFX_quantile=np.percentile(ZFX['TF_TG_SCORE'], 99)
    ZFX_quantile_95=np.percentile(ZFX['TF_TG_SCORE'], 95)
    ZFX_median=sts.median(ZFX['TF_TG_SCORE'])
if(ZFX_median > background_median and ZFX_mannwhitneyU['pvalue']<0.01):
    background_ZFX = pd.concat([background,ZFX],axis=0)
    ax_ZFX = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ZFX,palette="Set3")
    ax_ZFX.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ZFX.png')
    del background_ZFX
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ZFX',ZFX_sum,ZFX_length,ZFX_mean,ZFX_median,ZFX_quantile_95,ZFX_quantile]
    row_counter=row_counter+1
del ZFX
plt.figure(figsize=(20, 17))


ZEB1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ZEB1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ZEB1.columns=['TF_TG_SCORE','label']
ZEB1_sum = sum(ZEB1['TF_TG_SCORE'])
ZEB1_length = len(ZEB1)
ZEB1_mean=0
ZEB1_quantile=0
ZEB1_quantile_95=0
ZEB1_median=0
ZEB1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ZEB1['TF_TG_SCORE'])
if(ZEB1_length>0):
    ZEB1_mean=ZEB1_sum/ZEB1_length
    ZEB1_quantile=np.percentile(ZEB1['TF_TG_SCORE'], 99)
    ZEB1_quantile_95=np.percentile(ZEB1['TF_TG_SCORE'], 95)
    ZEB1_median=sts.median(ZEB1['TF_TG_SCORE'])
if(ZEB1_median > background_median and ZEB1_mannwhitneyU['pvalue']<0.01):
    background_ZEB1 = pd.concat([background,ZEB1],axis=0)
    ax_ZEB1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ZEB1,palette="Set3")
    ax_ZEB1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ZEB1.png')
    del background_ZEB1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ZEB1',ZEB1_sum,ZEB1_length,ZEB1_mean,ZEB1_median,ZEB1_quantile_95,ZEB1_quantile]
    row_counter=row_counter+1
del ZEB1
plt.figure(figsize=(20, 17))


SP1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SP1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SP1.columns=['TF_TG_SCORE','label']
SP1_sum = sum(SP1['TF_TG_SCORE'])
SP1_length = len(SP1)
SP1_mean=0
SP1_quantile=0
SP1_quantile_95=0
SP1_median=0
SP1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SP1['TF_TG_SCORE'])
if(SP1_length>0):
    SP1_mean=SP1_sum/SP1_length
    SP1_quantile=np.percentile(SP1['TF_TG_SCORE'], 99)
    SP1_quantile_95=np.percentile(SP1['TF_TG_SCORE'], 95)
    SP1_median=sts.median(SP1['TF_TG_SCORE'])
if(SP1_median > background_median and SP1_mannwhitneyU['pvalue']<0.01):
    background_SP1 = pd.concat([background,SP1],axis=0)
    ax_SP1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SP1,palette="Set3")
    ax_SP1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SP1.png')
    del background_SP1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SP1',SP1_sum,SP1_length,SP1_mean,SP1_median,SP1_quantile_95,SP1_quantile]
    row_counter=row_counter+1
del SP1
plt.figure(figsize=(20, 17))


BRCA1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BRCA1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BRCA1.columns=['TF_TG_SCORE','label']
BRCA1_sum = sum(BRCA1['TF_TG_SCORE'])
BRCA1_length = len(BRCA1)
BRCA1_mean=0
BRCA1_quantile=0
BRCA1_quantile_95=0
BRCA1_median=0
BRCA1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BRCA1['TF_TG_SCORE'])
if(BRCA1_length>0):
    BRCA1_mean=BRCA1_sum/BRCA1_length
    BRCA1_quantile=np.percentile(BRCA1['TF_TG_SCORE'], 99)
    BRCA1_quantile_95=np.percentile(BRCA1['TF_TG_SCORE'], 95)
    BRCA1_median=sts.median(BRCA1['TF_TG_SCORE'])
if(BRCA1_median > background_median and BRCA1_mannwhitneyU['pvalue']<0.01):
    background_BRCA1 = pd.concat([background,BRCA1],axis=0)
    ax_BRCA1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BRCA1,palette="Set3")
    ax_BRCA1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BRCA1.png')
    del background_BRCA1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BRCA1',BRCA1_sum,BRCA1_length,BRCA1_mean,BRCA1_median,BRCA1_quantile_95,BRCA1_quantile]
    row_counter=row_counter+1
del BRCA1
plt.figure(figsize=(20, 17))


CHD2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CHD2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CHD2.columns=['TF_TG_SCORE','label']
CHD2_sum = sum(CHD2['TF_TG_SCORE'])
CHD2_length = len(CHD2)
CHD2_mean=0
CHD2_quantile=0
CHD2_quantile_95=0
CHD2_median=0
CHD2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CHD2['TF_TG_SCORE'])
if(CHD2_length>0):
    CHD2_mean=CHD2_sum/CHD2_length
    CHD2_quantile=np.percentile(CHD2['TF_TG_SCORE'], 99)
    CHD2_quantile_95=np.percentile(CHD2['TF_TG_SCORE'], 95)
    CHD2_median=sts.median(CHD2['TF_TG_SCORE'])
if(CHD2_median > background_median and CHD2_mannwhitneyU['pvalue']<0.01):
    background_CHD2 = pd.concat([background,CHD2],axis=0)
    ax_CHD2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CHD2,palette="Set3")
    ax_CHD2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CHD2.png')
    del background_CHD2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CHD2',CHD2_sum,CHD2_length,CHD2_mean,CHD2_median,CHD2_quantile_95,CHD2_quantile]
    row_counter=row_counter+1
del CHD2
plt.figure(figsize=(20, 17))


SRF=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SRF_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SRF.columns=['TF_TG_SCORE','label']
SRF_sum = sum(SRF['TF_TG_SCORE'])
SRF_length = len(SRF)
SRF_mean=0
SRF_quantile=0
SRF_quantile_95=0
SRF_median=0
SRF_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SRF['TF_TG_SCORE'])
if(SRF_length>0):
    SRF_mean=SRF_sum/SRF_length
    SRF_quantile=np.percentile(SRF['TF_TG_SCORE'], 99)
    SRF_quantile_95=np.percentile(SRF['TF_TG_SCORE'], 95)
    SRF_median=sts.median(SRF['TF_TG_SCORE'])
if(SRF_median > background_median and SRF_mannwhitneyU['pvalue']<0.01):
    background_SRF = pd.concat([background,SRF],axis=0)
    ax_SRF = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SRF,palette="Set3")
    ax_SRF.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SRF.png')
    del background_SRF
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SRF',SRF_sum,SRF_length,SRF_mean,SRF_median,SRF_quantile_95,SRF_quantile]
    row_counter=row_counter+1
del SRF
plt.figure(figsize=(20, 17))


ID2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ID2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ID2.columns=['TF_TG_SCORE','label']
ID2_sum = sum(ID2['TF_TG_SCORE'])
ID2_length = len(ID2)
ID2_mean=0
ID2_quantile=0
ID2_quantile_95=0
ID2_median=0
ID2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ID2['TF_TG_SCORE'])
if(ID2_length>0):
    ID2_mean=ID2_sum/ID2_length
    ID2_quantile=np.percentile(ID2['TF_TG_SCORE'], 99)
    ID2_quantile_95=np.percentile(ID2['TF_TG_SCORE'], 95)
    ID2_median=sts.median(ID2['TF_TG_SCORE'])
if(ID2_median > background_median and ID2_mannwhitneyU['pvalue']<0.01):
    background_ID2 = pd.concat([background,ID2],axis=0)
    ax_ID2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ID2,palette="Set3")
    ax_ID2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ID2.png')
    del background_ID2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ID2',ID2_sum,ID2_length,ID2_mean,ID2_median,ID2_quantile_95,ID2_quantile]
    row_counter=row_counter+1
del ID2
plt.figure(figsize=(20, 17))


SMC3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SMC3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SMC3.columns=['TF_TG_SCORE','label']
SMC3_sum = sum(SMC3['TF_TG_SCORE'])
SMC3_length = len(SMC3)
SMC3_mean=0
SMC3_quantile=0
SMC3_quantile_95=0
SMC3_median=0
SMC3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SMC3['TF_TG_SCORE'])
if(SMC3_length>0):
    SMC3_mean=SMC3_sum/SMC3_length
    SMC3_quantile=np.percentile(SMC3['TF_TG_SCORE'], 99)
    SMC3_quantile_95=np.percentile(SMC3['TF_TG_SCORE'], 95)
    SMC3_median=sts.median(SMC3['TF_TG_SCORE'])
if(SMC3_median > background_median and SMC3_mannwhitneyU['pvalue']<0.01):
    background_SMC3 = pd.concat([background,SMC3],axis=0)
    ax_SMC3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SMC3,palette="Set3")
    ax_SMC3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SMC3.png')
    del background_SMC3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SMC3',SMC3_sum,SMC3_length,SMC3_mean,SMC3_median,SMC3_quantile_95,SMC3_quantile]
    row_counter=row_counter+1
del SMC3
plt.figure(figsize=(20, 17))


ESR1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ESR1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ESR1.columns=['TF_TG_SCORE','label']
ESR1_sum = sum(ESR1['TF_TG_SCORE'])
ESR1_length = len(ESR1)
ESR1_mean=0
ESR1_quantile=0
ESR1_quantile_95=0
ESR1_median=0
ESR1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ESR1['TF_TG_SCORE'])
if(ESR1_length>0):
    ESR1_mean=ESR1_sum/ESR1_length
    ESR1_quantile=np.percentile(ESR1['TF_TG_SCORE'], 99)
    ESR1_quantile_95=np.percentile(ESR1['TF_TG_SCORE'], 95)
    ESR1_median=sts.median(ESR1['TF_TG_SCORE'])
if(ESR1_median > background_median and ESR1_mannwhitneyU['pvalue']<0.01):
    background_ESR1 = pd.concat([background,ESR1],axis=0)
    ax_ESR1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ESR1,palette="Set3")
    ax_ESR1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ESR1.png')
    del background_ESR1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ESR1',ESR1_sum,ESR1_length,ESR1_mean,ESR1_median,ESR1_quantile_95,ESR1_quantile]
    row_counter=row_counter+1
del ESR1
plt.figure(figsize=(20, 17))


LYL1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/LYL1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
LYL1.columns=['TF_TG_SCORE','label']
LYL1_sum = sum(LYL1['TF_TG_SCORE'])
LYL1_length = len(LYL1)
LYL1_mean=0
LYL1_quantile=0
LYL1_quantile_95=0
LYL1_median=0
LYL1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],LYL1['TF_TG_SCORE'])
if(LYL1_length>0):
    LYL1_mean=LYL1_sum/LYL1_length
    LYL1_quantile=np.percentile(LYL1['TF_TG_SCORE'], 99)
    LYL1_quantile_95=np.percentile(LYL1['TF_TG_SCORE'], 95)
    LYL1_median=sts.median(LYL1['TF_TG_SCORE'])
if(LYL1_median > background_median and LYL1_mannwhitneyU['pvalue']<0.01):
    background_LYL1 = pd.concat([background,LYL1],axis=0)
    ax_LYL1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_LYL1,palette="Set3")
    ax_LYL1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/LYL1.png')
    del background_LYL1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['LYL1',LYL1_sum,LYL1_length,LYL1_mean,LYL1_median,LYL1_quantile_95,LYL1_quantile]
    row_counter=row_counter+1
del LYL1
plt.figure(figsize=(20, 17))


CREM=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CREM_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CREM.columns=['TF_TG_SCORE','label']
CREM_sum = sum(CREM['TF_TG_SCORE'])
CREM_length = len(CREM)
CREM_mean=0
CREM_quantile=0
CREM_quantile_95=0
CREM_median=0
CREM_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CREM['TF_TG_SCORE'])
if(CREM_length>0):
    CREM_mean=CREM_sum/CREM_length
    CREM_quantile=np.percentile(CREM['TF_TG_SCORE'], 99)
    CREM_quantile_95=np.percentile(CREM['TF_TG_SCORE'], 95)
    CREM_median=sts.median(CREM['TF_TG_SCORE'])
if(CREM_median > background_median and CREM_mannwhitneyU['pvalue']<0.01):
    background_CREM = pd.concat([background,CREM],axis=0)
    ax_CREM = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CREM,palette="Set3")
    ax_CREM.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CREM.png')
    del background_CREM
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CREM',CREM_sum,CREM_length,CREM_mean,CREM_median,CREM_quantile_95,CREM_quantile]
    row_counter=row_counter+1
del CREM
plt.figure(figsize=(20, 17))


BHLHE40=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BHLHE40_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BHLHE40.columns=['TF_TG_SCORE','label']
BHLHE40_sum = sum(BHLHE40['TF_TG_SCORE'])
BHLHE40_length = len(BHLHE40)
BHLHE40_mean=0
BHLHE40_quantile=0
BHLHE40_quantile_95=0
BHLHE40_median=0
BHLHE40_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BHLHE40['TF_TG_SCORE'])
if(BHLHE40_length>0):
    BHLHE40_mean=BHLHE40_sum/BHLHE40_length
    BHLHE40_quantile=np.percentile(BHLHE40['TF_TG_SCORE'], 99)
    BHLHE40_quantile_95=np.percentile(BHLHE40['TF_TG_SCORE'], 95)
    BHLHE40_median=sts.median(BHLHE40['TF_TG_SCORE'])
if(BHLHE40_median > background_median and BHLHE40_mannwhitneyU['pvalue']<0.01):
    background_BHLHE40 = pd.concat([background,BHLHE40],axis=0)
    ax_BHLHE40 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BHLHE40,palette="Set3")
    ax_BHLHE40.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BHLHE40.png')
    del background_BHLHE40
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BHLHE40',BHLHE40_sum,BHLHE40_length,BHLHE40_mean,BHLHE40_median,BHLHE40_quantile_95,BHLHE40_quantile]
    row_counter=row_counter+1
del BHLHE40
plt.figure(figsize=(20, 17))


REST=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/REST_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
REST.columns=['TF_TG_SCORE','label']
REST_sum = sum(REST['TF_TG_SCORE'])
REST_length = len(REST)
REST_mean=0
REST_quantile=0
REST_quantile_95=0
REST_median=0
REST_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],REST['TF_TG_SCORE'])
if(REST_length>0):
    REST_mean=REST_sum/REST_length
    REST_quantile=np.percentile(REST['TF_TG_SCORE'], 99)
    REST_quantile_95=np.percentile(REST['TF_TG_SCORE'], 95)
    REST_median=sts.median(REST['TF_TG_SCORE'])
if(REST_median > background_median and REST_mannwhitneyU['pvalue']<0.01):
    background_REST = pd.concat([background,REST],axis=0)
    ax_REST = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_REST,palette="Set3")
    ax_REST.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/REST.png')
    del background_REST
    plt.clf()
    df_interesting_stats.loc[row_counter]=['REST',REST_sum,REST_length,REST_mean,REST_median,REST_quantile_95,REST_quantile]
    row_counter=row_counter+1
del REST
plt.figure(figsize=(20, 17))


MXI1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MXI1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MXI1.columns=['TF_TG_SCORE','label']
MXI1_sum = sum(MXI1['TF_TG_SCORE'])
MXI1_length = len(MXI1)
MXI1_mean=0
MXI1_quantile=0
MXI1_quantile_95=0
MXI1_median=0
MXI1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MXI1['TF_TG_SCORE'])
if(MXI1_length>0):
    MXI1_mean=MXI1_sum/MXI1_length
    MXI1_quantile=np.percentile(MXI1['TF_TG_SCORE'], 99)
    MXI1_quantile_95=np.percentile(MXI1['TF_TG_SCORE'], 95)
    MXI1_median=sts.median(MXI1['TF_TG_SCORE'])
if(MXI1_median > background_median and MXI1_mannwhitneyU['pvalue']<0.01):
    background_MXI1 = pd.concat([background,MXI1],axis=0)
    ax_MXI1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MXI1,palette="Set3")
    ax_MXI1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MXI1.png')
    del background_MXI1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MXI1',MXI1_sum,MXI1_length,MXI1_mean,MXI1_median,MXI1_quantile_95,MXI1_quantile]
    row_counter=row_counter+1
del MXI1
plt.figure(figsize=(20, 17))


FOXO3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOXO3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOXO3.columns=['TF_TG_SCORE','label']
FOXO3_sum = sum(FOXO3['TF_TG_SCORE'])
FOXO3_length = len(FOXO3)
FOXO3_mean=0
FOXO3_quantile=0
FOXO3_quantile_95=0
FOXO3_median=0
FOXO3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOXO3['TF_TG_SCORE'])
if(FOXO3_length>0):
    FOXO3_mean=FOXO3_sum/FOXO3_length
    FOXO3_quantile=np.percentile(FOXO3['TF_TG_SCORE'], 99)
    FOXO3_quantile_95=np.percentile(FOXO3['TF_TG_SCORE'], 95)
    FOXO3_median=sts.median(FOXO3['TF_TG_SCORE'])
if(FOXO3_median > background_median and FOXO3_mannwhitneyU['pvalue']<0.01):
    background_FOXO3 = pd.concat([background,FOXO3],axis=0)
    ax_FOXO3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOXO3,palette="Set3")
    ax_FOXO3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOXO3.png')
    del background_FOXO3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOXO3',FOXO3_sum,FOXO3_length,FOXO3_mean,FOXO3_median,FOXO3_quantile_95,FOXO3_quantile]
    row_counter=row_counter+1
del FOXO3
plt.figure(figsize=(20, 17))


USF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/USF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
USF1.columns=['TF_TG_SCORE','label']
USF1_sum = sum(USF1['TF_TG_SCORE'])
USF1_length = len(USF1)
USF1_mean=0
USF1_quantile=0
USF1_quantile_95=0
USF1_median=0
USF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],USF1['TF_TG_SCORE'])
if(USF1_length>0):
    USF1_mean=USF1_sum/USF1_length
    USF1_quantile=np.percentile(USF1['TF_TG_SCORE'], 99)
    USF1_quantile_95=np.percentile(USF1['TF_TG_SCORE'], 95)
    USF1_median=sts.median(USF1['TF_TG_SCORE'])
if(USF1_median > background_median and USF1_mannwhitneyU['pvalue']<0.01):
    background_USF1 = pd.concat([background,USF1],axis=0)
    ax_USF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_USF1,palette="Set3")
    ax_USF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/USF1.png')
    del background_USF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['USF1',USF1_sum,USF1_length,USF1_mean,USF1_median,USF1_quantile_95,USF1_quantile]
    row_counter=row_counter+1
del USF1
plt.figure(figsize=(20, 17))


GABPA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/GABPA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
GABPA.columns=['TF_TG_SCORE','label']
GABPA_sum = sum(GABPA['TF_TG_SCORE'])
GABPA_length = len(GABPA)
GABPA_mean=0
GABPA_quantile=0
GABPA_quantile_95=0
GABPA_median=0
GABPA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],GABPA['TF_TG_SCORE'])
if(GABPA_length>0):
    GABPA_mean=GABPA_sum/GABPA_length
    GABPA_quantile=np.percentile(GABPA['TF_TG_SCORE'], 99)
    GABPA_quantile_95=np.percentile(GABPA['TF_TG_SCORE'], 95)
    GABPA_median=sts.median(GABPA['TF_TG_SCORE'])
if(GABPA_median > background_median and GABPA_mannwhitneyU['pvalue']<0.01):
    background_GABPA = pd.concat([background,GABPA],axis=0)
    ax_GABPA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_GABPA,palette="Set3")
    ax_GABPA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/GABPA.png')
    del background_GABPA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['GABPA',GABPA_sum,GABPA_length,GABPA_mean,GABPA_median,GABPA_quantile_95,GABPA_quantile]
    row_counter=row_counter+1
del GABPA
plt.figure(figsize=(20, 17))


IRF4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF4.columns=['TF_TG_SCORE','label']
IRF4_sum = sum(IRF4['TF_TG_SCORE'])
IRF4_length = len(IRF4)
IRF4_mean=0
IRF4_quantile=0
IRF4_quantile_95=0
IRF4_median=0
IRF4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF4['TF_TG_SCORE'])
if(IRF4_length>0):
    IRF4_mean=IRF4_sum/IRF4_length
    IRF4_quantile=np.percentile(IRF4['TF_TG_SCORE'], 99)
    IRF4_quantile_95=np.percentile(IRF4['TF_TG_SCORE'], 95)
    IRF4_median=sts.median(IRF4['TF_TG_SCORE'])
if(IRF4_median > background_median and IRF4_mannwhitneyU['pvalue']<0.01):
    background_IRF4 = pd.concat([background,IRF4],axis=0)
    ax_IRF4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF4,palette="Set3")
    ax_IRF4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF4.png')
    del background_IRF4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF4',IRF4_sum,IRF4_length,IRF4_mean,IRF4_median,IRF4_quantile_95,IRF4_quantile]
    row_counter=row_counter+1
del IRF4
plt.figure(figsize=(20, 17))


HEY1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HEY1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HEY1.columns=['TF_TG_SCORE','label']
HEY1_sum = sum(HEY1['TF_TG_SCORE'])
HEY1_length = len(HEY1)
HEY1_mean=0
HEY1_quantile=0
HEY1_quantile_95=0
HEY1_median=0
HEY1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HEY1['TF_TG_SCORE'])
if(HEY1_length>0):
    HEY1_mean=HEY1_sum/HEY1_length
    HEY1_quantile=np.percentile(HEY1['TF_TG_SCORE'], 99)
    HEY1_quantile_95=np.percentile(HEY1['TF_TG_SCORE'], 95)
    HEY1_median=sts.median(HEY1['TF_TG_SCORE'])
if(HEY1_median > background_median and HEY1_mannwhitneyU['pvalue']<0.01):
    background_HEY1 = pd.concat([background,HEY1],axis=0)
    ax_HEY1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HEY1,palette="Set3")
    ax_HEY1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HEY1.png')
    del background_HEY1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HEY1',HEY1_sum,HEY1_length,HEY1_mean,HEY1_median,HEY1_quantile_95,HEY1_quantile]
    row_counter=row_counter+1
del HEY1
plt.figure(figsize=(20, 17))


NFIA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFIA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFIA.columns=['TF_TG_SCORE','label']
NFIA_sum = sum(NFIA['TF_TG_SCORE'])
NFIA_length = len(NFIA)
NFIA_mean=0
NFIA_quantile=0
NFIA_quantile_95=0
NFIA_median=0
NFIA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFIA['TF_TG_SCORE'])
if(NFIA_length>0):
    NFIA_mean=NFIA_sum/NFIA_length
    NFIA_quantile=np.percentile(NFIA['TF_TG_SCORE'], 99)
    NFIA_quantile_95=np.percentile(NFIA['TF_TG_SCORE'], 95)
    NFIA_median=sts.median(NFIA['TF_TG_SCORE'])
if(NFIA_median > background_median and NFIA_mannwhitneyU['pvalue']<0.01):
    background_NFIA = pd.concat([background,NFIA],axis=0)
    ax_NFIA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFIA,palette="Set3")
    ax_NFIA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFIA.png')
    del background_NFIA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFIA',NFIA_sum,NFIA_length,NFIA_mean,NFIA_median,NFIA_quantile_95,NFIA_quantile]
    row_counter=row_counter+1
del NFIA
plt.figure(figsize=(20, 17))


JUND=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/JUND_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
JUND.columns=['TF_TG_SCORE','label']
JUND_sum = sum(JUND['TF_TG_SCORE'])
JUND_length = len(JUND)
JUND_mean=0
JUND_quantile=0
JUND_quantile_95=0
JUND_median=0
JUND_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],JUND['TF_TG_SCORE'])
if(JUND_length>0):
    JUND_mean=JUND_sum/JUND_length
    JUND_quantile=np.percentile(JUND['TF_TG_SCORE'], 99)
    JUND_quantile_95=np.percentile(JUND['TF_TG_SCORE'], 95)
    JUND_median=sts.median(JUND['TF_TG_SCORE'])
if(JUND_median > background_median and JUND_mannwhitneyU['pvalue']<0.01):
    background_JUND = pd.concat([background,JUND],axis=0)
    ax_JUND = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_JUND,palette="Set3")
    ax_JUND.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/JUND.png')
    del background_JUND
    plt.clf()
    df_interesting_stats.loc[row_counter]=['JUND',JUND_sum,JUND_length,JUND_mean,JUND_median,JUND_quantile_95,JUND_quantile]
    row_counter=row_counter+1
del JUND
plt.figure(figsize=(20, 17))


ATF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ATF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ATF1.columns=['TF_TG_SCORE','label']
ATF1_sum = sum(ATF1['TF_TG_SCORE'])
ATF1_length = len(ATF1)
ATF1_mean=0
ATF1_quantile=0
ATF1_quantile_95=0
ATF1_median=0
ATF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ATF1['TF_TG_SCORE'])
if(ATF1_length>0):
    ATF1_mean=ATF1_sum/ATF1_length
    ATF1_quantile=np.percentile(ATF1['TF_TG_SCORE'], 99)
    ATF1_quantile_95=np.percentile(ATF1['TF_TG_SCORE'], 95)
    ATF1_median=sts.median(ATF1['TF_TG_SCORE'])
if(ATF1_median > background_median and ATF1_mannwhitneyU['pvalue']<0.01):
    background_ATF1 = pd.concat([background,ATF1],axis=0)
    ax_ATF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ATF1,palette="Set3")
    ax_ATF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ATF1.png')
    del background_ATF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ATF1',ATF1_sum,ATF1_length,ATF1_mean,ATF1_median,ATF1_quantile_95,ATF1_quantile]
    row_counter=row_counter+1
del ATF1
plt.figure(figsize=(20, 17))


RFX5=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RFX5_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RFX5.columns=['TF_TG_SCORE','label']
RFX5_sum = sum(RFX5['TF_TG_SCORE'])
RFX5_length = len(RFX5)
RFX5_mean=0
RFX5_quantile=0
RFX5_quantile_95=0
RFX5_median=0
RFX5_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RFX5['TF_TG_SCORE'])
if(RFX5_length>0):
    RFX5_mean=RFX5_sum/RFX5_length
    RFX5_quantile=np.percentile(RFX5['TF_TG_SCORE'], 99)
    RFX5_quantile_95=np.percentile(RFX5['TF_TG_SCORE'], 95)
    RFX5_median=sts.median(RFX5['TF_TG_SCORE'])
if(RFX5_median > background_median and RFX5_mannwhitneyU['pvalue']<0.01):
    background_RFX5 = pd.concat([background,RFX5],axis=0)
    ax_RFX5 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RFX5,palette="Set3")
    ax_RFX5.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RFX5.png')
    del background_RFX5
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RFX5',RFX5_sum,RFX5_length,RFX5_mean,RFX5_median,RFX5_quantile_95,RFX5_quantile]
    row_counter=row_counter+1
del RFX5
plt.figure(figsize=(20, 17))


ELF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ELF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ELF1.columns=['TF_TG_SCORE','label']
ELF1_sum = sum(ELF1['TF_TG_SCORE'])
ELF1_length = len(ELF1)
ELF1_mean=0
ELF1_quantile=0
ELF1_quantile_95=0
ELF1_median=0
ELF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ELF1['TF_TG_SCORE'])
if(ELF1_length>0):
    ELF1_mean=ELF1_sum/ELF1_length
    ELF1_quantile=np.percentile(ELF1['TF_TG_SCORE'], 99)
    ELF1_quantile_95=np.percentile(ELF1['TF_TG_SCORE'], 95)
    ELF1_median=sts.median(ELF1['TF_TG_SCORE'])
if(ELF1_median > background_median and ELF1_mannwhitneyU['pvalue']<0.01):
    background_ELF1 = pd.concat([background,ELF1],axis=0)
    ax_ELF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ELF1,palette="Set3")
    ax_ELF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ELF1.png')
    del background_ELF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ELF1',ELF1_sum,ELF1_length,ELF1_mean,ELF1_median,ELF1_quantile_95,ELF1_quantile]
    row_counter=row_counter+1
del ELF1
plt.figure(figsize=(20, 17))


HLF=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HLF_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HLF.columns=['TF_TG_SCORE','label']
HLF_sum = sum(HLF['TF_TG_SCORE'])
HLF_length = len(HLF)
HLF_mean=0
HLF_quantile=0
HLF_quantile_95=0
HLF_median=0
HLF_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HLF['TF_TG_SCORE'])
if(HLF_length>0):
    HLF_mean=HLF_sum/HLF_length
    HLF_quantile=np.percentile(HLF['TF_TG_SCORE'], 99)
    HLF_quantile_95=np.percentile(HLF['TF_TG_SCORE'], 95)
    HLF_median=sts.median(HLF['TF_TG_SCORE'])
if(HLF_median > background_median and HLF_mannwhitneyU['pvalue']<0.01):
    background_HLF = pd.concat([background,HLF],axis=0)
    ax_HLF = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HLF,palette="Set3")
    ax_HLF.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HLF.png')
    del background_HLF
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HLF',HLF_sum,HLF_length,HLF_mean,HLF_median,HLF_quantile_95,HLF_quantile]
    row_counter=row_counter+1
del HLF
plt.figure(figsize=(20, 17))


STAT3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/STAT3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
STAT3.columns=['TF_TG_SCORE','label']
STAT3_sum = sum(STAT3['TF_TG_SCORE'])
STAT3_length = len(STAT3)
STAT3_mean=0
STAT3_quantile=0
STAT3_quantile_95=0
STAT3_median=0
STAT3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],STAT3['TF_TG_SCORE'])
if(STAT3_length>0):
    STAT3_mean=STAT3_sum/STAT3_length
    STAT3_quantile=np.percentile(STAT3['TF_TG_SCORE'], 99)
    STAT3_quantile_95=np.percentile(STAT3['TF_TG_SCORE'], 95)
    STAT3_median=sts.median(STAT3['TF_TG_SCORE'])
if(STAT3_median > background_median and STAT3_mannwhitneyU['pvalue']<0.01):
    background_STAT3 = pd.concat([background,STAT3],axis=0)
    ax_STAT3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_STAT3,palette="Set3")
    ax_STAT3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/STAT3.png')
    del background_STAT3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['STAT3',STAT3_sum,STAT3_length,STAT3_mean,STAT3_median,STAT3_quantile_95,STAT3_quantile]
    row_counter=row_counter+1
del STAT3
plt.figure(figsize=(20, 17))


ELF5=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ELF5_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ELF5.columns=['TF_TG_SCORE','label']
ELF5_sum = sum(ELF5['TF_TG_SCORE'])
ELF5_length = len(ELF5)
ELF5_mean=0
ELF5_quantile=0
ELF5_quantile_95=0
ELF5_median=0
ELF5_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ELF5['TF_TG_SCORE'])
if(ELF5_length>0):
    ELF5_mean=ELF5_sum/ELF5_length
    ELF5_quantile=np.percentile(ELF5['TF_TG_SCORE'], 99)
    ELF5_quantile_95=np.percentile(ELF5['TF_TG_SCORE'], 95)
    ELF5_median=sts.median(ELF5['TF_TG_SCORE'])
if(ELF5_median > background_median and ELF5_mannwhitneyU['pvalue']<0.01):
    background_ELF5 = pd.concat([background,ELF5],axis=0)
    ax_ELF5 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ELF5,palette="Set3")
    ax_ELF5.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ELF5.png')
    del background_ELF5
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ELF5',ELF5_sum,ELF5_length,ELF5_mean,ELF5_median,ELF5_quantile_95,ELF5_quantile]
    row_counter=row_counter+1
del ELF5
plt.figure(figsize=(20, 17))


HSF2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/HSF2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
HSF2.columns=['TF_TG_SCORE','label']
HSF2_sum = sum(HSF2['TF_TG_SCORE'])
HSF2_length = len(HSF2)
HSF2_mean=0
HSF2_quantile=0
HSF2_quantile_95=0
HSF2_median=0
HSF2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],HSF2['TF_TG_SCORE'])
if(HSF2_length>0):
    HSF2_mean=HSF2_sum/HSF2_length
    HSF2_quantile=np.percentile(HSF2['TF_TG_SCORE'], 99)
    HSF2_quantile_95=np.percentile(HSF2['TF_TG_SCORE'], 95)
    HSF2_median=sts.median(HSF2['TF_TG_SCORE'])
if(HSF2_median > background_median and HSF2_mannwhitneyU['pvalue']<0.01):
    background_HSF2 = pd.concat([background,HSF2],axis=0)
    ax_HSF2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_HSF2,palette="Set3")
    ax_HSF2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/HSF2.png')
    del background_HSF2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['HSF2',HSF2_sum,HSF2_length,HSF2_mean,HSF2_median,HSF2_quantile_95,HSF2_quantile]
    row_counter=row_counter+1
del HSF2
plt.figure(figsize=(20, 17))


ATF4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ATF4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ATF4.columns=['TF_TG_SCORE','label']
ATF4_sum = sum(ATF4['TF_TG_SCORE'])
ATF4_length = len(ATF4)
ATF4_mean=0
ATF4_quantile=0
ATF4_quantile_95=0
ATF4_median=0
ATF4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ATF4['TF_TG_SCORE'])
if(ATF4_length>0):
    ATF4_mean=ATF4_sum/ATF4_length
    ATF4_quantile=np.percentile(ATF4['TF_TG_SCORE'], 99)
    ATF4_quantile_95=np.percentile(ATF4['TF_TG_SCORE'], 95)
    ATF4_median=sts.median(ATF4['TF_TG_SCORE'])
if(ATF4_median > background_median and ATF4_mannwhitneyU['pvalue']<0.01):
    background_ATF4 = pd.concat([background,ATF4],axis=0)
    ax_ATF4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ATF4,palette="Set3")
    ax_ATF4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ATF4.png')
    del background_ATF4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ATF4',ATF4_sum,ATF4_length,ATF4_mean,ATF4_median,ATF4_quantile_95,ATF4_quantile]
    row_counter=row_counter+1
del ATF4
plt.figure(figsize=(20, 17))


ETV6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ETV6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ETV6.columns=['TF_TG_SCORE','label']
ETV6_sum = sum(ETV6['TF_TG_SCORE'])
ETV6_length = len(ETV6)
ETV6_mean=0
ETV6_quantile=0
ETV6_quantile_95=0
ETV6_median=0
ETV6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ETV6['TF_TG_SCORE'])
if(ETV6_length>0):
    ETV6_mean=ETV6_sum/ETV6_length
    ETV6_quantile=np.percentile(ETV6['TF_TG_SCORE'], 99)
    ETV6_quantile_95=np.percentile(ETV6['TF_TG_SCORE'], 95)
    ETV6_median=sts.median(ETV6['TF_TG_SCORE'])
if(ETV6_median > background_median and ETV6_mannwhitneyU['pvalue']<0.01):
    background_ETV6 = pd.concat([background,ETV6],axis=0)
    ax_ETV6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ETV6,palette="Set3")
    ax_ETV6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ETV6.png')
    del background_ETV6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ETV6',ETV6_sum,ETV6_length,ETV6_mean,ETV6_median,ETV6_quantile_95,ETV6_quantile]
    row_counter=row_counter+1
del ETV6
plt.figure(figsize=(20, 17))


EBF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/EBF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
EBF1.columns=['TF_TG_SCORE','label']
EBF1_sum = sum(EBF1['TF_TG_SCORE'])
EBF1_length = len(EBF1)
EBF1_mean=0
EBF1_quantile=0
EBF1_quantile_95=0
EBF1_median=0
EBF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],EBF1['TF_TG_SCORE'])
if(EBF1_length>0):
    EBF1_mean=EBF1_sum/EBF1_length
    EBF1_quantile=np.percentile(EBF1['TF_TG_SCORE'], 99)
    EBF1_quantile_95=np.percentile(EBF1['TF_TG_SCORE'], 95)
    EBF1_median=sts.median(EBF1['TF_TG_SCORE'])
if(EBF1_median > background_median and EBF1_mannwhitneyU['pvalue']<0.01):
    background_EBF1 = pd.concat([background,EBF1],axis=0)
    ax_EBF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_EBF1,palette="Set3")
    ax_EBF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/EBF1.png')
    del background_EBF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['EBF1',EBF1_sum,EBF1_length,EBF1_mean,EBF1_median,EBF1_quantile_95,EBF1_quantile]
    row_counter=row_counter+1
del EBF1
plt.figure(figsize=(20, 17))


ARID3A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ARID3A_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ARID3A.columns=['TF_TG_SCORE','label']
ARID3A_sum = sum(ARID3A['TF_TG_SCORE'])
ARID3A_length = len(ARID3A)
ARID3A_mean=0
ARID3A_quantile=0
ARID3A_quantile_95=0
ARID3A_median=0
ARID3A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ARID3A['TF_TG_SCORE'])
if(ARID3A_length>0):
    ARID3A_mean=ARID3A_sum/ARID3A_length
    ARID3A_quantile=np.percentile(ARID3A['TF_TG_SCORE'], 99)
    ARID3A_quantile_95=np.percentile(ARID3A['TF_TG_SCORE'], 95)
    ARID3A_median=sts.median(ARID3A['TF_TG_SCORE'])
if(ARID3A_median > background_median and ARID3A_mannwhitneyU['pvalue']<0.01):
    background_ARID3A = pd.concat([background,ARID3A],axis=0)
    ax_ARID3A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ARID3A,palette="Set3")
    ax_ARID3A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ARID3A.png')
    del background_ARID3A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ARID3A',ARID3A_sum,ARID3A_length,ARID3A_mean,ARID3A_median,ARID3A_quantile_95,ARID3A_quantile]
    row_counter=row_counter+1
del ARID3A
plt.figure(figsize=(20, 17))


E2F7=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/E2F7_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
E2F7.columns=['TF_TG_SCORE','label']
E2F7_sum = sum(E2F7['TF_TG_SCORE'])
E2F7_length = len(E2F7)
E2F7_mean=0
E2F7_quantile=0
E2F7_quantile_95=0
E2F7_median=0
E2F7_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],E2F7['TF_TG_SCORE'])
if(E2F7_length>0):
    E2F7_mean=E2F7_sum/E2F7_length
    E2F7_quantile=np.percentile(E2F7['TF_TG_SCORE'], 99)
    E2F7_quantile_95=np.percentile(E2F7['TF_TG_SCORE'], 95)
    E2F7_median=sts.median(E2F7['TF_TG_SCORE'])
if(E2F7_median > background_median and E2F7_mannwhitneyU['pvalue']<0.01):
    background_E2F7 = pd.concat([background,E2F7],axis=0)
    ax_E2F7 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_E2F7,palette="Set3")
    ax_E2F7.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/E2F7.png')
    del background_E2F7
    plt.clf()
    df_interesting_stats.loc[row_counter]=['E2F7',E2F7_sum,E2F7_length,E2F7_mean,E2F7_median,E2F7_quantile_95,E2F7_quantile]
    row_counter=row_counter+1
del E2F7
plt.figure(figsize=(20, 17))


KLF15=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/KLF15_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
KLF15.columns=['TF_TG_SCORE','label']
KLF15_sum = sum(KLF15['TF_TG_SCORE'])
KLF15_length = len(KLF15)
KLF15_mean=0
KLF15_quantile=0
KLF15_quantile_95=0
KLF15_median=0
KLF15_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],KLF15['TF_TG_SCORE'])
if(KLF15_length>0):
    KLF15_mean=KLF15_sum/KLF15_length
    KLF15_quantile=np.percentile(KLF15['TF_TG_SCORE'], 99)
    KLF15_quantile_95=np.percentile(KLF15['TF_TG_SCORE'], 95)
    KLF15_median=sts.median(KLF15['TF_TG_SCORE'])
if(KLF15_median > background_median and KLF15_mannwhitneyU['pvalue']<0.01):
    background_KLF15 = pd.concat([background,KLF15],axis=0)
    ax_KLF15 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_KLF15,palette="Set3")
    ax_KLF15.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/KLF15.png')
    del background_KLF15
    plt.clf()
    df_interesting_stats.loc[row_counter]=['KLF15',KLF15_sum,KLF15_length,KLF15_mean,KLF15_median,KLF15_quantile_95,KLF15_quantile]
    row_counter=row_counter+1
del KLF15
plt.figure(figsize=(20, 17))


CREB3L2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CREB3L2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CREB3L2.columns=['TF_TG_SCORE','label']
CREB3L2_sum = sum(CREB3L2['TF_TG_SCORE'])
CREB3L2_length = len(CREB3L2)
CREB3L2_mean=0
CREB3L2_quantile=0
CREB3L2_quantile_95=0
CREB3L2_median=0
CREB3L2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CREB3L2['TF_TG_SCORE'])
if(CREB3L2_length>0):
    CREB3L2_mean=CREB3L2_sum/CREB3L2_length
    CREB3L2_quantile=np.percentile(CREB3L2['TF_TG_SCORE'], 99)
    CREB3L2_quantile_95=np.percentile(CREB3L2['TF_TG_SCORE'], 95)
    CREB3L2_median=sts.median(CREB3L2['TF_TG_SCORE'])
if(CREB3L2_median > background_median and CREB3L2_mannwhitneyU['pvalue']<0.01):
    background_CREB3L2 = pd.concat([background,CREB3L2],axis=0)
    ax_CREB3L2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CREB3L2,palette="Set3")
    ax_CREB3L2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CREB3L2.png')
    del background_CREB3L2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CREB3L2',CREB3L2_sum,CREB3L2_length,CREB3L2_mean,CREB3L2_median,CREB3L2_quantile_95,CREB3L2_quantile]
    row_counter=row_counter+1
del CREB3L2
plt.figure(figsize=(20, 17))


SIN3A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SIN3A_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SIN3A.columns=['TF_TG_SCORE','label']
SIN3A_sum = sum(SIN3A['TF_TG_SCORE'])
SIN3A_length = len(SIN3A)
SIN3A_mean=0
SIN3A_quantile=0
SIN3A_quantile_95=0
SIN3A_median=0
SIN3A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SIN3A['TF_TG_SCORE'])
if(SIN3A_length>0):
    SIN3A_mean=SIN3A_sum/SIN3A_length
    SIN3A_quantile=np.percentile(SIN3A['TF_TG_SCORE'], 99)
    SIN3A_quantile_95=np.percentile(SIN3A['TF_TG_SCORE'], 95)
    SIN3A_median=sts.median(SIN3A['TF_TG_SCORE'])
if(SIN3A_median > background_median and SIN3A_mannwhitneyU['pvalue']<0.01):
    background_SIN3A = pd.concat([background,SIN3A],axis=0)
    ax_SIN3A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SIN3A,palette="Set3")
    ax_SIN3A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SIN3A.png')
    del background_SIN3A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SIN3A',SIN3A_sum,SIN3A_length,SIN3A_mean,SIN3A_median,SIN3A_quantile_95,SIN3A_quantile]
    row_counter=row_counter+1
del SIN3A
plt.figure(figsize=(20, 17))


RUNX1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RUNX1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RUNX1.columns=['TF_TG_SCORE','label']
RUNX1_sum = sum(RUNX1['TF_TG_SCORE'])
RUNX1_length = len(RUNX1)
RUNX1_mean=0
RUNX1_quantile=0
RUNX1_quantile_95=0
RUNX1_median=0
RUNX1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RUNX1['TF_TG_SCORE'])
if(RUNX1_length>0):
    RUNX1_mean=RUNX1_sum/RUNX1_length
    RUNX1_quantile=np.percentile(RUNX1['TF_TG_SCORE'], 99)
    RUNX1_quantile_95=np.percentile(RUNX1['TF_TG_SCORE'], 95)
    RUNX1_median=sts.median(RUNX1['TF_TG_SCORE'])
if(RUNX1_median > background_median and RUNX1_mannwhitneyU['pvalue']<0.01):
    background_RUNX1 = pd.concat([background,RUNX1],axis=0)
    ax_RUNX1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RUNX1,palette="Set3")
    ax_RUNX1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RUNX1.png')
    del background_RUNX1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RUNX1',RUNX1_sum,RUNX1_length,RUNX1_mean,RUNX1_median,RUNX1_quantile_95,RUNX1_quantile]
    row_counter=row_counter+1
del RUNX1
plt.figure(figsize=(20, 17))


DDIT3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/DDIT3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
DDIT3.columns=['TF_TG_SCORE','label']
DDIT3_sum = sum(DDIT3['TF_TG_SCORE'])
DDIT3_length = len(DDIT3)
DDIT3_mean=0
DDIT3_quantile=0
DDIT3_quantile_95=0
DDIT3_median=0
DDIT3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],DDIT3['TF_TG_SCORE'])
if(DDIT3_length>0):
    DDIT3_mean=DDIT3_sum/DDIT3_length
    DDIT3_quantile=np.percentile(DDIT3['TF_TG_SCORE'], 99)
    DDIT3_quantile_95=np.percentile(DDIT3['TF_TG_SCORE'], 95)
    DDIT3_median=sts.median(DDIT3['TF_TG_SCORE'])
if(DDIT3_median > background_median and DDIT3_mannwhitneyU['pvalue']<0.01):
    background_DDIT3 = pd.concat([background,DDIT3],axis=0)
    ax_DDIT3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_DDIT3,palette="Set3")
    ax_DDIT3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/DDIT3.png')
    del background_DDIT3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['DDIT3',DDIT3_sum,DDIT3_length,DDIT3_mean,DDIT3_median,DDIT3_quantile_95,DDIT3_quantile]
    row_counter=row_counter+1
del DDIT3
plt.figure(figsize=(20, 17))


NFKB2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFKB2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFKB2.columns=['TF_TG_SCORE','label']
NFKB2_sum = sum(NFKB2['TF_TG_SCORE'])
NFKB2_length = len(NFKB2)
NFKB2_mean=0
NFKB2_quantile=0
NFKB2_quantile_95=0
NFKB2_median=0
NFKB2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFKB2['TF_TG_SCORE'])
if(NFKB2_length>0):
    NFKB2_mean=NFKB2_sum/NFKB2_length
    NFKB2_quantile=np.percentile(NFKB2['TF_TG_SCORE'], 99)
    NFKB2_quantile_95=np.percentile(NFKB2['TF_TG_SCORE'], 95)
    NFKB2_median=sts.median(NFKB2['TF_TG_SCORE'])
if(NFKB2_median > background_median and NFKB2_mannwhitneyU['pvalue']<0.01):
    background_NFKB2 = pd.concat([background,NFKB2],axis=0)
    ax_NFKB2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFKB2,palette="Set3")
    ax_NFKB2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFKB2.png')
    del background_NFKB2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFKB2',NFKB2_sum,NFKB2_length,NFKB2_mean,NFKB2_median,NFKB2_quantile_95,NFKB2_quantile]
    row_counter=row_counter+1
del NFKB2
plt.figure(figsize=(20, 17))


TCF7L2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TCF7L2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TCF7L2.columns=['TF_TG_SCORE','label']
TCF7L2_sum = sum(TCF7L2['TF_TG_SCORE'])
TCF7L2_length = len(TCF7L2)
TCF7L2_mean=0
TCF7L2_quantile=0
TCF7L2_quantile_95=0
TCF7L2_median=0
TCF7L2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TCF7L2['TF_TG_SCORE'])
if(TCF7L2_length>0):
    TCF7L2_mean=TCF7L2_sum/TCF7L2_length
    TCF7L2_quantile=np.percentile(TCF7L2['TF_TG_SCORE'], 99)
    TCF7L2_quantile_95=np.percentile(TCF7L2['TF_TG_SCORE'], 95)
    TCF7L2_median=sts.median(TCF7L2['TF_TG_SCORE'])
if(TCF7L2_median > background_median and TCF7L2_mannwhitneyU['pvalue']<0.01):
    background_TCF7L2 = pd.concat([background,TCF7L2],axis=0)
    ax_TCF7L2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TCF7L2,palette="Set3")
    ax_TCF7L2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TCF7L2.png')
    del background_TCF7L2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TCF7L2',TCF7L2_sum,TCF7L2_length,TCF7L2_mean,TCF7L2_median,TCF7L2_quantile_95,TCF7L2_quantile]
    row_counter=row_counter+1
del TCF7L2
plt.figure(figsize=(20, 17))


SETDB1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SETDB1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SETDB1.columns=['TF_TG_SCORE','label']
SETDB1_sum = sum(SETDB1['TF_TG_SCORE'])
SETDB1_length = len(SETDB1)
SETDB1_mean=0
SETDB1_quantile=0
SETDB1_quantile_95=0
SETDB1_median=0
SETDB1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SETDB1['TF_TG_SCORE'])
if(SETDB1_length>0):
    SETDB1_mean=SETDB1_sum/SETDB1_length
    SETDB1_quantile=np.percentile(SETDB1['TF_TG_SCORE'], 99)
    SETDB1_quantile_95=np.percentile(SETDB1['TF_TG_SCORE'], 95)
    SETDB1_median=sts.median(SETDB1['TF_TG_SCORE'])
if(SETDB1_median > background_median and SETDB1_mannwhitneyU['pvalue']<0.01):
    background_SETDB1 = pd.concat([background,SETDB1],axis=0)
    ax_SETDB1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SETDB1,palette="Set3")
    ax_SETDB1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SETDB1.png')
    del background_SETDB1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SETDB1',SETDB1_sum,SETDB1_length,SETDB1_mean,SETDB1_median,SETDB1_quantile_95,SETDB1_quantile]
    row_counter=row_counter+1
del SETDB1
plt.figure(figsize=(20, 17))


PBX1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PBX1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PBX1.columns=['TF_TG_SCORE','label']
PBX1_sum = sum(PBX1['TF_TG_SCORE'])
PBX1_length = len(PBX1)
PBX1_mean=0
PBX1_quantile=0
PBX1_quantile_95=0
PBX1_median=0
PBX1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PBX1['TF_TG_SCORE'])
if(PBX1_length>0):
    PBX1_mean=PBX1_sum/PBX1_length
    PBX1_quantile=np.percentile(PBX1['TF_TG_SCORE'], 99)
    PBX1_quantile_95=np.percentile(PBX1['TF_TG_SCORE'], 95)
    PBX1_median=sts.median(PBX1['TF_TG_SCORE'])
if(PBX1_median > background_median and PBX1_mannwhitneyU['pvalue']<0.01):
    background_PBX1 = pd.concat([background,PBX1],axis=0)
    ax_PBX1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PBX1,palette="Set3")
    ax_PBX1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PBX1.png')
    del background_PBX1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PBX1',PBX1_sum,PBX1_length,PBX1_mean,PBX1_median,PBX1_quantile_95,PBX1_quantile]
    row_counter=row_counter+1
del PBX1
plt.figure(figsize=(20, 17))


YY1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/YY1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
YY1.columns=['TF_TG_SCORE','label']
YY1_sum = sum(YY1['TF_TG_SCORE'])
YY1_length = len(YY1)
YY1_mean=0
YY1_quantile=0
YY1_quantile_95=0
YY1_median=0
YY1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],YY1['TF_TG_SCORE'])
if(YY1_length>0):
    YY1_mean=YY1_sum/YY1_length
    YY1_quantile=np.percentile(YY1['TF_TG_SCORE'], 99)
    YY1_quantile_95=np.percentile(YY1['TF_TG_SCORE'], 95)
    YY1_median=sts.median(YY1['TF_TG_SCORE'])
if(YY1_median > background_median and YY1_mannwhitneyU['pvalue']<0.01):
    background_YY1 = pd.concat([background,YY1],axis=0)
    ax_YY1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_YY1,palette="Set3")
    ax_YY1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/YY1.png')
    del background_YY1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['YY1',YY1_sum,YY1_length,YY1_mean,YY1_median,YY1_quantile_95,YY1_quantile]
    row_counter=row_counter+1
del YY1
plt.figure(figsize=(20, 17))


ARID5A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ARID5A_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ARID5A.columns=['TF_TG_SCORE','label']
ARID5A_sum = sum(ARID5A['TF_TG_SCORE'])
ARID5A_length = len(ARID5A)
ARID5A_mean=0
ARID5A_quantile=0
ARID5A_quantile_95=0
ARID5A_median=0
ARID5A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ARID5A['TF_TG_SCORE'])
if(ARID5A_length>0):
    ARID5A_mean=ARID5A_sum/ARID5A_length
    ARID5A_quantile=np.percentile(ARID5A['TF_TG_SCORE'], 99)
    ARID5A_quantile_95=np.percentile(ARID5A['TF_TG_SCORE'], 95)
    ARID5A_median=sts.median(ARID5A['TF_TG_SCORE'])
if(ARID5A_median > background_median and ARID5A_mannwhitneyU['pvalue']<0.01):
    background_ARID5A = pd.concat([background,ARID5A],axis=0)
    ax_ARID5A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ARID5A,palette="Set3")
    ax_ARID5A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ARID5A.png')
    del background_ARID5A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ARID5A',ARID5A_sum,ARID5A_length,ARID5A_mean,ARID5A_median,ARID5A_quantile_95,ARID5A_quantile]
    row_counter=row_counter+1
del ARID5A
plt.figure(figsize=(20, 17))


PPARG=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PPARG..RXRA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PPARG.columns=['TF_TG_SCORE','label']
PPARG_sum = sum(PPARG['TF_TG_SCORE'])
PPARG_length = len(PPARG)
PPARG_mean=0
PPARG_quantile=0
PPARG_quantile_95=0
PPARG_median=0
PPARG_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PPARG['TF_TG_SCORE'])
if(PPARG_length>0):
    PPARG_mean=PPARG_sum/PPARG_length
    PPARG_quantile=np.percentile(PPARG['TF_TG_SCORE'], 99)
    PPARG_quantile_95=np.percentile(PPARG['TF_TG_SCORE'], 95)
    PPARG_median=sts.median(PPARG['TF_TG_SCORE'])
if(PPARG_median > background_median and PPARG_mannwhitneyU['pvalue']<0.01):
    background_PPARG = pd.concat([background,PPARG],axis=0)
    ax_PPARG = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PPARG,palette="Set3")
    ax_PPARG.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PPARG..RXRA.png')
    del background_PPARG
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PPARG..RXRA',PPARG_sum,PPARG_length,PPARG_mean,PPARG_median,PPARG_quantile_95,PPARG_quantile]
    row_counter=row_counter+1
del PPARG
plt.figure(figsize=(20, 17))


E2F4=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/E2F4_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
E2F4.columns=['TF_TG_SCORE','label']
E2F4_sum = sum(E2F4['TF_TG_SCORE'])
E2F4_length = len(E2F4)
E2F4_mean=0
E2F4_quantile=0
E2F4_quantile_95=0
E2F4_median=0
E2F4_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],E2F4['TF_TG_SCORE'])
if(E2F4_length>0):
    E2F4_mean=E2F4_sum/E2F4_length
    E2F4_quantile=np.percentile(E2F4['TF_TG_SCORE'], 99)
    E2F4_quantile_95=np.percentile(E2F4['TF_TG_SCORE'], 95)
    E2F4_median=sts.median(E2F4['TF_TG_SCORE'])
if(E2F4_median > background_median and E2F4_mannwhitneyU['pvalue']<0.01):
    background_E2F4 = pd.concat([background,E2F4],axis=0)
    ax_E2F4 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_E2F4,palette="Set3")
    ax_E2F4.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/E2F4.png')
    del background_E2F4
    plt.clf()
    df_interesting_stats.loc[row_counter]=['E2F4',E2F4_sum,E2F4_length,E2F4_mean,E2F4_median,E2F4_quantile_95,E2F4_quantile]
    row_counter=row_counter+1
del E2F4
plt.figure(figsize=(20, 17))


RARA=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/RARA_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
RARA.columns=['TF_TG_SCORE','label']
RARA_sum = sum(RARA['TF_TG_SCORE'])
RARA_length = len(RARA)
RARA_mean=0
RARA_quantile=0
RARA_quantile_95=0
RARA_median=0
RARA_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],RARA['TF_TG_SCORE'])
if(RARA_length>0):
    RARA_mean=RARA_sum/RARA_length
    RARA_quantile=np.percentile(RARA['TF_TG_SCORE'], 99)
    RARA_quantile_95=np.percentile(RARA['TF_TG_SCORE'], 95)
    RARA_median=sts.median(RARA['TF_TG_SCORE'])
if(RARA_median > background_median and RARA_mannwhitneyU['pvalue']<0.01):
    background_RARA = pd.concat([background,RARA],axis=0)
    ax_RARA = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_RARA,palette="Set3")
    ax_RARA.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/RARA.png')
    del background_RARA
    plt.clf()
    df_interesting_stats.loc[row_counter]=['RARA',RARA_sum,RARA_length,RARA_mean,RARA_median,RARA_quantile_95,RARA_quantile]
    row_counter=row_counter+1
del RARA
plt.figure(figsize=(20, 17))


IRF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF1.columns=['TF_TG_SCORE','label']
IRF1_sum = sum(IRF1['TF_TG_SCORE'])
IRF1_length = len(IRF1)
IRF1_mean=0
IRF1_quantile=0
IRF1_quantile_95=0
IRF1_median=0
IRF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF1['TF_TG_SCORE'])
if(IRF1_length>0):
    IRF1_mean=IRF1_sum/IRF1_length
    IRF1_quantile=np.percentile(IRF1['TF_TG_SCORE'], 99)
    IRF1_quantile_95=np.percentile(IRF1['TF_TG_SCORE'], 95)
    IRF1_median=sts.median(IRF1['TF_TG_SCORE'])
if(IRF1_median > background_median and IRF1_mannwhitneyU['pvalue']<0.01):
    background_IRF1 = pd.concat([background,IRF1],axis=0)
    ax_IRF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF1,palette="Set3")
    ax_IRF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF1.png')
    del background_IRF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF1',IRF1_sum,IRF1_length,IRF1_mean,IRF1_median,IRF1_quantile_95,IRF1_quantile]
    row_counter=row_counter+1
del IRF1
plt.figure(figsize=(20, 17))


CUX1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CUX1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CUX1.columns=['TF_TG_SCORE','label']
CUX1_sum = sum(CUX1['TF_TG_SCORE'])
CUX1_length = len(CUX1)
CUX1_mean=0
CUX1_quantile=0
CUX1_quantile_95=0
CUX1_median=0
CUX1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CUX1['TF_TG_SCORE'])
if(CUX1_length>0):
    CUX1_mean=CUX1_sum/CUX1_length
    CUX1_quantile=np.percentile(CUX1['TF_TG_SCORE'], 99)
    CUX1_quantile_95=np.percentile(CUX1['TF_TG_SCORE'], 95)
    CUX1_median=sts.median(CUX1['TF_TG_SCORE'])
if(CUX1_median > background_median and CUX1_mannwhitneyU['pvalue']<0.01):
    background_CUX1 = pd.concat([background,CUX1],axis=0)
    ax_CUX1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CUX1,palette="Set3")
    ax_CUX1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CUX1.png')
    del background_CUX1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CUX1',CUX1_sum,CUX1_length,CUX1_mean,CUX1_median,CUX1_quantile_95,CUX1_quantile]
    row_counter=row_counter+1
del CUX1
plt.figure(figsize=(20, 17))


IRF2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/IRF2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
IRF2.columns=['TF_TG_SCORE','label']
IRF2_sum = sum(IRF2['TF_TG_SCORE'])
IRF2_length = len(IRF2)
IRF2_mean=0
IRF2_quantile=0
IRF2_quantile_95=0
IRF2_median=0
IRF2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],IRF2['TF_TG_SCORE'])
if(IRF2_length>0):
    IRF2_mean=IRF2_sum/IRF2_length
    IRF2_quantile=np.percentile(IRF2['TF_TG_SCORE'], 99)
    IRF2_quantile_95=np.percentile(IRF2['TF_TG_SCORE'], 95)
    IRF2_median=sts.median(IRF2['TF_TG_SCORE'])
if(IRF2_median > background_median and IRF2_mannwhitneyU['pvalue']<0.01):
    background_IRF2 = pd.concat([background,IRF2],axis=0)
    ax_IRF2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_IRF2,palette="Set3")
    ax_IRF2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/IRF2.png')
    del background_IRF2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['IRF2',IRF2_sum,IRF2_length,IRF2_mean,IRF2_median,IRF2_quantile_95,IRF2_quantile]
    row_counter=row_counter+1
del IRF2
plt.figure(figsize=(20, 17))


ETS2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ETS2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ETS2.columns=['TF_TG_SCORE','label']
ETS2_sum = sum(ETS2['TF_TG_SCORE'])
ETS2_length = len(ETS2)
ETS2_mean=0
ETS2_quantile=0
ETS2_quantile_95=0
ETS2_median=0
ETS2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ETS2['TF_TG_SCORE'])
if(ETS2_length>0):
    ETS2_mean=ETS2_sum/ETS2_length
    ETS2_quantile=np.percentile(ETS2['TF_TG_SCORE'], 99)
    ETS2_quantile_95=np.percentile(ETS2['TF_TG_SCORE'], 95)
    ETS2_median=sts.median(ETS2['TF_TG_SCORE'])
if(ETS2_median > background_median and ETS2_mannwhitneyU['pvalue']<0.01):
    background_ETS2 = pd.concat([background,ETS2],axis=0)
    ax_ETS2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ETS2,palette="Set3")
    ax_ETS2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ETS2.png')
    del background_ETS2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ETS2',ETS2_sum,ETS2_length,ETS2_mean,ETS2_median,ETS2_quantile_95,ETS2_quantile]
    row_counter=row_counter+1
del ETS2
plt.figure(figsize=(20, 17))


TAF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TAF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TAF1.columns=['TF_TG_SCORE','label']
TAF1_sum = sum(TAF1['TF_TG_SCORE'])
TAF1_length = len(TAF1)
TAF1_mean=0
TAF1_quantile=0
TAF1_quantile_95=0
TAF1_median=0
TAF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TAF1['TF_TG_SCORE'])
if(TAF1_length>0):
    TAF1_mean=TAF1_sum/TAF1_length
    TAF1_quantile=np.percentile(TAF1['TF_TG_SCORE'], 99)
    TAF1_quantile_95=np.percentile(TAF1['TF_TG_SCORE'], 95)
    TAF1_median=sts.median(TAF1['TF_TG_SCORE'])
if(TAF1_median > background_median and TAF1_mannwhitneyU['pvalue']<0.01):
    background_TAF1 = pd.concat([background,TAF1],axis=0)
    ax_TAF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TAF1,palette="Set3")
    ax_TAF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TAF1.png')
    del background_TAF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TAF1',TAF1_sum,TAF1_length,TAF1_mean,TAF1_median,TAF1_quantile_95,TAF1_quantile]
    row_counter=row_counter+1
del TAF1
plt.figure(figsize=(20, 17))


XBP1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/XBP1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
XBP1.columns=['TF_TG_SCORE','label']
XBP1_sum = sum(XBP1['TF_TG_SCORE'])
XBP1_length = len(XBP1)
XBP1_mean=0
XBP1_quantile=0
XBP1_quantile_95=0
XBP1_median=0
XBP1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],XBP1['TF_TG_SCORE'])
if(XBP1_length>0):
    XBP1_mean=XBP1_sum/XBP1_length
    XBP1_quantile=np.percentile(XBP1['TF_TG_SCORE'], 99)
    XBP1_quantile_95=np.percentile(XBP1['TF_TG_SCORE'], 95)
    XBP1_median=sts.median(XBP1['TF_TG_SCORE'])
if(XBP1_median > background_median and XBP1_mannwhitneyU['pvalue']<0.01):
    background_XBP1 = pd.concat([background,XBP1],axis=0)
    ax_XBP1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_XBP1,palette="Set3")
    ax_XBP1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/XBP1.png')
    del background_XBP1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['XBP1',XBP1_sum,XBP1_length,XBP1_mean,XBP1_median,XBP1_quantile_95,XBP1_quantile]
    row_counter=row_counter+1
del XBP1
plt.figure(figsize=(20, 17))


STAT5A=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/STAT5A..STAT5B_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
STAT5A.columns=['TF_TG_SCORE','label']
STAT5A_sum = sum(STAT5A['TF_TG_SCORE'])
STAT5A_length = len(STAT5A)
STAT5A_mean=0
STAT5A_quantile=0
STAT5A_quantile_95=0
STAT5A_median=0
STAT5A_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],STAT5A['TF_TG_SCORE'])
if(STAT5A_length>0):
    STAT5A_mean=STAT5A_sum/STAT5A_length
    STAT5A_quantile=np.percentile(STAT5A['TF_TG_SCORE'], 99)
    STAT5A_quantile_95=np.percentile(STAT5A['TF_TG_SCORE'], 95)
    STAT5A_median=sts.median(STAT5A['TF_TG_SCORE'])
if(STAT5A_median > background_median and STAT5A_mannwhitneyU['pvalue']<0.01):
    background_STAT5A = pd.concat([background,STAT5A],axis=0)
    ax_STAT5A = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_STAT5A,palette="Set3")
    ax_STAT5A.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/STAT5A..STAT5B.png')
    del background_STAT5A
    plt.clf()
    df_interesting_stats.loc[row_counter]=['STAT5A..STAT5B',STAT5A_sum,STAT5A_length,STAT5A_mean,STAT5A_median,STAT5A_quantile_95,STAT5A_quantile]
    row_counter=row_counter+1
del STAT5A
plt.figure(figsize=(20, 17))


TFEB=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TFEB_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TFEB.columns=['TF_TG_SCORE','label']
TFEB_sum = sum(TFEB['TF_TG_SCORE'])
TFEB_length = len(TFEB)
TFEB_mean=0
TFEB_quantile=0
TFEB_quantile_95=0
TFEB_median=0
TFEB_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TFEB['TF_TG_SCORE'])
if(TFEB_length>0):
    TFEB_mean=TFEB_sum/TFEB_length
    TFEB_quantile=np.percentile(TFEB['TF_TG_SCORE'], 99)
    TFEB_quantile_95=np.percentile(TFEB['TF_TG_SCORE'], 95)
    TFEB_median=sts.median(TFEB['TF_TG_SCORE'])
if(TFEB_median > background_median and TFEB_mannwhitneyU['pvalue']<0.01):
    background_TFEB = pd.concat([background,TFEB],axis=0)
    ax_TFEB = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TFEB,palette="Set3")
    ax_TFEB.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TFEB.png')
    del background_TFEB
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TFEB',TFEB_sum,TFEB_length,TFEB_mean,TFEB_median,TFEB_quantile_95,TFEB_quantile]
    row_counter=row_counter+1
del TFEB
plt.figure(figsize=(20, 17))


TBX3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TBX3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TBX3.columns=['TF_TG_SCORE','label']
TBX3_sum = sum(TBX3['TF_TG_SCORE'])
TBX3_length = len(TBX3)
TBX3_mean=0
TBX3_quantile=0
TBX3_quantile_95=0
TBX3_median=0
TBX3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TBX3['TF_TG_SCORE'])
if(TBX3_length>0):
    TBX3_mean=TBX3_sum/TBX3_length
    TBX3_quantile=np.percentile(TBX3['TF_TG_SCORE'], 99)
    TBX3_quantile_95=np.percentile(TBX3['TF_TG_SCORE'], 95)
    TBX3_median=sts.median(TBX3['TF_TG_SCORE'])
if(TBX3_median > background_median and TBX3_mannwhitneyU['pvalue']<0.01):
    background_TBX3 = pd.concat([background,TBX3],axis=0)
    ax_TBX3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TBX3,palette="Set3")
    ax_TBX3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TBX3.png')
    del background_TBX3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TBX3',TBX3_sum,TBX3_length,TBX3_mean,TBX3_median,TBX3_quantile_95,TBX3_quantile]
    row_counter=row_counter+1
del TBX3
plt.figure(figsize=(20, 17))


GATA2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/GATA2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
GATA2.columns=['TF_TG_SCORE','label']
GATA2_sum = sum(GATA2['TF_TG_SCORE'])
GATA2_length = len(GATA2)
GATA2_mean=0
GATA2_quantile=0
GATA2_quantile_95=0
GATA2_median=0
GATA2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],GATA2['TF_TG_SCORE'])
if(GATA2_length>0):
    GATA2_mean=GATA2_sum/GATA2_length
    GATA2_quantile=np.percentile(GATA2['TF_TG_SCORE'], 99)
    GATA2_quantile_95=np.percentile(GATA2['TF_TG_SCORE'], 95)
    GATA2_median=sts.median(GATA2['TF_TG_SCORE'])
if(GATA2_median > background_median and GATA2_mannwhitneyU['pvalue']<0.01):
    background_GATA2 = pd.concat([background,GATA2],axis=0)
    ax_GATA2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_GATA2,palette="Set3")
    ax_GATA2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/GATA2.png')
    del background_GATA2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['GATA2',GATA2_sum,GATA2_length,GATA2_mean,GATA2_median,GATA2_quantile_95,GATA2_quantile]
    row_counter=row_counter+1
del GATA2
plt.figure(figsize=(20, 17))


STAT2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/STAT2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
STAT2.columns=['TF_TG_SCORE','label']
STAT2_sum = sum(STAT2['TF_TG_SCORE'])
STAT2_length = len(STAT2)
STAT2_mean=0
STAT2_quantile=0
STAT2_quantile_95=0
STAT2_median=0
STAT2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],STAT2['TF_TG_SCORE'])
if(STAT2_length>0):
    STAT2_mean=STAT2_sum/STAT2_length
    STAT2_quantile=np.percentile(STAT2['TF_TG_SCORE'], 99)
    STAT2_quantile_95=np.percentile(STAT2['TF_TG_SCORE'], 95)
    STAT2_median=sts.median(STAT2['TF_TG_SCORE'])
if(STAT2_median > background_median and STAT2_mannwhitneyU['pvalue']<0.01):
    background_STAT2 = pd.concat([background,STAT2],axis=0)
    ax_STAT2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_STAT2,palette="Set3")
    ax_STAT2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/STAT2.png')
    del background_STAT2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['STAT2',STAT2_sum,STAT2_length,STAT2_mean,STAT2_median,STAT2_quantile_95,STAT2_quantile]
    row_counter=row_counter+1
del STAT2
plt.figure(figsize=(20, 17))


STAT1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/STAT1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
STAT1.columns=['TF_TG_SCORE','label']
STAT1_sum = sum(STAT1['TF_TG_SCORE'])
STAT1_length = len(STAT1)
STAT1_mean=0
STAT1_quantile=0
STAT1_quantile_95=0
STAT1_median=0
STAT1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],STAT1['TF_TG_SCORE'])
if(STAT1_length>0):
    STAT1_mean=STAT1_sum/STAT1_length
    STAT1_quantile=np.percentile(STAT1['TF_TG_SCORE'], 99)
    STAT1_quantile_95=np.percentile(STAT1['TF_TG_SCORE'], 95)
    STAT1_median=sts.median(STAT1['TF_TG_SCORE'])
if(STAT1_median > background_median and STAT1_mannwhitneyU['pvalue']<0.01):
    background_STAT1 = pd.concat([background,STAT1],axis=0)
    ax_STAT1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_STAT1,palette="Set3")
    ax_STAT1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/STAT1.png')
    del background_STAT1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['STAT1',STAT1_sum,STAT1_length,STAT1_mean,STAT1_median,STAT1_quantile_95,STAT1_quantile]
    row_counter=row_counter+1
del STAT1
plt.figure(figsize=(20, 17))


E2F6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/E2F6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
E2F6.columns=['TF_TG_SCORE','label']
E2F6_sum = sum(E2F6['TF_TG_SCORE'])
E2F6_length = len(E2F6)
E2F6_mean=0
E2F6_quantile=0
E2F6_quantile_95=0
E2F6_median=0
E2F6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],E2F6['TF_TG_SCORE'])
if(E2F6_length>0):
    E2F6_mean=E2F6_sum/E2F6_length
    E2F6_quantile=np.percentile(E2F6['TF_TG_SCORE'], 99)
    E2F6_quantile_95=np.percentile(E2F6['TF_TG_SCORE'], 95)
    E2F6_median=sts.median(E2F6['TF_TG_SCORE'])
if(E2F6_median > background_median and E2F6_mannwhitneyU['pvalue']<0.01):
    background_E2F6 = pd.concat([background,E2F6],axis=0)
    ax_E2F6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_E2F6,palette="Set3")
    ax_E2F6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/E2F6.png')
    del background_E2F6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['E2F6',E2F6_sum,E2F6_length,E2F6_mean,E2F6_median,E2F6_quantile_95,E2F6_quantile]
    row_counter=row_counter+1
del E2F6
plt.figure(figsize=(20, 17))


EGR1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/EGR1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
EGR1.columns=['TF_TG_SCORE','label']
EGR1_sum = sum(EGR1['TF_TG_SCORE'])
EGR1_length = len(EGR1)
EGR1_mean=0
EGR1_quantile=0
EGR1_quantile_95=0
EGR1_median=0
EGR1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],EGR1['TF_TG_SCORE'])
if(EGR1_length>0):
    EGR1_mean=EGR1_sum/EGR1_length
    EGR1_quantile=np.percentile(EGR1['TF_TG_SCORE'], 99)
    EGR1_quantile_95=np.percentile(EGR1['TF_TG_SCORE'], 95)
    EGR1_median=sts.median(EGR1['TF_TG_SCORE'])
if(EGR1_median > background_median and EGR1_mannwhitneyU['pvalue']<0.01):
    background_EGR1 = pd.concat([background,EGR1],axis=0)
    ax_EGR1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_EGR1,palette="Set3")
    ax_EGR1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/EGR1.png')
    del background_EGR1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['EGR1',EGR1_sum,EGR1_length,EGR1_mean,EGR1_median,EGR1_quantile_95,EGR1_quantile]
    row_counter=row_counter+1
del EGR1
plt.figure(figsize=(20, 17))


SIX5=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SIX5_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SIX5.columns=['TF_TG_SCORE','label']
SIX5_sum = sum(SIX5['TF_TG_SCORE'])
SIX5_length = len(SIX5)
SIX5_mean=0
SIX5_quantile=0
SIX5_quantile_95=0
SIX5_median=0
SIX5_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SIX5['TF_TG_SCORE'])
if(SIX5_length>0):
    SIX5_mean=SIX5_sum/SIX5_length
    SIX5_quantile=np.percentile(SIX5['TF_TG_SCORE'], 99)
    SIX5_quantile_95=np.percentile(SIX5['TF_TG_SCORE'], 95)
    SIX5_median=sts.median(SIX5['TF_TG_SCORE'])
if(SIX5_median > background_median and SIX5_mannwhitneyU['pvalue']<0.01):
    background_SIX5 = pd.concat([background,SIX5],axis=0)
    ax_SIX5 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SIX5,palette="Set3")
    ax_SIX5.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SIX5.png')
    del background_SIX5
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SIX5',SIX5_sum,SIX5_length,SIX5_mean,SIX5_median,SIX5_quantile_95,SIX5_quantile]
    row_counter=row_counter+1
del SIX5
plt.figure(figsize=(20, 17))


CEBPG=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/CEBPG_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
CEBPG.columns=['TF_TG_SCORE','label']
CEBPG_sum = sum(CEBPG['TF_TG_SCORE'])
CEBPG_length = len(CEBPG)
CEBPG_mean=0
CEBPG_quantile=0
CEBPG_quantile_95=0
CEBPG_median=0
CEBPG_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],CEBPG['TF_TG_SCORE'])
if(CEBPG_length>0):
    CEBPG_mean=CEBPG_sum/CEBPG_length
    CEBPG_quantile=np.percentile(CEBPG['TF_TG_SCORE'], 99)
    CEBPG_quantile_95=np.percentile(CEBPG['TF_TG_SCORE'], 95)
    CEBPG_median=sts.median(CEBPG['TF_TG_SCORE'])
if(CEBPG_median > background_median and CEBPG_mannwhitneyU['pvalue']<0.01):
    background_CEBPG = pd.concat([background,CEBPG],axis=0)
    ax_CEBPG = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_CEBPG,palette="Set3")
    ax_CEBPG.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/CEBPG.png')
    del background_CEBPG
    plt.clf()
    df_interesting_stats.loc[row_counter]=['CEBPG',CEBPG_sum,CEBPG_length,CEBPG_mean,CEBPG_median,CEBPG_quantile_95,CEBPG_quantile]
    row_counter=row_counter+1
del CEBPG
plt.figure(figsize=(20, 17))


TGIF1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TGIF1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TGIF1.columns=['TF_TG_SCORE','label']
TGIF1_sum = sum(TGIF1['TF_TG_SCORE'])
TGIF1_length = len(TGIF1)
TGIF1_mean=0
TGIF1_quantile=0
TGIF1_quantile_95=0
TGIF1_median=0
TGIF1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TGIF1['TF_TG_SCORE'])
if(TGIF1_length>0):
    TGIF1_mean=TGIF1_sum/TGIF1_length
    TGIF1_quantile=np.percentile(TGIF1['TF_TG_SCORE'], 99)
    TGIF1_quantile_95=np.percentile(TGIF1['TF_TG_SCORE'], 95)
    TGIF1_median=sts.median(TGIF1['TF_TG_SCORE'])
if(TGIF1_median > background_median and TGIF1_mannwhitneyU['pvalue']<0.01):
    background_TGIF1 = pd.concat([background,TGIF1],axis=0)
    ax_TGIF1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TGIF1,palette="Set3")
    ax_TGIF1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TGIF1.png')
    del background_TGIF1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TGIF1',TGIF1_sum,TGIF1_length,TGIF1_mean,TGIF1_median,TGIF1_quantile_95,TGIF1_quantile]
    row_counter=row_counter+1
del TGIF1
plt.figure(figsize=(20, 17))


PBX2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PBX2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PBX2.columns=['TF_TG_SCORE','label']
PBX2_sum = sum(PBX2['TF_TG_SCORE'])
PBX2_length = len(PBX2)
PBX2_mean=0
PBX2_quantile=0
PBX2_quantile_95=0
PBX2_median=0
PBX2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PBX2['TF_TG_SCORE'])
if(PBX2_length>0):
    PBX2_mean=PBX2_sum/PBX2_length
    PBX2_quantile=np.percentile(PBX2['TF_TG_SCORE'], 99)
    PBX2_quantile_95=np.percentile(PBX2['TF_TG_SCORE'], 95)
    PBX2_median=sts.median(PBX2['TF_TG_SCORE'])
if(PBX2_median > background_median and PBX2_mannwhitneyU['pvalue']<0.01):
    background_PBX2 = pd.concat([background,PBX2],axis=0)
    ax_PBX2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PBX2,palette="Set3")
    ax_PBX2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PBX2.png')
    del background_PBX2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PBX2',PBX2_sum,PBX2_length,PBX2_mean,PBX2_median,PBX2_quantile_95,PBX2_quantile]
    row_counter=row_counter+1
del PBX2
plt.figure(figsize=(20, 17))


ZFP57=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ZFP57_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ZFP57.columns=['TF_TG_SCORE','label']
ZFP57_sum = sum(ZFP57['TF_TG_SCORE'])
ZFP57_length = len(ZFP57)
ZFP57_mean=0
ZFP57_quantile=0
ZFP57_quantile_95=0
ZFP57_median=0
ZFP57_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ZFP57['TF_TG_SCORE'])
if(ZFP57_length>0):
    ZFP57_mean=ZFP57_sum/ZFP57_length
    ZFP57_quantile=np.percentile(ZFP57['TF_TG_SCORE'], 99)
    ZFP57_quantile_95=np.percentile(ZFP57['TF_TG_SCORE'], 95)
    ZFP57_median=sts.median(ZFP57['TF_TG_SCORE'])
if(ZFP57_median > background_median and ZFP57_mannwhitneyU['pvalue']<0.01):
    background_ZFP57 = pd.concat([background,ZFP57],axis=0)
    ax_ZFP57 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ZFP57,palette="Set3")
    ax_ZFP57.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ZFP57.png')
    del background_ZFP57
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ZFP57',ZFP57_sum,ZFP57_length,ZFP57_mean,ZFP57_median,ZFP57_quantile_95,ZFP57_quantile]
    row_counter=row_counter+1
del ZFP57
plt.figure(figsize=(20, 17))


BCL6=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/BCL6_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
BCL6.columns=['TF_TG_SCORE','label']
BCL6_sum = sum(BCL6['TF_TG_SCORE'])
BCL6_length = len(BCL6)
BCL6_mean=0
BCL6_quantile=0
BCL6_quantile_95=0
BCL6_median=0
BCL6_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],BCL6['TF_TG_SCORE'])
if(BCL6_length>0):
    BCL6_mean=BCL6_sum/BCL6_length
    BCL6_quantile=np.percentile(BCL6['TF_TG_SCORE'], 99)
    BCL6_quantile_95=np.percentile(BCL6['TF_TG_SCORE'], 95)
    BCL6_median=sts.median(BCL6['TF_TG_SCORE'])
if(BCL6_median > background_median and BCL6_mannwhitneyU['pvalue']<0.01):
    background_BCL6 = pd.concat([background,BCL6],axis=0)
    ax_BCL6 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_BCL6,palette="Set3")
    ax_BCL6.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/BCL6.png')
    del background_BCL6
    plt.clf()
    df_interesting_stats.loc[row_counter]=['BCL6',BCL6_sum,BCL6_length,BCL6_mean,BCL6_median,BCL6_quantile_95,BCL6_quantile]
    row_counter=row_counter+1
del BCL6
plt.figure(figsize=(20, 17))


SOX17=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/SOX17_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
SOX17.columns=['TF_TG_SCORE','label']
SOX17_sum = sum(SOX17['TF_TG_SCORE'])
SOX17_length = len(SOX17)
SOX17_mean=0
SOX17_quantile=0
SOX17_quantile_95=0
SOX17_median=0
SOX17_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],SOX17['TF_TG_SCORE'])
if(SOX17_length>0):
    SOX17_mean=SOX17_sum/SOX17_length
    SOX17_quantile=np.percentile(SOX17['TF_TG_SCORE'], 99)
    SOX17_quantile_95=np.percentile(SOX17['TF_TG_SCORE'], 95)
    SOX17_median=sts.median(SOX17['TF_TG_SCORE'])
if(SOX17_median > background_median and SOX17_mannwhitneyU['pvalue']<0.01):
    background_SOX17 = pd.concat([background,SOX17],axis=0)
    ax_SOX17 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_SOX17,palette="Set3")
    ax_SOX17.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/SOX17.png')
    del background_SOX17
    plt.clf()
    df_interesting_stats.loc[row_counter]=['SOX17',SOX17_sum,SOX17_length,SOX17_mean,SOX17_median,SOX17_quantile_95,SOX17_quantile]
    row_counter=row_counter+1
del SOX17
plt.figure(figsize=(20, 17))


NR4A2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR4A2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR4A2.columns=['TF_TG_SCORE','label']
NR4A2_sum = sum(NR4A2['TF_TG_SCORE'])
NR4A2_length = len(NR4A2)
NR4A2_mean=0
NR4A2_quantile=0
NR4A2_quantile_95=0
NR4A2_median=0
NR4A2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR4A2['TF_TG_SCORE'])
if(NR4A2_length>0):
    NR4A2_mean=NR4A2_sum/NR4A2_length
    NR4A2_quantile=np.percentile(NR4A2['TF_TG_SCORE'], 99)
    NR4A2_quantile_95=np.percentile(NR4A2['TF_TG_SCORE'], 95)
    NR4A2_median=sts.median(NR4A2['TF_TG_SCORE'])
if(NR4A2_median > background_median and NR4A2_mannwhitneyU['pvalue']<0.01):
    background_NR4A2 = pd.concat([background,NR4A2],axis=0)
    ax_NR4A2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR4A2,palette="Set3")
    ax_NR4A2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR4A2.png')
    del background_NR4A2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR4A2',NR4A2_sum,NR4A2_length,NR4A2_mean,NR4A2_median,NR4A2_quantile_95,NR4A2_quantile]
    row_counter=row_counter+1
del NR4A2
plt.figure(figsize=(20, 17))


FOS=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/FOS..JUN_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
FOS.columns=['TF_TG_SCORE','label']
FOS_sum = sum(FOS['TF_TG_SCORE'])
FOS_length = len(FOS)
FOS_mean=0
FOS_quantile=0
FOS_quantile_95=0
FOS_median=0
FOS_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],FOS['TF_TG_SCORE'])
if(FOS_length>0):
    FOS_mean=FOS_sum/FOS_length
    FOS_quantile=np.percentile(FOS['TF_TG_SCORE'], 99)
    FOS_quantile_95=np.percentile(FOS['TF_TG_SCORE'], 95)
    FOS_median=sts.median(FOS['TF_TG_SCORE'])
if(FOS_median > background_median and FOS_mannwhitneyU['pvalue']<0.01):
    background_FOS = pd.concat([background,FOS],axis=0)
    ax_FOS = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_FOS,palette="Set3")
    ax_FOS.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/FOS..JUN.png')
    del background_FOS
    plt.clf()
    df_interesting_stats.loc[row_counter]=['FOS..JUN',FOS_sum,FOS_length,FOS_mean,FOS_median,FOS_quantile_95,FOS_quantile]
    row_counter=row_counter+1
del FOS
plt.figure(figsize=(20, 17))


TCF12=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/TCF12_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
TCF12.columns=['TF_TG_SCORE','label']
TCF12_sum = sum(TCF12['TF_TG_SCORE'])
TCF12_length = len(TCF12)
TCF12_mean=0
TCF12_quantile=0
TCF12_quantile_95=0
TCF12_median=0
TCF12_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],TCF12['TF_TG_SCORE'])
if(TCF12_length>0):
    TCF12_mean=TCF12_sum/TCF12_length
    TCF12_quantile=np.percentile(TCF12['TF_TG_SCORE'], 99)
    TCF12_quantile_95=np.percentile(TCF12['TF_TG_SCORE'], 95)
    TCF12_median=sts.median(TCF12['TF_TG_SCORE'])
if(TCF12_median > background_median and TCF12_mannwhitneyU['pvalue']<0.01):
    background_TCF12 = pd.concat([background,TCF12],axis=0)
    ax_TCF12 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_TCF12,palette="Set3")
    ax_TCF12.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/TCF12.png')
    del background_TCF12
    plt.clf()
    df_interesting_stats.loc[row_counter]=['TCF12',TCF12_sum,TCF12_length,TCF12_mean,TCF12_median,TCF12_quantile_95,TCF12_quantile]
    row_counter=row_counter+1
del TCF12
plt.figure(figsize=(20, 17))


PRDM5=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PRDM5_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PRDM5.columns=['TF_TG_SCORE','label']
PRDM5_sum = sum(PRDM5['TF_TG_SCORE'])
PRDM5_length = len(PRDM5)
PRDM5_mean=0
PRDM5_quantile=0
PRDM5_quantile_95=0
PRDM5_median=0
PRDM5_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PRDM5['TF_TG_SCORE'])
if(PRDM5_length>0):
    PRDM5_mean=PRDM5_sum/PRDM5_length
    PRDM5_quantile=np.percentile(PRDM5['TF_TG_SCORE'], 99)
    PRDM5_quantile_95=np.percentile(PRDM5['TF_TG_SCORE'], 95)
    PRDM5_median=sts.median(PRDM5['TF_TG_SCORE'])
if(PRDM5_median > background_median and PRDM5_mannwhitneyU['pvalue']<0.01):
    background_PRDM5 = pd.concat([background,PRDM5],axis=0)
    ax_PRDM5 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PRDM5,palette="Set3")
    ax_PRDM5.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PRDM5.png')
    del background_PRDM5
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PRDM5',PRDM5_sum,PRDM5_length,PRDM5_mean,PRDM5_median,PRDM5_quantile_95,PRDM5_quantile]
    row_counter=row_counter+1
del PRDM5
plt.figure(figsize=(20, 17))


MEF2D=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MEF2D_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MEF2D.columns=['TF_TG_SCORE','label']
MEF2D_sum = sum(MEF2D['TF_TG_SCORE'])
MEF2D_length = len(MEF2D)
MEF2D_mean=0
MEF2D_quantile=0
MEF2D_quantile_95=0
MEF2D_median=0
MEF2D_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MEF2D['TF_TG_SCORE'])
if(MEF2D_length>0):
    MEF2D_mean=MEF2D_sum/MEF2D_length
    MEF2D_quantile=np.percentile(MEF2D['TF_TG_SCORE'], 99)
    MEF2D_quantile_95=np.percentile(MEF2D['TF_TG_SCORE'], 95)
    MEF2D_median=sts.median(MEF2D['TF_TG_SCORE'])
if(MEF2D_median > background_median and MEF2D_mannwhitneyU['pvalue']<0.01):
    background_MEF2D = pd.concat([background,MEF2D],axis=0)
    ax_MEF2D = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MEF2D,palette="Set3")
    ax_MEF2D.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MEF2D.png')
    del background_MEF2D
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MEF2D',MEF2D_sum,MEF2D_length,MEF2D_mean,MEF2D_median,MEF2D_quantile_95,MEF2D_quantile]
    row_counter=row_counter+1
del MEF2D
plt.figure(figsize=(20, 17))


ARNTL=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/ARNTL_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
ARNTL.columns=['TF_TG_SCORE','label']
ARNTL_sum = sum(ARNTL['TF_TG_SCORE'])
ARNTL_length = len(ARNTL)
ARNTL_mean=0
ARNTL_quantile=0
ARNTL_quantile_95=0
ARNTL_median=0
ARNTL_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],ARNTL['TF_TG_SCORE'])
if(ARNTL_length>0):
    ARNTL_mean=ARNTL_sum/ARNTL_length
    ARNTL_quantile=np.percentile(ARNTL['TF_TG_SCORE'], 99)
    ARNTL_quantile_95=np.percentile(ARNTL['TF_TG_SCORE'], 95)
    ARNTL_median=sts.median(ARNTL['TF_TG_SCORE'])
if(ARNTL_median > background_median and ARNTL_mannwhitneyU['pvalue']<0.01):
    background_ARNTL = pd.concat([background,ARNTL],axis=0)
    ax_ARNTL = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_ARNTL,palette="Set3")
    ax_ARNTL.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/ARNTL.png')
    del background_ARNTL
    plt.clf()
    df_interesting_stats.loc[row_counter]=['ARNTL',ARNTL_sum,ARNTL_length,ARNTL_mean,ARNTL_median,ARNTL_quantile_95,ARNTL_quantile]
    row_counter=row_counter+1
del ARNTL
plt.figure(figsize=(20, 17))


AHR=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/AHR..ARNT_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
AHR.columns=['TF_TG_SCORE','label']
AHR_sum = sum(AHR['TF_TG_SCORE'])
AHR_length = len(AHR)
AHR_mean=0
AHR_quantile=0
AHR_quantile_95=0
AHR_median=0
AHR_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],AHR['TF_TG_SCORE'])
if(AHR_length>0):
    AHR_mean=AHR_sum/AHR_length
    AHR_quantile=np.percentile(AHR['TF_TG_SCORE'], 99)
    AHR_quantile_95=np.percentile(AHR['TF_TG_SCORE'], 95)
    AHR_median=sts.median(AHR['TF_TG_SCORE'])
if(AHR_median > background_median and AHR_mannwhitneyU['pvalue']<0.01):
    background_AHR = pd.concat([background,AHR],axis=0)
    ax_AHR = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_AHR,palette="Set3")
    ax_AHR.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/AHR..ARNT.png')
    del background_AHR
    plt.clf()
    df_interesting_stats.loc[row_counter]=['AHR..ARNT',AHR_sum,AHR_length,AHR_mean,AHR_median,AHR_quantile_95,AHR_quantile]
    row_counter=row_counter+1
del AHR
plt.figure(figsize=(20, 17))


DBP=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/DBP_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
DBP.columns=['TF_TG_SCORE','label']
DBP_sum = sum(DBP['TF_TG_SCORE'])
DBP_length = len(DBP)
DBP_mean=0
DBP_quantile=0
DBP_quantile_95=0
DBP_median=0
DBP_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],DBP['TF_TG_SCORE'])
if(DBP_length>0):
    DBP_mean=DBP_sum/DBP_length
    DBP_quantile=np.percentile(DBP['TF_TG_SCORE'], 99)
    DBP_quantile_95=np.percentile(DBP['TF_TG_SCORE'], 95)
    DBP_median=sts.median(DBP['TF_TG_SCORE'])
if(DBP_median > background_median and DBP_mannwhitneyU['pvalue']<0.01):
    background_DBP = pd.concat([background,DBP],axis=0)
    ax_DBP = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_DBP,palette="Set3")
    ax_DBP.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/DBP.png')
    del background_DBP
    plt.clf()
    df_interesting_stats.loc[row_counter]=['DBP',DBP_sum,DBP_length,DBP_mean,DBP_median,DBP_quantile_95,DBP_quantile]
    row_counter=row_counter+1
del DBP
plt.figure(figsize=(20, 17))


PRDM9=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/PRDM9_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
PRDM9.columns=['TF_TG_SCORE','label']
PRDM9_sum = sum(PRDM9['TF_TG_SCORE'])
PRDM9_length = len(PRDM9)
PRDM9_mean=0
PRDM9_quantile=0
PRDM9_quantile_95=0
PRDM9_median=0
PRDM9_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],PRDM9['TF_TG_SCORE'])
if(PRDM9_length>0):
    PRDM9_mean=PRDM9_sum/PRDM9_length
    PRDM9_quantile=np.percentile(PRDM9['TF_TG_SCORE'], 99)
    PRDM9_quantile_95=np.percentile(PRDM9['TF_TG_SCORE'], 95)
    PRDM9_median=sts.median(PRDM9['TF_TG_SCORE'])
if(PRDM9_median > background_median and PRDM9_mannwhitneyU['pvalue']<0.01):
    background_PRDM9 = pd.concat([background,PRDM9],axis=0)
    ax_PRDM9 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_PRDM9,palette="Set3")
    ax_PRDM9.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/PRDM9.png')
    del background_PRDM9
    plt.clf()
    df_interesting_stats.loc[row_counter]=['PRDM9',PRDM9_sum,PRDM9_length,PRDM9_mean,PRDM9_median,PRDM9_quantile_95,PRDM9_quantile]
    row_counter=row_counter+1
del PRDM9
plt.figure(figsize=(20, 17))


NFIL3=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NFIL3_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NFIL3.columns=['TF_TG_SCORE','label']
NFIL3_sum = sum(NFIL3['TF_TG_SCORE'])
NFIL3_length = len(NFIL3)
NFIL3_mean=0
NFIL3_quantile=0
NFIL3_quantile_95=0
NFIL3_median=0
NFIL3_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NFIL3['TF_TG_SCORE'])
if(NFIL3_length>0):
    NFIL3_mean=NFIL3_sum/NFIL3_length
    NFIL3_quantile=np.percentile(NFIL3['TF_TG_SCORE'], 99)
    NFIL3_quantile_95=np.percentile(NFIL3['TF_TG_SCORE'], 95)
    NFIL3_median=sts.median(NFIL3['TF_TG_SCORE'])
if(NFIL3_median > background_median and NFIL3_mannwhitneyU['pvalue']<0.01):
    background_NFIL3 = pd.concat([background,NFIL3],axis=0)
    ax_NFIL3 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NFIL3,palette="Set3")
    ax_NFIL3.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NFIL3.png')
    del background_NFIL3
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NFIL3',NFIL3_sum,NFIL3_length,NFIL3_mean,NFIL3_median,NFIL3_quantile_95,NFIL3_quantile]
    row_counter=row_counter+1
del NFIL3
plt.figure(figsize=(20, 17))


MEIS2=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/MEIS2_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
MEIS2.columns=['TF_TG_SCORE','label']
MEIS2_sum = sum(MEIS2['TF_TG_SCORE'])
MEIS2_length = len(MEIS2)
MEIS2_mean=0
MEIS2_quantile=0
MEIS2_quantile_95=0
MEIS2_median=0
MEIS2_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],MEIS2['TF_TG_SCORE'])
if(MEIS2_length>0):
    MEIS2_mean=MEIS2_sum/MEIS2_length
    MEIS2_quantile=np.percentile(MEIS2['TF_TG_SCORE'], 99)
    MEIS2_quantile_95=np.percentile(MEIS2['TF_TG_SCORE'], 95)
    MEIS2_median=sts.median(MEIS2['TF_TG_SCORE'])
if(MEIS2_median > background_median and MEIS2_mannwhitneyU['pvalue']<0.01):
    background_MEIS2 = pd.concat([background,MEIS2],axis=0)
    ax_MEIS2 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_MEIS2,palette="Set3")
    ax_MEIS2.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/MEIS2.png')
    del background_MEIS2
    plt.clf()
    df_interesting_stats.loc[row_counter]=['MEIS2',MEIS2_sum,MEIS2_length,MEIS2_mean,MEIS2_median,MEIS2_quantile_95,MEIS2_quantile]
    row_counter=row_counter+1
del MEIS2
plt.figure(figsize=(20, 17))


NR1D1=pd.read_table('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/02_TF_TG_SCORES/02_TF_DISTRIBUTION/01_ALL/NR1D1_distribution.csv', comment="#", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)
NR1D1.columns=['TF_TG_SCORE','label']
NR1D1_sum = sum(NR1D1['TF_TG_SCORE'])
NR1D1_length = len(NR1D1)
NR1D1_mean=0
NR1D1_quantile=0
NR1D1_quantile_95=0
NR1D1_median=0
NR1D1_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE'],NR1D1['TF_TG_SCORE'])
if(NR1D1_length>0):
    NR1D1_mean=NR1D1_sum/NR1D1_length
    NR1D1_quantile=np.percentile(NR1D1['TF_TG_SCORE'], 99)
    NR1D1_quantile_95=np.percentile(NR1D1['TF_TG_SCORE'], 95)
    NR1D1_median=sts.median(NR1D1['TF_TG_SCORE'])
if(NR1D1_median > background_median and NR1D1_mannwhitneyU['pvalue']<0.01):
    background_NR1D1 = pd.concat([background,NR1D1],axis=0)
    ax_NR1D1 = sns.boxplot(x="label", y="TF_TG_SCORE",data=background_NR1D1,palette="Set3")
    ax_NR1D1.set_yscale("log")
    plt.savefig(f'/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/04_DISTR_PLOTS/01_ALL/NR1D1.png')
    del background_NR1D1
    plt.clf()
    df_interesting_stats.loc[row_counter]=['NR1D1',NR1D1_sum,NR1D1_length,NR1D1_mean,NR1D1_median,NR1D1_quantile_95,NR1D1_quantile]
    row_counter=row_counter+1
del NR1D1
plt.figure(figsize=(20, 17))


df_interesting_stats.to_csv('/home/markus/data/COM2POSE/01_COM2POSE/06_LIKE_03_COMBINED_TGENE_TEPIC_COPY/working_dir/07_DISTRIBUTION_ANALYSIS/05_STATS/01_ALL/stats.csv',sep='	')