#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("data.table", quietly = TRUE))
  BiocManager::install("data.table")

if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")

if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")

if (!requireNamespace("gplots", quietly = TRUE))
  BiocManager::install("gplots")

if (!requireNamespace("glmnet", quietly = TRUE))
  BiocManager::install("glmnet")

if (!requireNamespace("doMC", quietly = TRUE))
  BiocManager::install("doMC")

if (!requireNamespace("methods", quietly = TRUE))
  BiocManager::install("methods")

 if (!requireNamespace("DiffLogo", quietly = TRUE))
   BiocManager::install("DiffLogo")

 if (!requireNamespace("seqLogo", quietly = TRUE))
   BiocManager::install("seqLogo")

install.packages("pheatmap")

BiocManager::install(version = "3.13", force = TRUE, ask = FALSE)