#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
    BiocManager::install("biomaRt")

if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")