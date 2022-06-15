#!/usr/bin/env Rscript

install.packages(c("BiocManager", "pheatmap"))
BiocManager::install(c("biomaRt", "data.table", "DESeq2", "ggplot2", "dplyr", "gplots", "doMC", "methods",
                       "DiffLogo", "seqLogo"), version = "3.15")