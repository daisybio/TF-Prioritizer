#!/usr/bin/env Rscript

install.packages(c("BiocManager", "pheatmap", "argparse", "glmnet", "devtools"))
BiocManager::install(c("biomaRt", "data.table", "DESeq2", "ggplot2", "dplyr", "gplots", "doMC", "methods",
                       "DiffLogo", "seqLogo", "IRanges", "GenomicRanges", "bamsignals", "rtracklayer",
                       "Rsamtools", "edgeR", "affyPLM"), version = "3.14")
devtools::install_github("tobiaszehnder/ehmm")
