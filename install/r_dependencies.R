#!/usr/bin/env Rscript

install.packages(c("BiocManager", "pheatmap"))

BiocManager::install(c("biomaRt", "data.table", "DESeq2", "ggplot2", "dplyr", "gplots", "doMC", "methods",
                       "DiffLogo", "seqLogo"), version = "3.15")
install_version("glmnet", version = "4.1.3")
