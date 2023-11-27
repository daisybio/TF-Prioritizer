#!/usr/bin/Rscript

install.packages("devtools")
BiocManager::install(c("IRanges", "GenomicRanges", "bamsignals", "rtracklayer", "Rsamtools", "edgeR", "affyPLM"))
