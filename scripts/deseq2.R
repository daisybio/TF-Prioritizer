if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(DESeq2)

metadata_df <- data.frame(sample_id = c({ SAMPLES }), group = c({ GROUPS }))
input_groups <- "{COMBINATION}"
group_one <- strsplit(input_groups, "_")[[1]][1]
group_two <- strsplit(input_groups, "_")[[1]][2]
rownames(metadata_df) <- metadata_df$sample_id
metadata_df$sample_id <- NULL
count_path <- "{INPUTFILE}"
count_df <- read.csv(count_path, sep = "\t", header = T, row.names = NULL)
count_df <- count_df[!duplicated(count_df$Geneid),]
row.names(count_df) <- count_df[, 1]
count_df$Geneid <- NULL
dds <- DESeqDataSetFromMatrix(
  countData = count_df,
  colData = metadata_df,
  design = ~group
)
threshold <- 50
keep <- rowSums(counts(dds)) >= threshold
dds <- dds[keep,]
output_path <- "{OUTPUTFILE}"
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.table(res, file = output_path, sep = "\t")