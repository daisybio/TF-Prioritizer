library(DESeq2)

##TODO: one of these groups must be active
# L1_vs_p6
#metadata_df = data.frame(sample_id = c("L1_WT1", "L1_WT2", "L1_WT4", "p6_WT1", "p6_WT2", "p6_WT4"), group = c("L1", "L1", "L1", "p6", "p6", "p6"))
#input_groups = 'L1_vs_p6'

# L1_vs_p13
#metadata_df = data.frame(sample_id = c("L1_WT1", "L1_WT2", "L1_WT4", "p13_WT1", "p13_WT2", "p13_WT3", "p13_WT4"), group = c("L1", "L1", "L1", "p13", "p13", "p13", "p13"))
#input_groups = 'L1_vs_p13'

# L10_vs_p13
#metadata_df = data.frame(sample_id = c("L10_WT1", "L10_WT2", "L10_WT3", "L10_WT4", "p13_WT1", "p13_WT2", "p13_WT3", "p13_WT4"), group = c("L10", "L10", "L10", "L10", "p13", "p13", "p13", "p13"))
#input_groups = 'L10_vs_p13'

# L10_vs_p6
metadata_df = data.frame(sample_id = c("L10_WT1", "L10_WT2", "L10_WT3", "L10_WT4", "p6_WT1", "p6_WT2", "p6_WT4"), group = c("L10", "L10", "L10", "L10", "p6", "p6", "p6"))
input_groups = 'L10_vs_p6'

# Lactation vs pregnancy
#metadata_df = data.frame(sample_id = c("L10_WT1", "L10_WT2", "L10_WT3", "L10_WT4", "L1_WT1", "L1_WT2", "L1_WT4", "p13_WT1", "p13_WT2", "p13_WT3", "p13_WT4", "p6_WT1", "p6_WT2", "p6_WT4"), group = c("Lactation", "Lactation", "Lactation", "Lactation", "Lactation", "Lactation", "Lactation", "Pregnancy", "Pregnancy", "Pregnancy", "Pregnancy", "Pregnancy", "Pregnancy", "Pregnancy"))
#input_groups = 'Lactation_vs_Pregnancy'



group_one = strsplit(input_groups, "_")[[1]][1]
group_two = strsplit(input_groups, "_")[[1]][3]
rownames(metadata_df) <- metadata_df$sample_id
metadata_df$sample_id <- NULL

count_path = paste("C:\\Users\\Marku\\Dropbox\\UNI\\Promotion_Projekte\\PhD_work\\B_Documents\\KevinsBachelorThesis\\02_Analysis\\DESeq2\\output\\", input_groups, ".tsv", sep="")

count_df = read.csv(count_path, sep = "\t", header = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData=count_df, 
                              colData=metadata_df, 
                              design=~group)

##TODO: one of these must be active or another threshold must be added
#threshold = 4
#threshold = 100
threshold = 1000

keep <- rowSums(counts(dds)) >= threshold
dds <- dds[keep,]


output_path = paste("C:\\Users\\Marku\\Dropbox\\UNI\\Promotion_Projekte\\PhD_work\\second_step\\deseq2\\output_raw\\", input_groups,"_thresholdCounts_",threshold, ".tsv", sep="")

dds <- DESeq(dds)
res <- results(dds)
summary(res)
res
write.table(res,file=output_path,sep="\t")
