if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library('biomaRt')
httr::set_config(httr::config(ssl_verifypeer = FALSE))
df <- read.csv('{INPUTFILE}')

not_done <- TRUE

G_list <- data.frame()

while (not_done)
{
  tryCatch({
    mart <- useDataset("{DATASET_SPECIES}", useMart("ensembl"))
    df$id <- NA
    G_list_intern <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "{SYMBOL_COLUMN}",
                                                                       "description"),
                           values = df$Geneid, mart = mart)
    G_list <- rbind(G_list, G_list_intern)
    not_done <- FALSE
  }, warning = function(w) {
    print("WARNING SECTION")
    print(w)
  }, error = function(e) {
    print("ERROR SECTION")
    print(e)
  }, finally = {
  })
}

write.table(G_list, "{OUTPUTFILE}", row.names = FALSE, quote = F, sep = "\t")
