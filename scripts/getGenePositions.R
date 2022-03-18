if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library('biomaRt')
httr::set_config(httr::config(ssl_verifypeer = FALSE))

input <- read.csv('{INPUTFILE}', sep = '\t')
input_ensg <- input$ensembl_gene_id
not_done <- TRUE

G_list <- data.frame()

while (not_done)
{
  tryCatch({
    mart <- useDataset("{SPECIES}", useMart("ensembl"))
    G_list_intern <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "band"), values = input_ensg, mart = mart)
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

G_list <- merge(G_list, input, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
G_list <- G_list[, c(7, 2, 3, 4, 5, 6, 1)]
write.table(G_list, "{DATA_PREV_FILE}", row.names = FALSE, quote = F, sep = "\t")
not_done <- TRUE
while (not_done) {
  tryCatch({
    #check which genome version we have in biomart

    ensembl <- useEnsembl(biomart = "genes")
    datasets <- listDatasets(ensembl)
    used_dataset <- searchDatasets(mart = ensembl, pattern = "{SPECIES}")
    version <- used_dataset[, 3]
    version <- as.character(version)
    fileConn <- file("{VERSIONFILE}")
    writeLines(version, fileConn)

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

