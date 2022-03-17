library("data.table")
library("EDASeq", include.only = "getGeneLengthAndGCContent")

httr::set_config(httr::config(ssl_verifypeer = FALSE))

geneIds <- c({ GENEIDs })

result <- data.frame()
result <- rbind(result, getGeneLengthAndGCContent(geneIds, "{BIOMART_SPECIES}"))

result <- setDT(result, keep.rownames = TRUE)[]
colnames(result) <- c("ENSG", "length", "gc")
write.table(result, file = "{OUTPUT_FILE}", sep = '\t', quote = FALSE, row.names = FALSE)
q(save = TRUE)