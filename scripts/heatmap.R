library(DESeq2)
library(pheatmap)

require(dplyr)

target_genes_path <- "{TARGET_GENES_PATH}"
map_path <- "{MAP_PATH}"
preprocessing_path <- "{PREPROCESSING_PATH}"
heatmap_dir <- "{HEATMAP_DIR}"

messagesfile <- file(paste(heatmap_dir, "messages.Rout", sep = "/"), open = "wt")
sink(messagesfile, type = "message")

number_of_genes_to_select <- { NUMBER_OF_GENES }
map <- read.csv(map_path, sep = "\t")
map[, 2] <- toupper(map[, 2])

for (group in list.files(path = target_genes_path)) {
  hms <- list.files(paste(target_genes_path, group, sep = "/"))
  for (hm in hms) {
    genes <- list.files(paste(target_genes_path, group, hm, sep = "/"))
    for (geneFileName in genes) {
    gene <- tools::file_path_sans_ext(geneFileName)
      selected_genes <- read.csv(paste(target_genes_path, group, hm, geneFileName, sep = "/"),
                                 sep =
                                   "\t", nrows = number_of_genes_to_select)
      selected_genes <- transform(selected_genes, ENSEMBL_GENE_ID = map$ensembl_gene_id[match(TARGET_GENE, map[, 2])])

      group_pairings <- list.files(preprocessing_path)

      # Select only currently available groups
      group_pairings <- group_pairings[grepl(paste(paste("^", group, "_.*", sep = ""), paste(".*_", group, "$", sep = ""),
                                                   sep = "|"), group_pairings)]
      target_dir <- paste(heatmap_dir, gene, hm, sep = "/")

      if (!file.exists(target_dir))
      {
        dir.create(target_dir, recursive = TRUE)
      }

      for (group_pairing in group_pairings) {
        group_pairing <- tools::file_path_sans_ext(group_pairing)
        plot_file <- paste(target_dir, paste(group_pairing, "png", sep = "."), sep = "/")
        if (file.exists(plot_file)) next

        co_group <- sub(group, "", group_pairing)
        co_group <- sub("_", "", co_group)
        co_genes_file <- paste(target_genes_path, co_group, hm, geneFileName, sep = "/")
        if (!file.exists(co_genes_file)) next

        co_selected_genes <- read.csv(co_genes_file, sep = "\t", nrows = number_of_genes_to_select)
        co_selected_genes <- transform(co_selected_genes, ENSEMBL_GENE_ID = map$ensembl_gene_id[match(TARGET_GENE,
                                                                                                      map[, 2])])

        combined_selected_genes <- rbind(selected_genes, co_selected_genes)

        plot_file <- paste(target_dir, paste(group_pairing, "png", sep = "."), sep = "/")

        data <- read.csv(paste(preprocessing_path, paste(group_pairing, "tsv", sep = "."), sep = "/"),
                         sep = "\t")

        chosen <- subset(data, subset = data$Geneid %in% combined_selected_genes$ENSEMBL_GENE_ID)
        target_file <- paste(target_dir, paste(group_pairing, "csv", sep = "."), sep = "/")

        colnames(chosen) <- sub(".count.txt", "", colnames(chosen))

        if (!file.exists(target_file))
        {
          file.create(target_file)
        }

        # TODO: Find a cleaner solution
        # Drop the entries that contain zeros, since they make DESeq and pheatmap crash
        # Might result in a total number of chosen elements that is smaller than the defined one.
        read_counts <- subset(chosen, select = -c(Geneid))
        chosen <- chosen[apply(read_counts, 1, function(row) all(row > 0) & any(row > 1)),]


        read_counts <- subset(chosen, select = -c(Geneid))

        # Add geneSymbols to data
        chosen <- transform(chosen, geneSymbol = map[, 2][match(Geneid, map$ensembl_gene_id)])

        # Store data
        write.csv(chosen, target_file, row.names = FALSE)

        metaData <- data.frame(row.names = colnames(read_counts))

        # Remove redundance from colnames
        metaData$group <- gsub("_.*", "", colnames(read_counts))

        dds <- DESeqDataSetFromMatrix(countData = read_counts,
                                      colData = metaData,
                                      design = ~group)
        dds <- DESeq(dds, quiet = TRUE)

        geneCounts_normalized <- counts(dds)

        pheatmap(geneCounts_normalized, scale = "row", filename = plot_file, labels_row = chosen$geneSymbol, legend = TRUE, annotation_legend = TRUE, angle_col = 45, show_row_dend = FALSE, show_col_dend = FALSE, show_rownames = TRUE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_col = metaData)
      }
    }
  }
}
q()