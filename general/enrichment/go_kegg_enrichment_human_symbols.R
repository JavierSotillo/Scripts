# =============================================================================
# GO (BP, MF, CC) and KEGG enrichment for human genes
#
# Input: plain-text file with one HGNC gene symbol per line. Symbols are
# mapped to Entrez Gene IDs through Ensembl BioMart, and enrichment is run with
# clusterProfiler. BioMart and KEGG require an internet connection.
#
# Usage: edit the settings below, then run  Rscript go_kegg_enrichment_human_symbols.R
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("biomaRt", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) BiocManager::install(missing, ask = FALSE, update = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

# --- Settings ----------------------------------------------------------------
input_file    <- "gene_symbols.txt"
output_prefix <- "enrichment_human"

# --- Map gene symbols to Entrez Gene IDs -------------------------------------
symbols <- unique(trimws(readLines(input_file, warn = FALSE)))
symbols <- symbols[symbols != ""]

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                    filters = "hgnc_symbol",
                    values = symbols,
                    mart = ensembl)
entrez_ids <- unique(na.omit(conversion$entrezgene_id))
message(length(entrez_ids), " of ", length(symbols), " symbols mapped to Entrez IDs")
if (length(entrez_ids) == 0) stop("No symbols could be mapped. Check that they are human HGNC symbols.")

# --- Enrichment --------------------------------------------------------------
results <- list()
for (ont in c("BP", "MF", "CC")) {
  results[[ont]] <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID", ont = ont,
                             pAdjustMethod = "BH", pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2, readable = TRUE)
}
results[["KEGG"]] <- setReadable(
  enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05),
  OrgDb = org.Hs.eg.db, keyType = "ENTREZID"
)

titles <- c(BP = "GO: Biological Process", MF = "GO: Molecular Function",
            CC = "GO: Cellular Component", KEGG = "KEGG pathway enrichment")

# --- Save tables and dot plots -----------------------------------------------
for (name in names(results)) {
  write.csv(as.data.frame(results[[name]]),
            paste0(output_prefix, "_", name, ".csv"), row.names = FALSE)
  if (nrow(as.data.frame(results[[name]])) > 0) {
    p <- dotplot(results[[name]], showCategory = 30, font.size = 10) + ggtitle(titles[[name]])
    ggsave(paste0(output_prefix, "_", name, ".tiff"), p, dpi = 300,
           width = 8, height = 9, units = "in", compression = "lzw")
  }
}

message("GO and KEGG enrichment completed")
