# =============================================================================
# GO enrichment (Biological Process and Molecular Function) for mouse proteins
#
# Input: plain-text file with one UniProt/Swiss-Prot accession per line.
# Accessions are mapped to Entrez Gene IDs through Ensembl BioMart (requires
# an internet connection) and tested with clusterProfiler::enrichGO.
#
# For human gene symbols, see go_kegg_enrichment_human_symbols.R.
#
# Usage: edit the settings below, then run  Rscript go_enrichment_mouse_uniprot.R
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("biomaRt", "clusterProfiler", "org.Mm.eg.db", "enrichplot", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) BiocManager::install(missing, ask = FALSE, update = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

# --- Settings ----------------------------------------------------------------
input_file    <- "uniprot_ids.txt"
output_prefix <- "GO_mouse"

# --- Map UniProt accessions to Entrez Gene IDs -------------------------------
uniprot_ids <- unique(trimws(readLines(input_file, warn = FALSE)))
uniprot_ids <- uniprot_ids[uniprot_ids != ""]

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
conversion <- getBM(attributes = c("uniprotswissprot", "entrezgene_id"),
                    filters = "uniprotswissprot",
                    values = uniprot_ids,
                    mart = ensembl)
entrez_ids <- unique(na.omit(conversion$entrezgene_id))
message(length(entrez_ids), " of ", length(uniprot_ids), " accessions mapped to Entrez IDs")
if (length(entrez_ids) == 0) stop("No accessions could be mapped. Check that they are mouse Swiss-Prot accessions.")

# --- Enrichment --------------------------------------------------------------
run_go <- function(ont) {
  enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
           ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff = 0.2, readable = TRUE)
}
go_bp <- run_go("BP")
go_mf <- run_go("MF")

write.csv(as.data.frame(go_bp), paste0(output_prefix, "_BP.csv"), row.names = FALSE)
write.csv(as.data.frame(go_mf), paste0(output_prefix, "_MF.csv"), row.names = FALSE)

# --- Plots -------------------------------------------------------------------
save_plot <- function(plot, suffix, width = 10, height = 12) {
  ggsave(paste0(output_prefix, "_", suffix, ".tiff"), plot = plot,
         dpi = 300, width = width, height = height, units = "in",
         device = "tiff", compression = "lzw")
}

save_plot(dotplot(go_bp, showCategory = 30, font.size = 10) +
            ggtitle("GO enrichment: Biological Process"), "BP_dotplot")
save_plot(dotplot(go_mf, showCategory = 30, font.size = 10) +
            ggtitle("GO enrichment: Molecular Function"), "MF_dotplot")
save_plot(barplot(go_bp, showCategory = 20, font.size = 10) +
            ggtitle("GO enrichment: Biological Process"), "BP_barplot", 8, 6)

# Optional plots for Biological Process (need at least two enriched terms)
if (nrow(as.data.frame(go_bp)) >= 2) {
  save_plot(cnetplot(go_bp, showCategory = 10), "BP_network", 10, 10)
  save_plot(emapplot(pairwise_termsim(go_bp), showCategory = 15), "BP_enrichment_map", 10, 10)
  save_plot(goplot(go_bp), "BP_goplot", 12, 10)
  save_plot(heatplot(go_bp, showCategory = 10), "BP_heatplot", 12, 6)
  save_plot(upsetplot(go_bp), "BP_upset", 10, 6)
}

message("GO enrichment completed")
