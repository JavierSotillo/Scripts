# ============================
# Install & Load Packages
# ============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("biomaRt", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
BiocManager::install(pkgs, ask = FALSE, update = FALSE)
lapply(pkgs, library, character.only = TRUE)

# ============================
# Working Directory
# ============================
setwd("C:/Users/javier.sotillo/OneDrive - Instituto de Salud Carlos III/Ov_EV_omics/RNA_Seq/RNAseq_Ov_EVs_Matt/New_October/miRNA_targetting")

# ============================
# Read Input Gene List (HGNC symbols)
# ============================
MV_ids <- readLines("MVs.txt")
MV_ids <- unique(trimws(MV_ids)) # Clean input

# ============================
# Convert gene symbols → Entrez IDs
# ============================
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

conversion <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = MV_ids,
  mart = ensembl
)

entrez_ids <- unique(na.omit(conversion$entrezgene_id))

if(length(entrez_ids) == 0){
  stop("No se pudieron obtener IDs Entrez. Revisa si tus símbolos son humanos y están bien escritos.")
}

message("\n✅ Conversión completada. Genes convertidos: ", length(entrez_ids))

# ============================
# GO Enrichment: BP, MF, CC
# ============================
UPregBP <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2, readable = TRUE)

UPregMF <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2, readable = TRUE)

UPregCC <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2, readable = TRUE)

# ============================
# KEGG Enrichment
# ============================
UPregKEGG <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)
UPregKEGG <- setReadable(UPregKEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ============================
# Save Results as CSV
# ============================
write.csv(as.data.frame(UPregBP), "MVs_GO_BP_enrichment_human.csv", row.names = FALSE)
write.csv(as.data.frame(UPregMF), "MVs_GO_MF_enrichment_human.csv", row.names = FALSE)
write.csv(as.data.frame(UPregCC), "MVs_GO_CC_enrichment_human.csv", row.names = FALSE)
write.csv(as.data.frame(UPregKEGG), "MVs_KEGG_enrichment_human.csv", row.names = FALSE)

# ============================
# Plot Results (dotplot)
# ============================
Dot_BP <- dotplot(UPregBP, showCategory = 30, font.size = 10) + ggtitle("GO: Biological Process")
ggsave("MVs_GO_BP_enrichment_human.tiff", Dot_BP, dpi = 300, width = 8, height = 9, units = "in")

Dot_MF <- dotplot(UPregMF, showCategory = 30, font.size = 10) + ggtitle("GO: Molecular Function")
ggsave("MVs_GO_MF_enrichment_human.tiff", Dot_MF, dpi = 300, width = 8, height = 9, units = "in")

Dot_CC <- dotplot(UPregCC, showCategory = 30, font.size = 10) + ggtitle("GO: Cellular Component")
ggsave("MVs_GO_CC_enrichment_human.tiff", Dot_CC, dpi = 300, width = 8, height = 9, units = "in")

Dot_KEGG <- dotplot(UPregKEGG, showCategory = 30, font.size = 10) + ggtitle("KEGG Pathway Enrichment")
ggsave("MVs_KEGG_enrichment_human.tiff", Dot_KEGG, dpi = 300, width = 8, height = 9, units = "in")

cat("\n✅ GO enrichment completed and separated into BP / MF / CC.\n")
