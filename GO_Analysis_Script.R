# Instala y carga paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("biomaRt", "clusterProfiler", "org.Mm.eg.db", "enrichplot")
BiocManager::install(pkgs, ask = FALSE, update = FALSE)
lapply(pkgs, library, character.only = TRUE)


library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

#Set working directory
setwd("C:/Users/javier.sotillo/Desktop")


# 1. Leer tu lista de UniProt IDs desde un archivo de texto
Up_ids <- readLines("Up_ids.txt")

# In case the table is a tab-delimited file
#Read tab-delimited file and extract UniProt ID column
#df <- read.table("uniprot_ids.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#uniprot_ids <- df$uniprot_id  # adjust column name as needed

# 2. Conectar con el biomart de Ensembl para Mus musculus
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 3. Mapear UniProt IDs a Entrez Gene IDs
conversion <- getBM(
  attributes = c("uniprotswissprot", "entrezgene_id"),
  filters = "uniprotswissprot",
  values = Up_ids,
  mart = ensembl
)

# Filtrar los valores válidos
entrez_ids <- unique(na.omit(conversion$entrezgene_id))

# 4. Enriquecimiento GO BP
UPregBP <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",           # Puede ser "BP", "MF" o "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 4. Enriquecimiento GO MF
UPregMF <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",           # Puede ser "BP", "MF" o "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 5. Visualización
DotPlot_BP_Up <- dotplot(UPregBP, showCategory = 30, font.size = 10) + ggtitle("GO Enrichment: Biological Process")
ggsave("GO_BP_enrichment.tiff", plot = DotPlot_BP_Up, dpi = 300, width = 8, height = 9, units = "in", device = "tiff")

DotPlot_MF_Up <- dotplot(UPregMF, showCategory = 30, font.size = 10) + ggtitle("GO Enrichment: Molecular Function")
ggsave("GO_MF_enrichment.tiff", plot = DotPlot_MF_Up, dpi = 300, width = 8, height = 9, units = "in", device = "tiff")


# Barplot
barplot_BP_Up <- barplot(UPregBP, showCategory = 20, font.size = 10) + ggtitle("GO Enrichment: Biological Function")
ggsave("GO_BP_enrichment.tiff", plot = barplot_BP_Up, dpi = 300, width = 8, height = 6, units = "in", device = "tiff")

# Save as a CSV
write.csv(as.data.frame(UPregBP), file="GO_enrichment_results_UPregBP_mouse.csv", row.names=FALSE)
write.csv(as.data.frame(UPregMF), file="GO_enrichment_results_UPregMF_mouse.csv", row.names=FALSE)


###OTRAS FIGURAS OPCIONALES###
# Network plot: genes and their enriched GO terms
#cnetplot(ego, showCategory=10, foldChange=NULL, circular=FALSE, colorEdge=TRUE)

# Enrichment map
#emapplot(ego, showCategory=15, layout = "kk")  # layout can be "kk", "fr", "circle", etc.

# Combined GO + gene plot
#goplot(ego)

# Heatmap of genes per GO category
#heatplot(ego, showCategory=10)

# UpSet plot for gene overlap
#upsetplot(ego)