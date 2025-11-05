# GO enrichment for Opisthorchis viverrini using InterPro2GO mapping

# ------------------------
# 0. Set working directory
# ------------------------
setwd("C:/Users/javier.sotillo/OneDrive - Instituto de Salud Carlos III/Ov_EV_omics/RNA_Seq/RNAseq_Ov_EVs_Matt/mRNA_analysis")

# ------------------------
# 1. Install and load packages
# ------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("dplyr", "stringr", "tidyr", "clusterProfiler", "enrichplot", "ggplot2", "readr")
lapply(pkgs, function(x) {
  if (!require(x, character.only = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
})

if (!requireNamespace("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")
library(GO.db)


# ------------------------
# 2. Read GFF3 annotation file
# ------------------------
gff <- read.delim("opisthorchis_viverrini.PRJNA222628.WBPS19.annotations.gff3",
                  comment.char = "#", header = FALSE, quote = "")

colnames(gff) <- c("seqid", "source", "type", "start", "end",
                   "score", "strand", "phase", "attributes")

# ------------------------
# 3. Extract protein ID and InterPro domain (IPR ID)
# ------------------------
gene2ipr <- gff %>%
  dplyr::filter(stringr::str_detect(attributes, "InterPro")) %>%
  dplyr::mutate(
    ID = stringr::str_extract(attributes, "(?<=Parent=)[^;]+"),
    ID = stringr::str_replace(ID, "^(gene:|mRNA:|transcript:)", ""),
    ID = stringr::str_replace(ID, "\\.\\d+$", ""),
    IPR = stringr::str_extract(attributes, "IPR[0-9]+")
  ) %>%
  dplyr::filter(!is.na(IPR)) %>%
  dplyr::select(ID, IPR) %>%
  dplyr::distinct()

cat("Found", nrow(gene2ipr), "protein-to-InterPro mappings\n")

# ------------------------
# 4. Load InterPro2GO mapping file
# ------------------------
interpro2go <- read.delim("interpro2go.txt",
                          header = FALSE, comment.char = "!", quote = "")
colnames(interpro2go) <- c("raw")

interpro2go <- interpro2go %>%
  dplyr::mutate(
    IPR = stringr::str_extract(raw, "IPR[0-9]+"),
    GO = stringr::str_extract_all(raw, "GO:[0-9]+")
  ) %>%
  dplyr::select(IPR, GO) %>%
  tidyr::unnest(GO) %>%
  dplyr::filter(!is.na(IPR), !is.na(GO)) %>%
  dplyr::distinct()

# ------------------------
# 5. Join InterPro and GO mappings
# ------------------------
gene2go <- dplyr::inner_join(gene2ipr, interpro2go, by = "IPR", relationship = "many-to-many") %>%
  dplyr::distinct(ID, GO)

cat("Generated", nrow(gene2go), "protein-to-GO mappings\n")

# ------------------------
# 6. Read list of proteins of interest
# ------------------------
prot_list <- readLines("sEVs.txt")
prot_list <- gsub("\\.\\d+$", "", prot_list)

# ------------------------
# 7. Perform GO enrichment analysis
# ------------------------
ego <- clusterProfiler::enricher(
  gene = prot_list,
  TERM2GENE = gene2go[, c("GO", "ID")],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Add GO term names and ontology type
ego_res <- ego@result %>%
  dplyr::mutate(
    Description = mapIds(GO.db, keys = ID, column = "TERM", keytype = "GOID", multiVals = "first"),
    Ontology = mapIds(GO.db, keys = ID, column = "ONTOLOGY", keytype = "GOID", multiVals = "first")
  )

ego_res$Description[is.na(ego_res$Description)] <- ego_res$ID[is.na(ego_res$Description)]
ego@result <- ego_res


# ------------------------
# 8. Split into Biological Process / Molecular Function / Cellular Component
# ------------------------
ego_BP <- ego; ego_BP@result <- ego_res %>% dplyr::filter(Ontology == "BP")
ego_MF <- ego; ego_MF@result <- ego_res %>% dplyr::filter(Ontology == "MF")
ego_CC <- ego; ego_CC@result <- ego_res %>% dplyr::filter(Ontology == "CC")


# ------------------------
# 9. Save results + Plots
# ------------------------
write.csv(ego_res, "sEVs_mRNA_GO_enrichment_from_InterPro2GO_FULL.csv", row.names = FALSE)
write.csv(ego_BP@result, "sEVs_mRNA_GO_BP.csv", row.names = FALSE)
write.csv(ego_MF@result, "sEVs_mRNA_GO_MF.csv", row.names = FALSE)
write.csv(ego_CC@result, "sEVs_mRNA_GO_CC.csv", row.names = FALSE)

# Define a consistent theme for all dotplots
theme_custom <- theme(
  axis.text = element_text(size = 30),
  axis.title = element_text(size = 22),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 19),
  plot.title = element_text(size = 22, face = "bold")
)

# Biological Process
p1 <- dotplot(ego_BP, showCategory = 20, font.size = 16) +
  ggtitle("Biological Process sEVs mRNA") +
  theme_custom
ggsave("sEVs_mRNA_GO_BP_dotplot.tiff", p1, dpi = 300, width = 10, height = 6, units = "in")

# Molecular Function
p2 <- dotplot(ego_MF, showCategory = 20,  font.size = 16) +
  ggtitle("Molecular Function sEVs mRNA") +
  theme_custom
ggsave("sEVs_mRNA_GO_MF_dotplot.tiff", p2, dpi = 300, width = 10, height = 6, units = "in")

# Cellular Component
p3 <- dotplot(ego_CC, showCategory = 20,  font.size = 16) +
  ggtitle("Cellular Component sEVs mRNA") +
  theme_custom
ggsave("sEVs_mRNA_GO_CC_dotplot.tiff", p3, dpi = 300, width = 10, height = 6, units = "in")

cat("\n✅ GO enrichment completed and separated into BP / MF / CC.\n")
