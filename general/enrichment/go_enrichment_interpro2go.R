# =============================================================================
# GO enrichment for non-model organisms using InterPro2GO
#
# GO terms are assigned to each protein through its InterPro domains, taken
# from a WormBase ParaSite GFF3 annotation, and the official InterPro2GO
# mapping. Enrichment is then run with clusterProfiler::enricher, using all
# annotated proteins as the background.
#
# Inputs:
#   - GFF3 annotation file, e.g.
#     opisthorchis_viverrini.PRJNA222628.WBPS19.annotations.gff3
#     (https://parasite.wormbase.org)
#   - interpro2go mapping file
#     (https://current.geneontology.org/ontology/external2go/interpro2go)
#   - plain-text file with one protein/transcript ID per line
#
# Usage: edit the settings below, then run  Rscript go_enrichment_interpro2go.R
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran_pkgs <- c("dplyr", "stringr", "tidyr", "ggplot2")
bioc_pkgs <- c("clusterProfiler", "enrichplot", "GO.db")
for (pkg in cran_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
missing <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) BiocManager::install(missing, ask = FALSE, update = FALSE)
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# --- Settings ----------------------------------------------------------------
gff_file         <- "opisthorchis_viverrini.PRJNA222628.WBPS19.annotations.gff3"
interpro2go_file <- "interpro2go.txt"
proteins_file    <- "proteins_of_interest.txt"
output_prefix    <- "GO_InterPro2GO"
plot_label       <- ""          # added to plot titles, e.g. "sEVs mRNA"

# --- 1. Protein to InterPro mapping from the GFF3 ----------------------------
gff <- read.delim(gff_file, comment.char = "#", header = FALSE, quote = "")
colnames(gff) <- c("seqid", "source", "type", "start", "end",
                   "score", "strand", "phase", "attributes")

gene2ipr <- gff %>%
  filter(str_detect(attributes, "InterPro")) %>%
  mutate(ID  = str_extract(attributes, "(?<=Parent=)[^;]+"),
         ID  = str_replace(ID, "^(gene:|mRNA:|transcript:)", ""),
         ID  = str_replace(ID, "\\.\\d+$", ""),
         IPR = str_extract(attributes, "IPR[0-9]+")) %>%
  filter(!is.na(IPR)) %>%
  distinct(ID, IPR)
message("Found ", nrow(gene2ipr), " protein-to-InterPro mappings")

# --- 2. InterPro to GO mapping -----------------------------------------------
interpro2go <- read.delim(interpro2go_file, header = FALSE,
                          comment.char = "!", quote = "", col.names = "raw") %>%
  mutate(IPR = str_extract(raw, "IPR[0-9]+"),
         GO  = str_extract_all(raw, "GO:[0-9]+")) %>%
  select(IPR, GO) %>%
  unnest(GO) %>%
  filter(!is.na(IPR), !is.na(GO)) %>%
  distinct()

gene2go <- inner_join(gene2ipr, interpro2go, by = "IPR",
                      relationship = "many-to-many") %>%
  distinct(ID, GO)
message("Generated ", nrow(gene2go), " protein-to-GO mappings")

# --- 3. Enrichment -----------------------------------------------------------
prot_list <- gsub("\\.\\d+$", "", trimws(readLines(proteins_file, warn = FALSE)))
prot_list <- unique(prot_list[prot_list != ""])

ego <- enricher(gene = prot_list,
                TERM2GENE = gene2go[, c("GO", "ID")],
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
if (is.null(ego)) stop("None of the proteins of interest have GO annotations.")

ego_res <- ego@result %>%
  mutate(Description = AnnotationDbi::mapIds(GO.db, keys = ID, column = "TERM",
                                             keytype = "GOID", multiVals = "first"),
         Ontology    = AnnotationDbi::mapIds(GO.db, keys = ID, column = "ONTOLOGY",
                                             keytype = "GOID", multiVals = "first"))
ego_res$Description[is.na(ego_res$Description)] <- ego_res$ID[is.na(ego_res$Description)]
ego@result <- ego_res
write.csv(ego_res, paste0(output_prefix, "_all.csv"), row.names = FALSE)

# --- 4. Results and dot plots per ontology -----------------------------------
theme_custom <- theme(axis.text    = element_text(size = 30),
                      axis.title   = element_text(size = 22),
                      legend.text  = element_text(size = 18),
                      legend.title = element_text(size = 19),
                      plot.title   = element_text(size = 22, face = "bold"))

ontologies <- c(BP = "Biological Process", MF = "Molecular Function", CC = "Cellular Component")
for (ont in names(ontologies)) {
  sub <- ego
  sub@result <- filter(ego_res, Ontology == ont)
  write.csv(sub@result, paste0(output_prefix, "_", ont, ".csv"), row.names = FALSE)
  if (sum(sub@result$p.adjust < 0.05) > 0) {
    p <- dotplot(sub, showCategory = 20, font.size = 16) +
      ggtitle(trimws(paste(ontologies[[ont]], plot_label))) +
      theme_custom
    ggsave(paste0(output_prefix, "_", ont, "_dotplot.tiff"), p, dpi = 300,
           width = 10, height = 6, units = "in", compression = "lzw")
  }
}

message("GO enrichment completed")
