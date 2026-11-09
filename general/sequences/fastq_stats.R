# =============================================================================
# Number of reads and read-length summary of a FASTQ file
#
# Usage: edit the settings below, then run  Rscript fastq_stats.R
# =============================================================================

if (!requireNamespace("ShortRead", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ShortRead", ask = FALSE, update = FALSE)
}
library(ShortRead)

# --- Settings ----------------------------------------------------------------
fastq_file <- "sample_R1.fastq.gz"

# --- Statistics --------------------------------------------------------------
fq <- readFastq(fastq_file)
read_lengths <- width(sread(fq))

cat("File:", fastq_file, "\n")
cat("Number of reads:", length(read_lengths), "\n")
cat("Read length:\n")
print(summary(read_lengths))
