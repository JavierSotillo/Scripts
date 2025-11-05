if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

setwd("C:/Users/javier.sotillo/Desktop")

library(ShortRead)

# Load a FASTQ file
fq <- readFastq("MV1_R1.fastq.gz")

# Get read length
read_lengths <- width(sread(fq))
summary(read_lengths)

# Number of reads
length(read_lengths)
