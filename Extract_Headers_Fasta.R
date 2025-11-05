# R script to extract headers from a FASTA file

#Set working Directory
setwd("C:/Users/javier.sotillo/Desktop")

# Path to your FASTA file
fasta_file <- "OvEVs_mRNA_seqs.fasta"

# Read all lines
lines <- readLines(fasta_file)

# Extract header lines (those starting with ">")
headers <- grep("^>", lines, value = TRUE)

# Remove ">" and keep only those starting with "T265_"
headers_clean <- sub("^>", "", headers)
headers_filtered <- grep("^T265_", headers_clean, value = TRUE)

# Remove everything after the first space
headers_final <- sub(" .*", "", headers_filtered)

# Print headers to console
print(headers_final)

# Save headers to a text file (optional)
writeLines(headers_final, "fasta_headers.txt")

cat("Headers have been extracted and saved to 'fasta_headers.txt'\n")
