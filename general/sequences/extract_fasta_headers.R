# =============================================================================
# Extract sequence identifiers from a FASTA file
#
# Keeps the first word of each header (the text after ">" up to the first
# space) and, optionally, only those starting with a given prefix.
#
# Usage: edit the settings below, then run  Rscript extract_fasta_headers.R
# =============================================================================

# --- Settings ----------------------------------------------------------------
fasta_file  <- "sequences.fasta"
prefix      <- ""                  # e.g. "T265_"; leave empty to keep all
output_file <- "fasta_headers.txt"

# --- Extract -----------------------------------------------------------------
lines   <- readLines(fasta_file, warn = FALSE)
headers <- sub("^>", "", grep("^>", lines, value = TRUE))
if (nzchar(prefix)) headers <- headers[startsWith(headers, prefix)]
ids <- sub(" .*", "", headers)

writeLines(ids, output_file)
message(length(ids), " identifiers saved to ", output_file)
