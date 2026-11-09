# =============================================================================
# Latent class analysis (LCA) of binary diagnostic test results
#
# Fits one- and two-class models with randomLCA, in three forms:
#   - standard latent class models
#   - models with a random effect (conditional dependence between tests)
#   - random-effect models whose loadings vary by test (constload = FALSE)
# and compares them by BIC (lowest is best). Posterior class probabilities of
# the chosen model are exported to Excel.
#
# Input: Excel file with one row per individual and one column per test,
# coded 0 (negative) / 1 (positive). Other columns must be removed.
#
# Usage: edit the settings below, then run  Rscript random_lca.R
# =============================================================================

for (pkg in c("randomLCA", "readxl", "writexl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(randomLCA)
library(readxl)
library(writexl)

# --- Settings ----------------------------------------------------------------
input_file   <- "LCA_data.xlsx"
output_file  <- "postClassProbs.xlsx"
export_model <- "lca2"        # model whose posterior probabilities are exported
seed         <- 1234          # fixed seed so that fits are reproducible

# --- Data --------------------------------------------------------------------
mydata <- as.data.frame(read_excel(input_file))
if (!all(unlist(mydata) %in% c(0, 1, NA))) stop("All columns must be coded 0/1.")

# --- Models ------------------------------------------------------------------
fit <- function(...) {
  set.seed(seed)
  randomLCA(mydata, cores = 1, ...)
}

models <- list(
  lca1         = fit(nclass = 1),
  lca2         = fit(nclass = 2),
  lca1_random  = fit(nclass = 1, random = TRUE),
  lca2_random  = fit(nclass = 2, random = TRUE),
  lca1_random2 = fit(nclass = 1, random = TRUE, constload = FALSE, quadpoints = 61),
  lca2_random2 = fit(nclass = 2, random = TRUE, constload = FALSE)
)

bic <- data.frame(model = names(models), BIC = sapply(models, BIC), row.names = NULL)
print(bic[order(bic$BIC), ], row.names = FALSE)
message("Lowest BIC: ", bic$model[which.min(bic$BIC)])

# --- Posterior class probabilities -------------------------------------------
probs <- postClassProbs(models[[export_model]])
write_xlsx(as.data.frame(probs), output_file)
message("Posterior probabilities of ", export_model, " saved to ", output_file)
