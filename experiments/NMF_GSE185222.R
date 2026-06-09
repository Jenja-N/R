################################################################################
# NMF + TF ENRICHMENT ANALYSIS ON GSE185222
# Dental pulp: sound vs enamel caries vs deep caries
# 6,582 cells x 20,293 genes (log-normalized)
################################################################################

setwd("D:/!jenja/ARC/DATA3/")
library(Matrix)

# ============================================================
# STEP 0: LOAD DATA (reuse from previous analysis)
# ============================================================

cat("=== STEP 0: LOADING GSE185222 ===\n")

extract_dir <- "GSE185222_extracted"
csv_files <- list.files(extract_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

sample_map <- list(
  "GSM5608427" = list(condition = "sound", id = "DTP_01"),
  "GSM5608428" = list(condition = "enamel_caries", id = "DTP_02"),
  "GSM5608429" = list(condition = "deep_caries", id = "DTP_04"),
  "GSM5608430" = list(condition = "deep_caries", id = "DTP_05")
)

all_expr <- NULL
all_meta <- data.frame(
  cell = character(),
  sample = character(),
  condition = character(),
  stringsAsFactors = FALSE
)

for (f in csv_files) {
  con <- gzfile(f, "r")
  df <- read.csv(con, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  close(con)

  gsm <- gsub("_.*$", "", basename(f))
  info <- sample_map[[gsm]]

  expr <- as.matrix(df)

  if (is.null(all_expr)) {
    all_expr <- expr
  } else {
    common_genes <- intersect(rownames(all_expr), rownames(expr))
    all_expr <- all_expr[common_genes, , drop = FALSE]
    expr <- expr[common_genes, , drop = FALSE]
    all_expr <- cbind(all_expr, expr)
  }

  meta_add <- data.frame(
    cell = colnames(expr),
    sample = info$id,
    condition = info$condition,
    stringsAsFactors = FALSE
  )
  all_meta <- rbind(all_meta, meta_add)
}

cat("Loaded:", nrow(all_expr), "genes x", ncol(all_expr), "cells\n")
cat("Conditions:", paste(table(all_meta$condition), collapse = ", "), "\n")

# ============================================================
# STEP 1: SELECT VARIABLE GENES FOR NMF
# ============================================================

cat("\n=== STEP 1: VARIABLE GENE SELECTION ===\n")

# Calculate variance for each gene
# For log-normalized data, use row variance
row_var <- apply(all_expr, 1, var, na.rm = TRUE)
row_mean <- rowMeans(all_expr, na.rm = TRUE)

# Select top variable genes (top 3000 by variance)
top_n <- 3000
var_genes <- names(sort(row_var, decreasing = TRUE))[1:top_n]
cat("Selected top", top_n, "variable genes\n")
cat("Variance range:", round(min(row_var[var_genes]), 4), "to", round(max(row_var[var_genes]), 4), "\n")

# Subset expression matrix
expr_nmf <- all_expr[var_genes, , drop = FALSE]
cat("NMF input:", nrow(expr_nmf), "genes x", ncol(expr_nmf), "cells\n")

# Ensure non-negative (NMF requires non-negative input)
# Data is log-normalized, so shift to make all values >= 0
min_val <- min(expr_nmf, na.rm = TRUE)
if (min_val < 0) {
  expr_nmf <- expr_nmf - min_val + 0.01
  cat("Shifted data by", -min_val + 0.01, "to ensure non-negativity\n")
}

# ============================================================
# STEP 2: NMF WITH MULTIPLE K VALUES
# ============================================================

cat("\n=== STEP 2: NMF ANALYSIS ===\n")

# Try k = 5, 8, 10, 12, 15
k_values <- c(5, 8, 10, 12, 15)

# Simple NMF using multiplicative updates
# W: genes x k (basis matrix)
# H: k x cells (coefficient matrix)
# V ~ W * H

run_nmf <- function(V, k, max_iter = 200, tol = 1e-4) {
  n_genes <- nrow(V)
  n_cells <- ncol(V)

  # Initialize W and H with random non-negative values
  set.seed(42)
  W <- matrix(runif(n_genes * k, 0.01, 1), nrow = n_genes)
  H <- matrix(runif(k * n_cells, 0.01, 1), nrow = k)

  prev_err <- Inf

  for (iter in 1:max_iter) {
    # Update H
    WH <- W %*% H
    WH[WH == 0] <- 1e-10
    H <- H * (t(W) %*% (V / WH)) / (t(W) %*% matrix(1, nrow = n_genes, ncol = n_cells) + 1e-10)
    H[H < 0] <- 0

    # Update W
    WH <- W %*% H
    WH[WH == 0] <- 1e-10
    W <- W * ((V / WH) %*% t(H)) / (matrix(1, nrow = n_genes, ncol = n_cells) %*% t(H) + 1e-10)
    W[W < 0] <- 0

    # Calculate error
    WH <- W %*% H
    err <- sum((V - WH)^2)

    if (abs(prev_err - err) / prev_err < tol && iter > 50) {
      cat("    Converged at iteration", iter, "\n")
      break
    }
    prev_err <- err

    if (iter %% 50 == 0) {
      cat("    Iteration", iter, "error:", round(err, 2), "\n")
    }
  }

  # Calculate reconstruction error
  final_WH <- W %*% H
  recon_error <- sqrt(mean((V - final_WH)^2))

  list(W = W, H = H, error = recon_error, iterations = iter)
}

# Run NMF for each k
nmf_results <- list()

for (k in k_values) {
  cat("\n  Running NMF with k =", k, "\n")
  result <- run_nmf(expr_nmf, k)
  nmf_results[[as.character(k)]] <- result
  cat("  k =", k, "| Reconstruction RMSE:", round(result$error, 4), "\n")

  # Calculate module scores per cell
  H <- result$H
  module_scores <- t(H)  # cells x k

  # Check which module contains PTN/MDK/COL1A1
  W <- result$W
  rownames(W) <- var_genes

  target_genes <- c("PTN", "MDK", "COL1A1", "IGF1R", "SOX2", "IGFBP7", "NFKBIA", "IL6")

  cat("  Module assignments for target genes:\n")
  for (gene in target_genes) {
    if (gene %in% var_genes) {
      module <- which.max(W[gene, ])
      weight <- max(W[gene, ])
      cat("    ", gene, "-> Module", module, "(weight:", round(weight, 4), ")\n")
    } else {
      cat("    ", gene, "-> NOT in variable gene set\n")
    }
  }

  # Test module association with condition
  cat("  Module scores by condition:\n")
  for (mod in 1:k) {
    mod_name <- paste0("Module_", mod)
    scores <- module_scores[, mod]

    cat("\n    ", mod_name, ":\n", sep="")
    for (cond in c("sound", "enamel_caries", "deep_caries")) {
      idx <- which(all_meta$condition == cond)
      if (length(idx) > 0) {
        mean_score <- mean(scores[idx])
        cat("      ", cond, ": mean =", round(mean_score, 4), "\n", sep="")
      }
    }

    # ANOVA across conditions
    cond_vec <- all_meta$condition
    score_vec <- scores

    # Simple F-test approximation
    group_means <- tapply(score_vec, cond_vec, mean)
    group_vars <- tapply(score_vec, cond_vec, var)
    group_ns <- tapply(score_vec, cond_vec, length)

    overall_mean <- mean(score_vec)
    ss_between <- sum(group_ns * (group_means - overall_mean)^2)
    ss_within <- sum((group_ns - 1) * group_vars, na.rm = TRUE)

    df_between <- length(group_means) - 1
    df_within <- sum(group_ns) - length(group_means)

    if (df_within > 0 && ss_within > 0) {
      ms_between <- ss_between / df_between
      ms_within <- ss_within / df_within
      f_stat <- ms_between / ms_within
      p_val <- pf(f_stat, df_between, df_within, lower.tail = FALSE)

      cat("      ANOVA F =", round(f_stat, 3), ", p =", format(p_val, digits = 3), "\n")

      if (p_val < 0.05) {
        cat("      *** SIGNIFICANTLY ASSOCIATED WITH CONDITION ***\n")
      }
    }
  }
}

# ============================================================
# STEP 3: SELECT BEST K AND ANALYZE MODULES
# ============================================================

cat("\n=== STEP 3: MODULE ANALYSIS FOR BEST K ===\n")

# Select k=10 as balance between interpretability and complexity
best_k <- "10"
if (best_k %in% names(nmf_results)) {
  result <- nmf_results[[best_k]]
  W <- result$W
  H <- result$H
  rownames(W) <- var_genes

  cat("Analyzing k =", best_k, "\n")

  # Top genes per module
  cat("\n--- Top 20 genes per module ---\n")
  for (mod in 1:as.numeric(best_k)) {
    cat("\nModule", mod, "top genes:\n")
    gene_weights <- W[, mod]
    top_genes <- names(sort(gene_weights, decreasing = TRUE))[1:20]
    for (i in 1:20) {
      cat("  ", i, ". ", top_genes[i], " (", round(gene_weights[top_genes[i]], 4), ")\n", sep="")
    }
  }

  # Find module with PTN/MDK/COL1A1
  cat("\n--- PTN/MDK/COL1A1 MODULE IDENTIFICATION ---\n")

  ptn_module <- which.max(W["PTN", ])
  mdk_module <- which.max(W["MDK", ])
  col1a1_module <- which.max(W["COL1A1", ])

  cat("PTN max weight in Module", ptn_module, "\n")
  cat("MDK max weight in Module", mdk_module, "\n")
  cat("COL1A1 max weight in Module", col1a1_module, "\n")

  # Check if they co-occur in same module
  if (length(unique(c(ptn_module, mdk_module, col1a1_module))) == 1) {
    cat("\n*** PTN, MDK, COL1A1 ALL IN MODULE", ptn_module, "***\n")

    # Get all genes in this module with high weight
    module_genes <- var_genes[W[, ptn_module] > quantile(W[, ptn_module], 0.9)]
    cat("High-weight genes in PTN/MDK/COL1A1 module:", length(module_genes), "\n")

    # Check for transcription factors
    # Common TF markers
    tf_markers <- c("RUNX2", "SP7", "LEF1", "TCF7", "SOX9", "SOX5", "SOX6", 
                    "TWIST1", "TWIST2", "GLI1", "GLI2", "GLI3", "SMAD1", "SMAD5",
                    "STAT1", "STAT3", "REL", "RELA", "NFKB1", "JUN", "FOS",
                    "MYC", "MYCN", "EGR1", "EGR2", "KLF4", "KLF5", "SP1")

    tf_in_module <- intersect(tf_markers, module_genes)
    cat("\nTranscription factors in module:\n")
    for (tf in tf_in_module) {
      cat("  ", tf, " (weight:", round(W[tf, ptn_module], 4), ")\n")
    }

    if (length(tf_in_module) == 0) {
      cat("  No known TFs in high-weight subset. Checking all TFs in module...\n")
      # Broader TF search
      all_tfs <- c("RUNX", "SOX", "TCF", "LEF", "TWIST", "GLI", "SMAD", 
                   "STAT", "REL", "NFKB", "JUN", "FOS", "MYC", "EGR", 
                   "KLF", "SP", "IRF", "FOX", "GATA", "HOX", "PAX",
                   "MSX", "DLX", "BARX", "HAND", "MEIS", "PBX", "PKNOX")

      for (tf_prefix in all_tfs) {
        matching <- module_genes[grepl(tf_prefix, module_genes)]
        if (length(matching) > 0) {
          for (m in matching) {
            cat("  ", m, " (weight:", round(W[m, ptn_module], 4), ")\n")
          }
        }
      }
    }
  } else {
    cat("\nPTN/MDK/COL1A1 are in DIFFERENT modules\n")
    cat("  PTN -> Module", ptn_module, "\n")
    cat("  MDK -> Module", mdk_module, "\n")
    cat("  COL1A1 -> Module", col1a1_module, "\n")
  }
}

# ============================================================
# STEP 4: MODULE-CONDITION ASSOCIATION
# ============================================================

cat("\n=== STEP 4: MODULE-CONDITION ASSOCIATION (k=10) ===\n")

if (best_k %in% names(nmf_results)) {
  result <- nmf_results[[best_k]]
  H <- result$H
  module_scores <- t(H)

  # For each module, test sound vs deep caries
  for (mod in 1:as.numeric(best_k)) {
    scores <- module_scores[, mod]

    s_idx <- which(all_meta$condition == "sound")
    d_idx <- which(all_meta$condition == "deep_caries")

    if (length(s_idx) > 5 && length(d_idx) > 5) {
      tt <- t.test(scores[s_idx], scores[d_idx])

      cat("\nModule", mod, ":\n")
      cat("  Sound mean:", round(mean(scores[s_idx]), 4), "\n")
      cat("  Deep mean:", round(mean(scores[d_idx]), 4), "\n")
      cat("  Difference:", round(mean(scores[d_idx]) - mean(scores[s_idx]), 4), "\n")
      cat("  t-test p:", format(tt$p.value, digits = 3), "\n")

      if (tt$p.value < 0.05) {
        if (mean(scores[d_idx]) > mean(scores[s_idx])) {
          cat("  *** UP in deep caries ***\n")
        } else {
          cat("  *** DOWN in deep caries ***\n")
        }
      }
    }
  }
}

# ============================================================
# STEP 5: CORRELATION NETWORK AROUND PTN/MDK
# ============================================================

cat("\n=== STEP 5: PTN/MDK CORRELATION NETWORK ===\n")

# Find genes highly correlated with PTN and MDK
ptn_cor <- apply(all_expr, 1, function(x) cor(x, all_expr["PTN", ], method = "spearman", use = "complete.obs"))
mdk_cor <- apply(all_expr, 1, function(x) cor(x, all_expr["MDK", ], method = "spearman", use = "complete.obs"))

# Genes correlated with BOTH PTN and MDK (|rho| > 0.3)
both_cor <- data.frame(
  gene = names(ptn_cor),
  ptn_rho = ptn_cor,
  mdk_rho = mdk_cor,
  stringsAsFactors = FALSE
)
both_cor$mean_rho <- (abs(both_cor$ptn_rho) + abs(both_cor$mdk_rho)) / 2
both_cor <- both_cor[order(both_cor$mean_rho, decreasing = TRUE), ]

# Top correlated genes
cat("\nTop 30 genes correlated with BOTH PTN and MDK:\n")
print(head(both_cor, 30))

# Check for TFs in top correlated genes
tf_list <- c("RUNX2", "SP7", "LEF1", "TCF7", "SOX9", "SOX5", "SOX6", 
             "TWIST1", "TWIST2", "GLI1", "GLI2", "GLI3", "SMAD1", "SMAD5",
             "STAT1", "STAT3", "REL", "RELA", "NFKB1", "JUN", "FOS",
             "MYC", "MYCN", "EGR1", "EGR2", "KLF4", "KLF5", "SP1",
             "MSX1", "MSX2", "DLX5", "DLX6", "BARX1", "BARX2", "HAND1", "HAND2")

tf_cor <- both_cor[both_cor$gene %in% tf_list, ]
cat("\nTFs correlated with PTN/MDK:\n")
print(tf_cor)

# ============================================================
# STEP 6: SUMMARY
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("NMF ANALYSIS COMPLETE\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("\nKey outputs:\n")
cat("1. NMF modules for k = 5, 8, 10, 12, 15\n")
cat("2. PTN/MDK/COL1A1 module identification\n")
cat("3. Module-condition association (ANOVA + t-tests)\n")
cat("4. PTN/MDK correlation network\n")
cat("5. TF candidates in regeneration module\n")
