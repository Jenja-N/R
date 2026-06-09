################################################################################
# FINAL CROSS-VERIFICATION — GSE185222 CSV + COMPLETE SUMMARY
################################################################################

setwd("D:/!jenja/ARC/DATA3/")
library(Matrix)

################################################################################
# PART 1: GSE185222 — CSV FORMAT (NOT 10X)
################################################################################

cat("\n========================================\n")
cat("GSE185222: CSV ANALYSIS\n")
cat("========================================\n")

extract_dir <- "GSE185222_extracted"
csv_files <- list.files(extract_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
cat("CSV files:", length(csv_files), "\n")
print(csv_files)

# GSE185222 samples:
# GSM5608427_DTP_01 = sound (control)
# GSM5608428_DTP_02 = enamel caries  
# GSM5608429_DTP_04 = deep caries
# GSM5608430_DTP_05 = deep caries

sample_map <- list(
  "GSM5608427" = "sound",
  "GSM5608428" = "enamel_caries", 
  "GSM5608429" = "deep_caries",
  "GSM5608430" = "deep_caries"
)

all_expr <- NULL
all_cond <- c()
all_cells <- c()

for (f in csv_files) {
  cat("\nLoading:", basename(f), "\n")

  # Read CSV - check format first
  con <- gzfile(f, "r")
  first_line <- readLines(con, n = 1)
  close(con)

  cat("First line:", substr(first_line, 1, 100), "\n")

  # Read full file
  con <- gzfile(f, "r")
  df <- read.csv(con, stringsAsFactors = FALSE, row.names = 1)
  close(con)

  cat("Dim:", nrow(df), "genes x", ncol(df), "cells\n")
  cat("First genes:", paste(rownames(df)[1:5], collapse = ", "), "\n")

  # Determine condition
  gsm <- gsub("_.*$", "", basename(f))
  cond <- sample_map[[gsm]]
  if (is.null(cond)) cond <- gsm

  cat("Condition:", cond, "\n")

  # Check if data is counts or log-normalized
  max_val <- max(as.matrix(df), na.rm = TRUE)
  cat("Max value:", max_val, "\n")

  if (max_val > 100) {
    cat("Data appears to be COUNTS (max > 100)\n")
    # Log-normalize
    lib_size <- colSums(df)
    sf <- median(lib_size[lib_size > 0])
    expr <- as.matrix(df)
    for (i in 1:ncol(expr)) {
      if (lib_size[i] > 0) {
        expr[,i] <- log1p(expr[,i] / lib_size[i] * sf)
      }
    }
  } else {
    cat("Data appears to already be normalized (max <= 100)\n")
    expr <- as.matrix(df)
  }

  # Merge
  if (is.null(all_expr)) {
    all_expr <- expr
  } else {
    common_genes <- intersect(rownames(all_expr), rownames(expr))
    all_expr <- all_expr[common_genes,]
    expr <- expr[common_genes,]
    all_expr <- cbind(all_expr, expr)
  }

  all_cond <- c(all_cond, rep(cond, ncol(expr)))
  all_cells <- c(all_cells, colnames(expr))
}

cat("\n=== MERGED DATA ===\n")
cat("Genes:", nrow(all_expr), "| Cells:", ncol(all_expr), "\n")
cat("Conditions:\n")
print(table(all_cond))

# IGF1R analysis
if ("IGF1R" %in% rownames(all_expr)) {
  cat("\n--- IGF1R by condition ---\n")
  for (cond in c("sound", "enamel_caries", "deep_caries")) {
    idx <- which(all_cond == cond)
    if (length(idx) > 0) {
      vals <- all_expr["IGF1R", idx]
      cat(cond, ": mean =", round(mean(vals, na.rm = TRUE), 4), 
          "| % pos =", round(mean(vals > 0, na.rm = TRUE)*100, 1), 
          "| n =", length(vals), "\n")
    }
  }

  # T-test sound vs deep
  s_idx <- which(all_cond == "sound")
  d_idx <- which(all_cond == "deep_caries")
  if (length(s_idx) > 0 && length(d_idx) > 0) {
    tt <- t.test(all_expr["IGF1R", s_idx], all_expr["IGF1R", d_idx])
    cat("\n--- T-test: sound vs deep caries ---\n")
    cat("Sound mean:", round(mean(all_expr["IGF1R", s_idx], na.rm = TRUE), 4), "\n")
    cat("Deep mean:", round(mean(all_expr["IGF1R", d_idx], na.rm = TRUE), 4), "\n")
    cat("p-value:", format(tt$p.value, digits = 3), "\n")

    if (tt$p.value < 0.05) {
      if (mean(all_expr["IGF1R", d_idx], na.rm = TRUE) > mean(all_expr["IGF1R", s_idx], na.rm = TRUE)) {
        cat("\n>>> IGF1R SIGNIFICANTLY HIGHER in deep caries <<<\n")
      } else {
        cat("\n>>> IGF1R SIGNIFICANTLY LOWER in deep caries <<<\n")
      }
    } else {
      cat("\n>>> No significant difference <<<\n")
    }
  }
}

# SOX2 analysis
if ("SOX2" %in% rownames(all_expr)) {
  cat("\n--- SOX2 by condition ---\n")
  for (cond in c("sound", "enamel_caries", "deep_caries")) {
    idx <- which(all_cond == cond)
    if (length(idx) > 0) {
      vals <- all_expr["SOX2", idx]
      cat(cond, ": mean =", round(mean(vals, na.rm = TRUE), 4), 
          "| % pos =", round(mean(vals > 0, na.rm = TRUE)*100, 1), "\n")
    }
  }
}

# IGFBP7 analysis
if ("IGFBP7" %in% rownames(all_expr)) {
  cat("\n--- IGFBP7 by condition ---\n")
  for (cond in c("sound", "enamel_caries", "deep_caries")) {
    idx <- which(all_cond == cond)
    if (length(idx) > 0) {
      vals <- all_expr["IGFBP7", idx]
      cat(cond, ": mean =", round(mean(vals, na.rm = TRUE), 4), 
          "| % pos =", round(mean(vals > 0, na.rm = TRUE)*100, 1), "\n")
    }
  }
}


################################################################################
# PART 2: COMPLETE CROSS-VERIFICATION SUMMARY
################################################################################

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║           CROSS-VERIFICATION: FINAL RESULTS                                  ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                              ║\n")
cat("║  DATASET 1: GSE184749 (Dental Mesenchyme, human tooth germ)                  ║\n")
cat("║  ─────────────────────────────────────────────────────────                   ║\n")
cat("║  Clusters: DEM (progenitor), POB, SOB, SOBP, OB                              ║\n")
cat("║  IGF1R highest in: DEM (0.117) — PROGENITOR cluster                         ║\n")
cat("║  SOX2 highest in: POB (0.59% cells) — PRE-ODONTOBLAST                       ║\n")
cat("║  GLOBAL IGF1R+ SOX2+ co-expression: 4 / 17,510 cells (0.023%)                ║\n")
cat("║                                                                              ║\n")
cat("║  >>> VERDICT: IGF1R and SOX2 are NOT co-expressed in dental mesenchyme <<<   ║\n")
cat("║                                                                              ║\n")
cat("║  They are in DIFFERENT clusters:                                               ║\n")
cat("║    - IGF1R = DEM (common progenitor, early stage)                            ║\n")
cat("║    - SOX2  = POB (pre-odontoblast, later stage)                              ║\n")
cat("║                                                                              ║\n")
cat("║  This means: IGF1R is on EARLY progenitors, SOX2 on LATE progenitors         ║\n")
cat("║  They represent DIFFERENT stages of odontoblast differentiation                ║\n")
cat("║                                                                              ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
cat("║  DATASET 2: GSE185222 (Dental pulp sound vs caries)                          ║\n")
cat("║  ─────────────────────────────────────────────────                             ║\n")

# Print GSE185222 results
if ("IGF1R" %in% rownames(all_expr)) {
  s_mean <- mean(all_expr["IGF1R", all_cond == "sound"], na.rm = TRUE)
  d_mean <- mean(all_expr["IGF1R", all_cond == "deep_caries"], na.rm = TRUE)
  cat("║  IGF1R sound:", sprintf("%-10.4f", s_mean), "| deep caries:", sprintf("%-10.4f", d_mean), "║\n")

  s_idx <- which(all_cond == "sound")
  d_idx <- which(all_cond == "deep_caries")
  if (length(s_idx) > 0 && length(d_idx) > 0) {
    tt <- t.test(all_expr["IGF1R", s_idx], all_expr["IGF1R", d_idx])
    if (tt$p.value < 0.05 && d_mean > s_mean) {
      cat("║  >>> IGF1R HIGHER in deep caries (p < 0.05) <<<                            ║\n")
    } else if (tt$p.value < 0.05) {
      cat("║  >>> IGF1R LOWER in deep caries (p < 0.05) <<<                             ║\n")
    } else {
      cat("║  >>> No significant difference in IGF1R between sound and deep caries <<<  ║\n")
    }
  }
}

cat("║                                                                              ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
cat("║  DATASET 3: GSE167251 (Aging pulp) — from published data                     ║\n")
cat("║  ─────────────────────────────────────────────────                             ║\n")
cat("║  IGF signaling INCREASES with age via IGFBP7 compensation                    ║\n")
cat("║  IGF1R does NOT decrease — it is maintained/upregulated                       ║\n")
cat("║                                                                              ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
cat("║  DATASET 4: GSE146123 (Dental Cell Atlas) — from published data              ║\n")
cat("║  ─────────────────────────────────────────────────                             ║\n")
cat("║  IGF1R + Sox2 co-expression CONFIRMED in mouse incisor & human apical papilla║\n")
cat("║  BUT: This is a DIFFERENT tissue (apical papilla vs tooth germ pulp)         ║\n")
cat("║                                                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")


cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                    FALSIFICATION STATUS                                        ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                              ║\n")
cat("║  HYPOTHESIS: "IGF1R is co-expressed with Sox2 on dental stem cells"          ║\n")
cat("║                                                                              ║\n")
cat("║  GSE146123 (apical papilla):     ✓ CONFIRMED                                 ║\n")
cat("║  GSE184749 (tooth germ DM):      ✗ FALSIFIED — only 0.023% co-expression    ║\n")
cat("║  GSE184749: IGF1R and SOX2 are on DIFFERENT differentiation stages            ║\n")
cat("║                                                                              ║\n")
cat("║  CONCLUSION: The hypothesis is PARTIALLY FALSIFIED                           ║\n")
cat("║                                                                              ║\n")
cat("║  IGF1R is on EARLY progenitors (DEM), SOX2 on LATE progenitors (POB)         ║\n")
cat("║  They are NOT co-expressed in the same cells in dental mesenchyme             ║\n")
cat("║                                                                              ║\n")
cat("║  REVISED MODEL:                                                               ║\n")
cat("║    DEM (IGF1R+) → POB (SOX2+) → OB (mature odontoblast)                      ║\n")
cat("║    IGF1R drives EARLY progenitor maintenance                                  ║\n")
cat("║    SOX2 marks LATE progenitor/transition state                                 ║\n")
cat("║    IGF1R is DOWNREGULATED before SOX2 is UPREGULATED                         ║\n")
cat("║                                                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
