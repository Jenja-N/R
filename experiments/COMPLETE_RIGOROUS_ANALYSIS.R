################################################################################
# COMPLETE CROSS-VERIFICATION ANALYSIS
# Strict scientific rigor, absolute honesty
# All datasets, all target genes, all statistical tests
################################################################################

setwd("D:/!jenja/ARC/DATA3/")
library(Matrix)

# ============================================================
# SECTION 0: DATA INVENTORY AND QC
# ============================================================

cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 0: DATA INVENTORY AND QUALITY CONTROL\n")
cat("=" ,rep("=", 78), "\n", sep="")

# GSE184749 files
cat("\n--- GSE184749 Files ---\n")
gse184749_files <- list.files("GSE184749", full.names = TRUE)
print(gse184749_files)

# GSE185222 files  
cat("\n--- GSE185222 Files ---\n")
gse185222_files <- list.files("GSE185222_RAW", full.names = TRUE)
print(gse185222_files)

# GSE167251 files
cat("\n--- GSE167251 Files ---\n")
gse167251_files <- list.files("GSE167251", full.names = TRUE)
print(gse167251_files)

# GSE132472 files
cat("\n--- GSE132472 Files ---\n")
gse132472_files <- list.files(pattern = "GSE132472", full.names = TRUE)
print(gse132472_files)


# ============================================================
# SECTION 1: GSE184749 COMPLETE ANALYSIS
# Dental Mesenchyme, Human Tooth Germ, sci-RNA-seq
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 1: GSE184749 - DENTAL MESENCHYME\n")
cat("=" ,rep("=", 78), "\n", sep="")

# 1.1 Load data
counts_file <- "GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"
cell_file <- "GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"
gene_file <- "GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"

counts <- readMM(counts_file)
counts <- as(counts, "dgCMatrix")

cells <- read.csv(gzfile(cell_file), stringsAsFactors = FALSE)
genes <- read.csv(gzfile(gene_file), stringsAsFactors = FALSE)

colnames(counts) <- cells[[1]]
gene_col <- grep("gene|symbol|name", colnames(genes), ignore.case = TRUE, value = TRUE)[1]
if (is.na(gene_col)) gene_col <- colnames(genes)[1]
rownames(counts) <- genes[[gene_col]]

# 1.2 Basic QC
cat("\n--- 1.2 BASIC QC ---\n")
cat("Matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")
cat("Non-zero entries:", length(counts@x), "(", round(length(counts@x)/(nrow(counts)*ncol(counts))*100, 4), "% density)\n")

lib_size <- colSums(counts)
cat("Library size - min:", min(lib_size), "median:", median(lib_size), "max:", max(lib_size), "\n")
cat("Cells with < 500 UMIs:", sum(lib_size < 500), "\n")
cat("Cells with > 10000 UMIs:", sum(lib_size > 10000), "\n")

genes_per_cell <- colSums(counts > 0)
cat("Genes per cell - min:", min(genes_per_cell), "median:", median(genes_per_cell), "max:", max(genes_per_cell), "\n")

# 1.3 Annotation structure
cat("\n--- 1.3 ANNOTATION STRUCTURE ---\n")
cat("Cell annotation columns:", paste(colnames(cells), collapse = ", "), "\n")
for (col in colnames(cells)) {
  vals <- cells[[col]]
  if (is.character(vals) || is.factor(vals)) {
    uvals <- unique(vals)
    if (length(uvals) <= 20) {
      cat("\n", col, "(", length(uvals), "unique values):\n")
      print(table(vals, useNA = "ifany"))
    }
  } else {
    cat("\n", col, "(numeric): range =", min(vals, na.rm = TRUE), "to", max(vals, na.rm = TRUE), "\n")
  }
}

# 1.4 Target gene detection
cat("\n--- 1.4 TARGET GENE DETECTION ---\n")
target_genes <- c("IGF1R", "SOX2", "IGF1", "IGFBP7", "IGFBP3", "IGFBP5", "PTN", "MDK", "FGFR1", "FGFR2", "BMPR1A", "BMPR1B", "BMP2", "BMP4", "WNT5A", "WNT10A", "AXIN2", "LEF1", "RUNX2", "SP7", "DSPP", "DMP1", "COL1A1", "THY1", "NT5E", "ENG", "MCAM", "CD146")
avail_genes <- rownames(counts)

detection <- data.frame(
  gene = target_genes,
  present = target_genes %in% avail_genes,
  stringsAsFactors = FALSE
)
print(detection)

present_targets <- detection$gene[detection$present]
cat("\nPresent target genes:", length(present_targets), "\n")
cat("Missing target genes:", paste(detection$gene[!detection$present], collapse = ", "), "\n")

# 1.5 Normalization
sf <- median(lib_size[lib_size > 0])
norm_counts <- counts
norm_counts@x <- norm_counts@x / rep(lib_size, diff(norm_counts@p)) * sf
log_expr <- log1p(norm_counts)

# 1.6 Cluster-wise expression analysis
cat("\n--- 1.6 CLUSTER-WISE EXPRESSION ANALYSIS ---\n")
clusters <- cells$cluster
cluster_names <- sort(unique(clusters))

# Expression matrix for present targets
expr_matrix <- as.matrix(log_expr[present_targets[present_targets %in% avail_genes], , drop = FALSE])

# Mean expression by cluster
mean_expr <- sapply(cluster_names, function(cl) {
  idx <- which(clusters == cl)
  if (length(idx) > 0) rowMeans(expr_matrix[, idx, drop = FALSE])
  else rep(NA, nrow(expr_matrix))
})
rownames(mean_expr) <- rownames(expr_matrix)

cat("\nMean log-normalized expression by cluster:\n")
print(round(mean_expr, 4))

# Percentage of expressing cells by cluster
pct_expr <- sapply(cluster_names, function(cl) {
  idx <- which(clusters == cl)
  if (length(idx) > 0) rowMeans(expr_matrix[, idx, drop = FALSE] > 0) * 100
  else rep(NA, nrow(expr_matrix))
})
rownames(pct_expr) <- rownames(expr_matrix)

cat("\nPercentage of cells expressing gene (>0) by cluster:\n")
print(round(pct_expr, 2))

# 1.7 IGF1R and SOX2 detailed analysis
cat("\n--- 1.7 IGF1R AND SOX2 DETAILED ANALYSIS ---\n")

if (all(c("IGF1R", "SOX2") %in% present_targets)) {

  # Per-cluster statistics
  cat("\nPer-cluster IGF1R statistics:\n")
  for (cl in cluster_names) {
    idx <- which(clusters == cl)
    igf1r_vals <- expr_matrix["IGF1R", idx]
    sox2_vals <- expr_matrix["SOX2", idx]

    cat("\n", cl, "(n=", length(idx), "):\n", sep="")
    cat("  IGF1R: mean=", round(mean(igf1r_vals), 4), 
        " median=", round(median(igf1r_vals), 4),
        " sd=", round(sd(igf1r_vals), 4),
        " %pos=", round(mean(igf1r_vals > 0)*100, 2), "\n", sep="")
    cat("  SOX2:  mean=", round(mean(sox2_vals), 4), 
        " median=", round(median(sox2_vals), 4),
        " sd=", round(sd(sox2_vals), 4),
        " %pos=", round(mean(sox2_vals > 0)*100, 2), "\n", sep="")

    # Co-expression
    both_pos <- sum(igf1r_vals > 0 & sox2_vals > 0)
    igf1r_only <- sum(igf1r_vals > 0 & sox2_vals == 0)
    sox2_only <- sum(igf1r_vals == 0 & sox2_vals > 0)
    neither <- sum(igf1r_vals == 0 & sox2_vals == 0)

    cat("  Co-expression:\n")
    cat("    IGF1R+ SOX2+ :", both_pos, "(", round(both_pos/length(idx)*100, 3), "%)\n", sep="")
    cat("    IGF1R+ SOX2- :", igf1r_only, "(", round(igf1r_only/length(idx)*100, 3), "%)\n", sep="")
    cat("    IGF1R- SOX2+ :", sox2_only, "(", round(sox2_only/length(idx)*100, 3), "%)\n", sep="")
    cat("    IGF1R- SOX2- :", neither, "(", round(neither/length(idx)*100, 3), "%)\n", sep="")
  }

  # Global co-expression
  cat("\n--- GLOBAL CO-EXPRESSION ---\n")
  global_both <- sum(expr_matrix["IGF1R", ] > 0 & expr_matrix["SOX2", ] > 0)
  global_igf1r <- sum(expr_matrix["IGF1R", ] > 0)
  global_sox2 <- sum(expr_matrix["SOX2", ] > 0)
  global_total <- ncol(expr_matrix)

  cat("Total cells:", global_total, "\n")
  cat("IGF1R+ cells:", global_igf1r, "(", round(global_igf1r/global_total*100, 3), "%)\n", sep="")
  cat("SOX2+ cells:", global_sox2, "(", round(global_sox2/global_total*100, 3), "%)\n", sep="")
  cat("Both positive:", global_both, "(", round(global_both/global_total*100, 4), "%)\n", sep="")
  cat("Expected by chance (independent):", round(global_igf1r/global_total * global_sox2/global_total * global_total, 2), "cells\n")
  cat("Observed/Expected ratio:", round(global_both / (global_igf1r/global_total * global_sox2/global_total * global_total), 3), "\n")

  # Fisher's exact test for co-expression
  contingency <- matrix(c(
    global_both, global_igf1r - global_both,
    global_sox2 - global_both, global_total - global_igf1r - global_sox2 + global_both
  ), nrow = 2)

  cat("\nContingency table:\n")
  print(contingency)

  ft <- fisher.test(contingency)
  cat("\nFisher's exact test:\n")
  cat("  Odds ratio:", ft$estimate, "\n")
  cat("  p-value:", format(ft$p.value, digits = 3), "\n")
  cat("  95% CI:", paste(round(ft$conf.int, 3), collapse = " to "), "\n")

  if (ft$p.value < 0.05 && ft$estimate > 1) {
    cat("  Result: SIGNIFICANT POSITIVE ASSOCIATION (co-expression enriched)\n")
  } else if (ft$p.value < 0.05 && ft$estimate < 1) {
    cat("  Result: SIGNIFICANT NEGATIVE ASSOCIATION (mutual exclusion)\n")
  } else {
    cat("  Result: NO SIGNIFICANT ASSOCIATION\n")
  }
}

# 1.8 Age group analysis
cat("\n--- 1.8 AGE GROUP ANALYSIS ---\n")
age_groups <- cells$age_group
if (!is.null(age_groups)) {
  cat("Age groups:\n")
  print(table(age_groups, useNA = "ifany"))

  if ("IGF1R" %in% present_targets) {
    cat("\nIGF1R by age group:\n")
    for (ag in sort(unique(age_groups))) {
      idx <- which(age_groups == ag)
      if (length(idx) > 0) {
        vals <- expr_matrix["IGF1R", idx]
        cat(ag, ": mean=", round(mean(vals), 4), " %pos=", round(mean(vals > 0)*100, 2), " n=", length(idx), "\n", sep="")
      }
    }
  }

  if ("SOX2" %in% present_targets) {
    cat("\nSOX2 by age group:\n")
    for (ag in sort(unique(age_groups))) {
      idx <- which(age_groups == ag)
      if (length(idx) > 0) {
        vals <- expr_matrix["SOX2", idx]
        cat(ag, ": mean=", round(mean(vals), 4), " %pos=", round(mean(vals > 0)*100, 2), " n=", length(idx), "\n", sep="")
      }
    }
  }
}

# 1.9 Tissue origin analysis
cat("\n--- 1.9 TISSUE ORIGIN ANALYSIS ---\n")
tissue <- cells$tissue_origin
if (!is.null(tissue)) {
  cat("Tissue origins:\n")
  print(table(tissue, useNA = "ifany"))
}

# 1.10 Correlation analysis
cat("\n--- 1.10 PAIRWISE GENE CORRELATIONS ---\n")
if (length(present_targets) >= 2) {
  present_in_matrix <- present_targets[present_targets %in% rownames(expr_matrix)]
  if (length(present_in_matrix) >= 2) {
    cor_matrix <- cor(t(expr_matrix[present_in_matrix, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs")
    cat("\nSpearman correlation matrix (all cells):\n")
    print(round(cor_matrix, 3))

    if ("IGF1R" %in% present_in_matrix && "SOX2" %in% present_in_matrix) {
      cat("\nIGF1R vs SOX2 Spearman rho:", round(cor_matrix["IGF1R", "SOX2"], 4), "\n")
    }
    if ("IGF1R" %in% present_in_matrix && "IGFBP7" %in% present_in_matrix) {
      cat("IGF1R vs IGFBP7 Spearman rho:", round(cor_matrix["IGF1R", "IGFBP7"], 4), "\n")
    }
  }
}

rm(counts, norm_counts, log_expr, expr_matrix)
invisible(gc())


# ============================================================
# SECTION 2: GSE185222 COMPLETE ANALYSIS
# Dental Pulp Sound vs Caries, 10X Genomics
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 2: GSE185222 - DENTAL PULP SOUND VS CARIES\n")
cat("=" ,rep("=", 78), "\n", sep="")

# 2.1 Load data
cat("\n--- 2.1 DATA LOADING ---\n")
extract_dir <- "GSE185222_extracted"
csv_files <- list.files(extract_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
cat("CSV files found:", length(csv_files), "\n")

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
  cat("\nLoading:", basename(f), "\n")

  con <- gzfile(f, "r")
  df <- read.csv(con, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  close(con)

  cat("  Dimensions:", nrow(df), "genes x", ncol(df), "cells\n")
  cat("  First genes:", paste(rownames(df)[1:5], collapse = ", "), "\n")
  cat("  Max value:", max(as.matrix(df), na.rm = TRUE), "\n")

  gsm <- gsub("_.*$", "", basename(f))
  info <- sample_map[[gsm]]

  cat("  Sample:", info$id, "| Condition:", info$condition, "\n")

  # Data is already log-normalized (max ~8-9)
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

cat("\n--- 2.2 MERGED DATA QC ---\n")
cat("Total genes:", nrow(all_expr), "\n")
cat("Total cells:", ncol(all_expr), "\n")
cat("Condition distribution:\n")
print(table(all_meta$condition))

# 2.3 Target gene detection
cat("\n--- 2.3 TARGET GENE DETECTION ---\n")
target_genes_185 <- c("IGF1R", "SOX2", "IGF1", "IGFBP7", "IGFBP3", "IGFBP5", "PTN", "MDK", "FGFR1", "FGFR2", "BMPR1A", "BMPR1B", "BMP2", "BMP4", "WNT5A", "WNT10A", "AXIN2", "LEF1", "RUNX2", "SP7", "DSPP", "DMP1", "COL1A1", "THY1", "NT5E", "ENG", "MCAM", "CD146", "IL6", "CXCL8", "TNF", "NFKBIA")

avail_185 <- rownames(all_expr)
detection_185 <- data.frame(
  gene = target_genes_185,
  present = target_genes_185 %in% avail_185,
  stringsAsFactors = FALSE
)
print(detection_185)

present_185 <- detection_185$gene[detection_185$present]
cat("\nPresent:", length(present_185), "| Missing:", sum(!detection_185$present), "\n")

# 2.4 Expression by condition
cat("\n--- 2.4 EXPRESSION BY CONDITION ---\n")
conditions <- c("sound", "enamel_caries", "deep_caries")

for (gene in present_185) {
  cat("\n", gene, ":\n", sep="")

  cond_stats <- sapply(conditions, function(cond) {
    idx <- which(all_meta$condition == cond)
    if (length(idx) == 0) return(c(mean = NA, median = NA, sd = NA, pct_pos = NA, n = length(idx)))
    vals <- all_expr[gene, idx]
    c(
      mean = mean(vals, na.rm = TRUE),
      median = median(vals, na.rm = TRUE),
      sd = sd(vals, na.rm = TRUE),
      pct_pos = mean(vals > 0, na.rm = TRUE) * 100,
      n = length(vals)
    )
  })

  stats_df <- as.data.frame(t(cond_stats))
  stats_df$condition <- rownames(stats_df)
 num_cols <- sapply(stats_df, is.numeric)
stats_df[num_cols] <- round(stats_df[num_cols], 4)
print(stats_df)

  # Pairwise t-tests
  cat("  Pairwise t-tests (sound vs enamel, sound vs deep):\n")
  s_idx <- which(all_meta$condition == "sound")
  e_idx <- which(all_meta$condition == "enamel_caries")
  d_idx <- which(all_meta$condition == "deep_caries")

  if (length(s_idx) > 0 && length(e_idx) > 0) {
    tt_se <- t.test(all_expr[gene, s_idx], all_expr[gene, e_idx])
    cat("    sound vs enamel: p =", format(tt_se$p.value, digits = 3), 
        "| diff =", round(mean(all_expr[gene, e_idx]) - mean(all_expr[gene, s_idx]), 4), "\n")
  }
  if (length(s_idx) > 0 && length(d_idx) > 0) {
    tt_sd <- t.test(all_expr[gene, s_idx], all_expr[gene, d_idx])
    cat("    sound vs deep:   p =", format(tt_sd$p.value, digits = 3), 
        "| diff =", round(mean(all_expr[gene, d_idx]) - mean(all_expr[gene, s_idx]), 4), "\n")
  }
}

# 2.5 IGF1R vs IGFBP7 correlation by condition
cat("\n--- 2.5 IGF1R vs IGFBP7 CORRELATION ---\n")
if (all(c("IGF1R", "IGFBP7") %in% present_185)) {
  for (cond in conditions) {
    idx <- which(all_meta$condition == cond)
    if (length(idx) > 10) {
      cor_val <- cor(all_expr["IGF1R", idx], all_expr["IGFBP7", idx], method = "spearman", use = "complete.obs")
      cat(cond, ": Spearman rho =", round(cor_val, 4), " n =", length(idx), "\n")
    }
  }
}

# 2.6 SOX2 analysis
cat("\n--- 2.6 SOX2 ANALYSIS ---\n")
if ("SOX2" %in% present_185) {
  cat("SOX2 expression by condition:\n")
  for (cond in conditions) {
    idx <- which(all_meta$condition == cond)
    if (length(idx) > 0) {
      vals <- all_expr["SOX2", idx]
      cat(cond, ": mean=", round(mean(vals), 4), " %pos=", round(mean(vals > 0)*100, 2), " n=", length(idx), "\n", sep="")
    }
  }
}

rm(all_expr)
invisible(gc())


# ============================================================
# SECTION 3: GSE167251 SERIES MATRIX ANALYSIS
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 3: GSE167251 - SERIES MATRIX\n")
cat("=" ,rep("=", 78), "\n", sep="")

series_file <- "GSE167251/GSE167251_series_matrix.txt.gz"
if (file.exists(series_file)) {
  lines <- readLines(gzfile(series_file))

  cat("Total lines in series matrix:", length(lines), "\n")

  # Extract key metadata
  cat("\n--- Sample Characteristics ---\n")
  char_lines <- lines[grep("Sample_characteristics", lines)]
  for (l in char_lines) {
    cat(substr(l, 1, min(200, nchar(l))), "\n")
  }

  cat("\n--- Source Names ---\n")
  src_lines <- lines[grep("Sample_source_name", lines)]
  for (l in src_lines) {
    cat(substr(l, 1, min(200, nchar(l))), "\n")
  }

  cat("\n--- Titles ---\n")
  title_lines <- lines[grep("Sample_title", lines)]
  for (l in title_lines) {
    cat(substr(l, 1, min(200, nchar(l))), "\n")
  }
} else {
  cat("Series matrix not found\n")
}


# ============================================================
# SECTION 4: GSE132472 TAR ANALYSIS
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 4: GSE132472 - AMELOBLASTOMA BULK MICROARRAY\n")
cat("=" ,rep("=", 78), "\n", sep="")

tar_file <- "GSE132472_RAW.tar"
if (file.exists(tar_file)) {
  cat("TAR file size:", file.size(tar_file), "bytes\n")

  extract_132472 <- "GSE132472_extracted"
  if (!dir.exists(extract_132472)) {
    cat("Extracting...\n")
    untar(tar_file, exdir = extract_132472)
  }

  files_132472 <- list.files(extract_132472, recursive = TRUE, full.names = TRUE)
  cat("Extracted files:", length(files_132472), "\n")
  cat("First 10:\n")
  print(head(files_132472, 10))

  # Look for TXT files with expression
  txt_files <- files_132472[grepl("\\.txt$", files_132472)]
  cat("\nTXT files:", length(txt_files), "\n")

  if (length(txt_files) > 0) {
    # Read first file to check format
    first <- read.table(txt_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1, nrows = 5)
    cat("\nFirst file structure (first 5 rows):\n")
    print(first[, 1:min(4, ncol(first))])

    # Check for IGF1R
    full_first <- read.table(txt_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1)
    igf1r_rows <- grep("IGF1R", full_first[[1]], ignore.case = TRUE)
    cat("\nIGF1R probes in first sample:", length(igf1r_rows), "\n")
    if (length(igf1r_rows) > 0) {
      print(full_first[igf1r_rows, 1:min(5, ncol(full_first))])
    }
  }
} else {
  cat("TAR file not found\n")
}


# ============================================================
# SECTION 5: FINAL SYNTHESIS
# ============================================================

cat("\n")
cat("=" ,rep("=", 78), "\n", sep="")
cat("SECTION 5: FINAL SYNTHESIS\n")
cat("=" ,rep("=", 78), "\n", sep="")

cat("\nThis analysis is complete. All raw numbers are reported above.\n")
cat("No interpretations are made here - only the data.\n")
cat("\nKey datasets analyzed:\n")
cat("  1. GSE184749: 42,174 genes x 17,510 cells (dental mesenchyme)\n")
cat("  2. GSE185222: 20,293 genes x 6,582 cells (dental pulp sound/caries)\n")
cat("  3. GSE167251: Series matrix parsed (aging pulp)\n")
cat("  4. GSE132472: Bulk microarray TAR examined (ameloblastoma)\n")
cat("\nAll statistical tests, correlations, and expression values are in Sections 1-4.\n")
