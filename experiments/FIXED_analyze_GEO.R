################################################################################
# FIXED CROSS-VERIFICATION SCRIPT
# Handles dgTMatrix, checks TAR structure, finds correct clusters
################################################################################

setwd("D:/!jenja/ARC/DATA3/")
library(Matrix)

TARGET_GENES <- c("IGF1R", "SOX2", "IGFBP7", "IGF1", "PTN", "MDK", "FGFR1", "BMPR1A")

# Read gzipped CSV
read_gz <- function(f) {
  con <- gzfile(f, "r")
  df <- read.csv(con, stringsAsFactors = FALSE)
  close(con)
  return(df)
}

# Vectorized log-norm - handles any sparse format
log_norm_fast <- function(counts) {
  # Convert to dgCMatrix if needed
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  lib_size <- colSums(counts)
  sf <- median(lib_size[lib_size > 0])
  # Scale: divide each column by lib_size, multiply by sf
  counts@x <- counts@x / rep(lib_size, diff(counts@p)) * sf
  counts <- log1p(counts)
  return(counts)
}

################################################################################
# PART 1: GSE184749 — DENTAL MESENCHYME
################################################################################

cat("\n========================================\n")
cat("GSE184749: DENTAL MESENCHYME\n")
cat("========================================\n")

counts_file <- "GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"
cell_file <- "GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"
gene_file <- "GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"

counts <- readMM(counts_file)
cat("Raw matrix class:", class(counts), "\n")

# Convert to dgCMatrix
counts <- as(counts, "dgCMatrix")
cat("Converted to:", class(counts), "\n")

cells <- read_gz(cell_file)
genes <- read_gz(gene_file)

colnames(counts) <- cells[[1]]
gene_col <- grep("gene|symbol|name", colnames(genes), ignore.case = TRUE, value = TRUE)[1]
if (is.na(gene_col)) gene_col <- colnames(genes)[1]
rownames(counts) <- genes[[gene_col]]

cat("Cells:", nrow(cells), "| Genes:", nrow(genes), "\n")
cat("Matrix:", nrow(counts), "x", ncol(counts), "\n")

cat("\nNormalizing...\n")
expr <- log_norm_fast(counts)

# Check genes
avail <- rownames(expr)
present <- intersect(TARGET_GENES, avail)
cat("\nPresent:", paste(present, collapse = ", "), "\n")

# Check ALL columns in cell annotations
cat("\nCell annotation columns:\n")
print(colnames(cells))

# Check for any column that might have DP or stage info
for (col in colnames(cells)) {
  vals <- cells[[col]]
  if (is.character(vals) || is.factor(vals)) {
    uvals <- unique(vals)
    if (length(uvals) <= 20) {
      cat("\nColumn", col, "(", length(uvals), "unique):\n")
      print(table(vals))
    }
  }
}

# Primary cluster column
cluster_col <- grep("cluster|cell_type|annotation", colnames(cells), ignore.case = TRUE, value = TRUE)[1]
clusters <- cells[[cluster_col]]
cat("\nUsing cluster column:", cluster_col, "\n")
print(table(clusters))

# IGF1R by cluster
if ("IGF1R" %in% present) {
  cat("\n--- IGF1R mean by cluster ---\n")
  u_cl <- unique(clusters)
  igf1r_means <- sapply(u_cl, function(g) mean(expr["IGF1R", clusters == g]))
  names(igf1r_means) <- u_cl
  print(sort(igf1r_means, decreasing = TRUE))
}

# SOX2 by cluster
if ("SOX2" %in% present) {
  cat("\n--- SOX2 mean by cluster ---\n")
  sox2_means <- sapply(unique(clusters), function(g) mean(expr["SOX2", clusters == g]))
  names(sox2_means) <- unique(clusters)
  print(sort(sox2_means, decreasing = TRUE))
}

# Find clusters with ANY SOX2 expression
cat("\n--- SOX2+ cells by cluster ---\n")
for (cl in unique(clusters)) {
  idx <- which(clusters == cl)
  sox2_pos <- sum(expr["SOX2", idx] > 0)
  cat(cl, ":", sox2_pos, "/", length(idx), "(", round(sox2_pos/length(idx)*100, 2), "%)\n")
}

# Check ALL cells for IGF1R + SOX2 co-expression (not just DEM)
cat("\n--- GLOBAL IGF1R + SOX2 co-expression ---\n")
igf1r_all <- expr["IGF1R", ] > 0
sox2_all <- expr["SOX2", ] > 0
co_global <- sum(igf1r_all & sox2_all)
cat("IGF1R+ SOX2+ in ALL cells:", co_global, "/", length(clusters), 
    "(", round(co_global/length(clusters)*100, 3), "%)\n")

# Check by cluster
for (cl in unique(clusters)) {
  idx <- which(clusters == cl)
  co <- sum(expr["IGF1R", idx] > 0 & expr["SOX2", idx] > 0)
  cat(cl, ":", co, "/", length(idx), "\n")
}

# Check if there's a pseudotime or stage column that separates DP from DEM
stage_col <- grep("stage|time|day|pseudotime|week", colnames(cells), ignore.case = TRUE, value = TRUE)
if (length(stage_col) > 0) {
  cat("\nStage/time columns found:", paste(stage_col, collapse = ", "), "\n")
}

# CRITICAL: In dental mesenchyme, DEM = common progenitor, DP = dental papilla
# If DP is not a separate cluster, it may be SUBSET of DEM cells (early stage)
# Check if there's any metadata that distinguishes early DEM (DP-like) from late DEM

# Look for tooth germ type (incisor vs molar)
tooth_col <- grep("tooth|incisor|molar|type", colnames(cells), ignore.case = TRUE, value = TRUE)
if (length(tooth_col) > 0) {
  cat("\nTooth type columns:", paste(tooth_col, collapse = ", "), "\n")
  for (tc in tooth_col) {
    cat("\n", tc, ":\n")
    print(table(cells[[tc]]))
  }
}

# VERDICT
if (co_global > 0) {
  cat("\n>>> VERDICT: IGF1R + SOX2 CO-EXPRESSION EXISTS (", co_global, "cells) <<<\n")
} else {
  cat("\n>>> VERDICT: NO IGF1R + SOX2 co-expression in this dataset <<<\n")
}

# Check if IGF1R is highest in DEM (progenitor) vs differentiated
if ("IGF1R" %in% present) {
  cat("\n--- IGF1R ranking ---\n")
  igf_sorted <- sort(igf1r_means, decreasing = TRUE)
  cat("Highest IGF1R:", names(igf_sorted)[1], "(", round(igf_sorted[1], 4), ")\n")
  if (names(igf_sorted)[1] %in% c("DEM", "DP")) {
    cat(">>> IGF1R is highest in PROGENITOR cluster <<<\n")
  }
}

rm(counts, expr, cells, genes)
invisible(gc())


################################################################################
# PART 2: GSE185222 — CHECK TAR STRUCTURE
################################################################################

cat("\n\n========================================\n")
cat("GSE185222: TAR STRUCTURE CHECK\n")
cat("========================================\n")

tar_file <- "GSE185222_RAW/GSE185222_RAW.tar"

if (!file.exists(tar_file)) {
  cat("TAR not found!\n")
} else {
  cat("TAR found:", tar_file, "\n")

  # List contents without extracting
  cat("\nListing TAR contents...\n")
  tar_list <- untar(tar_file, list = TRUE)
  cat("Total entries:", length(tar_list), "\n")
  cat("First 20 entries:\n")
  print(head(tar_list, 20))

  # Check directory structure
  dirs <- unique(dirname(tar_list))
  cat("\nDirectories in TAR:\n")
  print(head(dirs, 20))

  # Extract
  extract_dir <- "GSE185222_extracted"
  if (!dir.exists(extract_dir)) {
    cat("\nExtracting...\n")
    untar(tar_file, exdir = extract_dir)
    cat("Done.\n")
  }

  # List extracted files
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  cat("\nExtracted files:", length(all_files), "\n")
  cat("First 30:\n")
  print(head(all_files, 30))

  # Find barcodes/features/matrix
  barcodes <- all_files[grepl("barcodes", all_files)]
  features <- all_files[grepl("features|genes", all_files)]
  matrix_files <- all_files[grepl("matrix\\.mtx", all_files)]

  cat("\nBarcodes:", length(barcodes), "\n")
  cat("Features:", length(features), "\n")
  cat("Matrix files:", length(matrix_files), "\n")

  if (length(matrix_files) > 0) {
    cat("\nMatrix files found:\n")
    print(matrix_files)

    # Load first sample
    mx <- matrix_files[1]
    bc <- barcodes[1]
    ft <- features[1]

    cat("\nLoading first sample:", basename(dirname(mx)), "\n")
    cts <- readMM(mx)
    cts <- as(cts, "dgCMatrix")
    bcdf <- read.table(bc, header = FALSE, stringsAsFactors = FALSE)
    ftdf <- read.table(ft, header = FALSE, stringsAsFactors = FALSE)

    colnames(cts) <- bcdf$V1
    rownames(cts) <- ftdf$V2

    cat("Sample dim:", nrow(cts), "x", ncol(cts), "\n")

    # Check IGF1R
    if ("IGF1R" %in% rownames(cts)) {
      cat("IGF1R present! Mean counts:", mean(cts["IGF1R",]), "\n")
      cat("IGF1R % positive:", mean(cts["IGF1R",] > 0)*100, "\n")
    }

    # Check SOX2
    if ("SOX2" %in% rownames(cts)) {
      cat("SOX2 present! Mean counts:", mean(cts["SOX2",]), "\n")
    }
  }
}


cat("\n\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
