################################################################################
# ULTRA-FAST BASE-R CROSS-VERIFICATION
# Vectorized, no cell-wise loops
# For: D:/!jenja/ARC/DATA3/
################################################################################

setwd("D:/!jenja/ARC/DATA3/")
library(Matrix)

TARGET_GENES <- c("IGF1R", "SOX2", "IGFBP7", "IGF1", "PTN", "MDK", "FGFR1", "BMPR1A")

# Vectorized log-normalization
log_norm_fast <- function(counts) {
  lib_size <- colSums(counts)
  sf <- median(lib_size[lib_size > 0])
  # Sparse-friendly: scale each column by its lib size, then multiply by sf
  norm <- counts
  norm@x <- norm@x / rep(lib_size, diff(norm@p)) * sf
  norm <- log1p(norm)
  return(norm)
}

# Read gzipped CSV
read_gz <- function(f) {
  con <- gzfile(f, "r")
  df <- read.csv(con, stringsAsFactors = FALSE)
  close(con)
  return(df)
}

################################################################################
# PART 1: GSE184749 — DENTAL MESENCHYME (THE CRITICAL ONE)
################################################################################

cat("\n========================================\n")
cat("GSE184749: DENTAL MESENCHYME\n")
cat("========================================\n")

counts_file <- "GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"
cell_file <- "GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"
gene_file <- "GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"

if (!all(file.exists(c(counts_file, cell_file, gene_file)))) {
  stop("Files missing!")
}

cat("Reading counts...\n")
counts <- readMM(counts_file)
cat("  Dim:", nrow(counts), "x", ncol(counts), "\n")

cat("Reading annotations...\n")
cells <- read_gz(cell_file)
genes <- read_gz(gene_file)
cat("  Cells:", nrow(cells), "| Genes:", nrow(genes), "\n")

# Set names
colnames(counts) <- cells[[1]]
gene_col <- grep("gene|symbol|name", colnames(genes), ignore.case = TRUE, value = TRUE)[1]
if (is.na(gene_col)) gene_col <- colnames(genes)[1]
rownames(counts) <- genes[[gene_col]]

cat("Normalizing (vectorized)...\n")
expr <- log_norm_fast(counts)

# Check genes
avail <- rownames(expr)
present <- intersect(TARGET_GENES, avail)
cat("\nPresent genes:", paste(present, collapse = ", "), "\n")
missing <- setdiff(TARGET_GENES, avail)
if (length(missing) > 0) cat("Missing:", paste(missing, collapse = ", "), "\n")

# Find cluster column
cluster_col <- grep("cluster|cell_type|annotation", colnames(cells), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(cluster_col)) {
  clusters <- cells[[cluster_col]]
  cat("\nCluster column:", cluster_col, "\n")
  tab <- table(clusters)
  print(tab)

  # IGF1R by cluster (vectorized)
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

  # Find DP cluster
  dp_patterns <- c("DP", "dental.papilla", "papilla", "DEM", "mesenchyme")
  dp_clusters <- c()
  for (pat in dp_patterns) {
    dp_clusters <- c(dp_clusters, grep(pat, unique(clusters), ignore.case = TRUE, value = TRUE))
  }
  dp_clusters <- unique(dp_clusters)
  cat("\nDP/progenitor clusters:", paste(dp_clusters, collapse = ", "), "\n")

  if (length(dp_clusters) > 0) {
    for (dp in dp_clusters) {
      idx <- which(clusters == dp)
      n <- length(idx)
      cat("\n=== CLUSTER", dp, "(", n, "cells) ===\n")

      if ("IGF1R" %in% present) {
        v <- expr["IGF1R", idx]
        cat("IGF1R  mean:", round(mean(v), 4), "| % > 0:", round(mean(v > 0)*100, 1), "\n")
      }
      if ("SOX2" %in% present) {
        v <- expr["SOX2", idx]
        cat("SOX2   mean:", round(mean(v), 4), "| % > 0:", round(mean(v > 0)*100, 1), "\n")
      }

      # CO-EXPRESSION (vectorized)
      if (all(c("IGF1R", "SOX2") %in% present)) {
        igf1r_pos <- expr["IGF1R", idx] > 0
        sox2_pos <- expr["SOX2", idx] > 0
        co_pos <- sum(igf1r_pos & sox2_pos)
        cat("\nCO-EXPRESSION:\n")
        cat("  IGF1R+ SOX2+ :", co_pos, "/", n, "(", round(co_pos/n*100, 1), "%)\n")
        cat("  IGF1R+ SOX2- :", sum(igf1r_pos & !sox2_pos), "\n")
        cat("  IGF1R- SOX2+ :", sum(!igf1r_pos & sox2_pos), "\n")

        if (co_pos > 0) {
          cat("\n>>> VERDICT: IGF1R + SOX2 CO-EXPRESSED in", dp, "<<<\n")
        } else {
          cat("\n>>> VERDICT: NO co-expression in", dp, "<<<\n")
        }
      }
    }
  }
} else {
  cat("No cluster column found!\n")
}

# Clean up memory
rm(counts, expr, cells, genes)
invisible(gc())


################################################################################
# PART 2: GSE185222 — SOUND VS CARIES
################################################################################

cat("\n\n========================================\n")
cat("GSE185222: SOUND VS CARIES\n")
cat("========================================\n")

tar_file <- "GSE185222_RAW/GSE185222_RAW.tar"

if (!file.exists(tar_file)) {
  cat("TAR not found\n")
} else {
  extract_dir <- "GSE185222_extracted"
  if (!dir.exists(extract_dir)) {
    cat("Extracting...\n")
    untar(tar_file, exdir = extract_dir)
  }

  sample_dirs <- list.dirs(extract_dir, recursive = TRUE)
  sample_dirs <- sample_dirs[sample_dirs != extract_dir]
  cat("Samples:", length(sample_dirs), "\n")

  # Read first sample to get gene list
  all_counts <- NULL
  all_cond <- c()

  for (sdir in sample_dirs) {
    bc <- list.files(sdir, pattern = "barcodes", full.names = TRUE)[1]
    ft <- list.files(sdir, pattern = "features|genes", full.names = TRUE)[1]
    mx <- list.files(sdir, pattern = "matrix\\.mtx", full.names = TRUE)[1]

    if (is.na(mx) || is.na(bc) || is.na(ft)) next

    cat("Loading:", basename(sdir), "\n")
    cts <- readMM(mx)
    bcdf <- read.table(bc, header = FALSE, stringsAsFactors = FALSE)
    ftdf <- read.table(ft, header = FALSE, stringsAsFactors = FALSE)

    colnames(cts) <- bcdf$V1
    rownames(cts) <- ftdf$V2

    gsm <- basename(sdir)
    cond <- ifelse(grepl("5593724", gsm), "sound",
            ifelse(grepl("5593725", gsm), "enamel_caries", "deep_caries"))

    if (is.null(all_counts)) {
      all_counts <- cts
    } else {
      cg <- intersect(rownames(all_counts), rownames(cts))
      all_counts <- all_counts[cg,]
      cts <- cts[cg,]
      all_counts <- cbind(all_counts, cts)
    }
    all_cond <- c(all_cond, rep(cond, ncol(cts)))
  }

  cat("Total:", ncol(all_counts), "cells | Conditions:", paste(table(all_cond), collapse = ", "), "\n")

  cat("Normalizing...\n")
  expr2 <- log_norm_fast(all_counts)
  rm(all_counts)
  invisible(gc())

  avail2 <- rownames(expr2)

  if ("IGF1R" %in% avail2) {
    cat("\n--- IGF1R by condition ---\n")
    for (cond in c("sound", "enamel_caries", "deep_caries")) {
      idx <- which(all_cond == cond)
      if (length(idx) > 0) {
        v <- expr2["IGF1R", idx]
        cat(cond, ": mean =", round(mean(v), 4), "| % pos =", round(mean(v > 0)*100, 1), "| n =", length(v), "\n")
      }
    }

    # T-test sound vs deep
    s_idx <- which(all_cond == "sound")
    d_idx <- which(all_cond == "deep_caries")
    if (length(s_idx) > 0 && length(d_idx) > 0) {
      tt <- t.test(expr2["IGF1R", s_idx], expr2["IGF1R", d_idx])
      cat("\nT-test sound vs deep:\n")
      cat("  Sound mean:", round(mean(expr2["IGF1R", s_idx]), 4), "\n")
      cat("  Deep mean:", round(mean(expr2["IGF1R", d_idx]), 4), "\n")
      cat("  p-value:", format(tt$p.value, digits = 3), "\n")
      if (tt$p.value < 0.05 && mean(expr2["IGF1R", d_idx]) > mean(expr2["IGF1R", s_idx])) {
        cat("\n>>> IGF1R HIGHER in deep caries (p < 0.05) <<<\n")
      } else if (tt$p.value < 0.05) {
        cat("\n>>> IGF1R LOWER in deep caries (p < 0.05) <<<\n")
      } else {
        cat("\n>>> No significant difference <<<\n")
      }
    }
  }
}


cat("\n\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
