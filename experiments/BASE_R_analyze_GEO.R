################################################################################
# BASE-R CROSS-VERIFICATION (NO Seurat, NO ggplot2)
# Works with clean R installation — only base packages needed
# For: D:/!jenja/ARC/DATA3/
################################################################################

setwd("D:/!jenja/ARC/DATA3/")

# Check if Matrix is available (should be, it's recommended package)
if (!requireNamespace("Matrix", quietly = TRUE)) {
  cat("Installing Matrix package...\n")
  install.packages("Matrix")
}
library(Matrix)

TARGET_GENES <- c("IGF1R", "SOX2", "IGFBP7", "IGF1", "PTN", "MDK", "FGFR1", "BMPR1A")

################################################################################
# HELPER FUNCTIONS
################################################################################

# Read gzipped CSV
read_gz_csv <- function(file) {
  con <- gzfile(file, "r")
  df <- read.csv(con, stringsAsFactors = FALSE)
  close(con)
  return(df)
}

# Simple log-normalization
log_norm <- function(counts) {
  lib_size <- colSums(counts)
  scale_factor <- median(lib_size)
  norm <- counts
  for (i in 1:ncol(counts)) {
    if (lib_size[i] > 0) {
      norm[,i] <- counts[,i] / lib_size[i] * scale_factor
    }
  }
  norm <- log1p(norm)
  return(norm)
}

# Mean expression by group
mean_by_group <- function(expr, groups) {
  u_groups <- unique(groups)
  result <- matrix(0, nrow = nrow(expr), ncol = length(u_groups))
  rownames(result) <- rownames(expr)
  colnames(result) <- u_groups

  for (g in u_groups) {
    idx <- which(groups == g)
    if (length(idx) > 0) {
      result[,g] <- rowMeans(expr[,idx, drop = FALSE])
    }
  }
  return(result)
}

# Percentage of cells expressing gene > 0
pct_express <- function(expr, groups) {
  u_groups <- unique(groups)
  result <- matrix(0, nrow = nrow(expr), ncol = length(u_groups))
  rownames(result) <- rownames(expr)
  colnames(result) <- u_groups

  for (g in u_groups) {
    idx <- which(groups == g)
    if (length(idx) > 0) {
      result[,g] <- rowMeans(expr[,idx, drop = FALSE] > 0) * 100
    }
  }
  return(result)
}

# T-test between two groups
ttest_gene <- function(expr, groups, g1, g2, gene) {
  idx1 <- which(groups == g1)
  idx2 <- which(groups == g2)

  if (length(idx1) == 0 || length(idx2) == 0) return(NULL)

  x1 <- expr[gene, idx1]
  x2 <- expr[gene, idx2]

  tt <- t.test(x1, x2)
  return(list(
    mean1 = mean(x1), mean2 = mean(x2),
    p_value = tt$p.value,
    log2fc = log2((mean(x2) + 0.01) / (mean(x1) + 0.01))
  ))
}

################################################################################
# PART 1: GSE184749 — DENTAL MESENCHYME
################################################################################

cat("\n========================================\n")
cat("GSE184749: DENTAL MESENCHYME\n")
cat("========================================\n")

counts_file <- "GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"
cell_annot_file <- "GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"
gene_annot_file <- "GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"

if (!all(file.exists(c(counts_file, cell_annot_file, gene_annot_file)))) {
  stop("GSE184749 files not found!")
}

cat("Loading counts...\n")
counts <- readMM(counts_file)

cat("Loading annotations...\n")
cells <- read_gz_csv(cell_annot_file)
genes <- read_gz_csv(gene_annot_file)

cat("Cells:", nrow(cells), "| Genes:", nrow(genes), "\n")
cat("Matrix:", nrow(counts), "x", ncol(counts), "\n")

# Set names
if (ncol(counts) == nrow(cells)) {
  colnames(counts) <- cells[[1]]
}
if (nrow(counts) == nrow(genes)) {
  gene_col <- grep("gene|symbol|name", colnames(genes), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(gene_col)) {
    rownames(counts) <- genes[[gene_col]]
  } else {
    rownames(counts) <- genes[[1]]
  }
}

# Normalize
expr <- log_norm(counts)

# Check target genes
avail <- rownames(expr)
cat("\n--- Gene Availability ---\n")
for (g in TARGET_GENES) {
  cat("  ", g, ifelse(g %in% avail, "✓", "✗"), "\n")
}

# Find cluster column
cluster_col <- grep("cluster|cell_type|annotation", colnames(cells), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(cluster_col)) {
  clusters <- cells[[cluster_col]]
  cat("\nCluster column:", cluster_col, "\n")
  print(table(clusters))

  # Mean expression by cluster
  if ("IGF1R" %in% avail) {
    cat("\n--- IGF1R mean expression by cluster ---\n")
    igf1r_means <- mean_by_group(expr["IGF1R",, drop = FALSE], clusters)
    print(sort(igf1r_means[1,], decreasing = TRUE))
  }

  if ("SOX2" %in% avail) {
    cat("\n--- SOX2 mean expression by cluster ---\n")
    sox2_means <- mean_by_group(expr["SOX2",, drop = FALSE], clusters)
    print(sort(sox2_means[1,], decreasing = TRUE))
  }

  # Find DP cluster
  dp_patterns <- c("DP", "dental papilla", "papilla", "DEM", "mesenchyme")
  dp_clusters <- c()
  for (pat in dp_patterns) {
    dp_clusters <- c(dp_clusters, grep(pat, unique(clusters), ignore.case = TRUE, value = TRUE))
  }
  dp_clusters <- unique(dp_clusters)

  cat("\n--- DP/Progenitor clusters found:", paste(dp_clusters, collapse = ", "), "---\n")

  if (length(dp_clusters) > 0) {
    for (dp in dp_clusters) {
      idx <- which(clusters == dp)
      cat("\nCluster", dp, "(" , length(idx), "cells):\n")

      if ("IGF1R" %in% avail) {
        igf1r_vals <- expr["IGF1R", idx]
        cat("  IGF1R mean:", mean(igf1r_vals), "| % > 0:", mean(igf1r_vals > 0)*100, "\n")
      }
      if ("SOX2" %in% avail) {
        sox2_vals <- expr["SOX2", idx]
        cat("  SOX2 mean:", mean(sox2_vals), "| % > 0:", mean(sox2_vals > 0)*100, "\n")
      }

      # CO-EXPRESSION
      if (all(c("IGF1R", "SOX2") %in% avail)) {
        co_pos <- sum(expr["IGF1R", idx] > 0 & expr["SOX2", idx] > 0)
        total <- length(idx)
        cat("  IGF1R+ SOX2+ co-expression:", co_pos, "/", total, 
            "(", round(co_pos/total*100, 1), "%)\n")

        if (co_pos > 0) {
          cat("\n  >>> VERDICT: IGF1R and SOX2 ARE CO-EXPRESSED in", dp, "<<<\n")
        } else {
          cat("\n  >>> VERDICT: NO co-expression in", dp, "<<<\n")
        }
      }
    }
  }
} else {
  cat("No cluster annotations found!\n")
}


################################################################################
# PART 2: GSE185222 — DENTAL PULP SOUND VS CARIES
################################################################################

cat("\n\n========================================\n")
cat("GSE185222: DENTAL PULP SOUND VS CARIES\n")
cat("========================================\n")

tar_file <- "GSE185222_RAW/GSE185222_RAW.tar"

if (!file.exists(tar_file)) {
  cat("ERROR: TAR not found\n")
} else {
  extract_dir <- "GSE185222_extracted"
  if (!dir.exists(extract_dir)) {
    cat("Extracting TAR...\n")
    untar(tar_file, exdir = extract_dir)
  }

  # Find all sample directories
  sample_dirs <- list.dirs(extract_dir, recursive = TRUE)
  sample_dirs <- sample_dirs[sample_dirs != extract_dir]

  cat("Samples found:", length(sample_dirs), "\n")

  # Merge all samples
  all_counts <- NULL
  all_conditions <- c()
  all_cells <- c()

  for (sdir in sample_dirs) {
    barcodes <- list.files(sdir, pattern = "barcodes", full.names = TRUE)[1]
    features <- list.files(sdir, pattern = "features|genes", full.names = TRUE)[1]
    matrix <- list.files(sdir, pattern = "matrix\\.mtx", full.names = TRUE)[1]

    if (!is.na(matrix) && !is.na(barcodes) && !is.na(features)) {
      cat("Loading:", basename(sdir), "\n")

      counts_s <- readMM(matrix)
      bc <- read.table(barcodes, header = FALSE, stringsAsFactors = FALSE)
      ft <- read.table(features, header = FALSE, stringsAsFactors = FALSE)

      colnames(counts_s) <- bc$V1
      rownames(counts_s) <- ft$V2

      # Determine condition
      gsm <- basename(sdir)
      cond <- ifelse(grepl("5593724", gsm), "sound",
              ifelse(grepl("5593725", gsm), "enamel_caries", "deep_caries"))

      if (is.null(all_counts)) {
        all_counts <- counts_s
      } else {
        # Match genes
        common_genes <- intersect(rownames(all_counts), rownames(counts_s))
        all_counts <- all_counts[common_genes,]
        counts_s <- counts_s[common_genes,]
        all_counts <- cbind(all_counts, counts_s)
      }

      all_conditions <- c(all_conditions, rep(cond, ncol(counts_s)))
      all_cells <- c(all_cells, colnames(counts_s))
    }
  }

  cat("Total cells:", ncol(all_counts), "\n")
  cat("Conditions:", table(all_conditions), "\n")

  # Normalize
  expr185 <- log_norm(all_counts)
  colnames(expr185) <- all_cells

  avail185 <- rownames(expr185)

  # IGF1R analysis
  if ("IGF1R" %in% avail185) {
    cat("\n--- IGF1R by condition ---\n")

    for (cond in c("sound", "enamel_caries", "deep_caries")) {
      idx <- which(all_conditions == cond)
      if (length(idx) > 0) {
        vals <- expr185["IGF1R", idx]
        cat(cond, ": mean =", round(mean(vals), 3), 
            "| % pos =", round(mean(vals > 0)*100, 1), 
            "| n =", length(vals), "\n")
      }
    }

    # T-test sound vs deep
    if (all(c("sound", "deep_caries") %in% all_conditions)) {
      tt_result <- ttest_gene(expr185, all_conditions, "sound", "deep_caries", "IGF1R")
      if (!is.null(tt_result)) {
        cat("\n--- T-test: sound vs deep caries ---\n")
        cat("Sound mean:", round(tt_result$mean1, 3), "\n")
        cat("Deep mean:", round(tt_result$mean2, 3), "\n")
        cat("Log2FC:", round(tt_result$log2fc, 3), "\n")
        cat("p-value:", format(tt_result$p_value, digits = 3), "\n")

        if (tt_result$p_value < 0.05) {
          if (tt_result$mean2 > tt_result$mean1) {
            cat("\n>>> IGF1R SIGNIFICANTLY HIGHER in deep caries <<<\n")
          } else {
            cat("\n>>> IGF1R SIGNIFICANTLY LOWER in deep caries <<<\n")
          }
        } else {
          cat("\n>>> No significant difference <<<\n")
        }
      }
    }
  }
}


################################################################################
# PART 3: GSE167251 — SERIES MATRIX ONLY
################################################################################

cat("\n\n========================================\n")
cat("GSE167251: SERIES MATRIX\n")
cat("========================================\n")

series_matrix <- "GSE167251/GSE167251_series_matrix.txt.gz"
if (file.exists(series_matrix)) {
  lines <- readLines(gzfile(series_matrix))

  cat("Sample characteristics:\n")
  char_lines <- lines[grep("Sample_characteristics", lines)]
  for (l in char_lines) {
    cat(substr(l, 1, min(150, nchar(l))), "\n")
  }

  cat("\nSource names:\n")
  src_lines <- lines[grep("Sample_source_name", lines)]
  for (l in src_lines) {
    cat(substr(l, 1, min(150, nchar(l))), "\n")
  }
} else {
  cat("Series matrix not found\n")
}


################################################################################
# PART 4: GSE132472 — BULK
################################################################################

cat("\n\n========================================\n")
cat("GSE132472: BULK MICROARRAY\n")
cat("========================================\n")

tar_132472 <- "GSE132472_RAW.tar"
if (file.exists(tar_132472)) {
  cat("TAR found. Extracting...\n")
  extract_132472 <- "GSE132472_extracted"
  if (!dir.exists(extract_132472)) {
    untar(tar_132472, exdir = extract_132472)
  }

  txt_files <- list.files(extract_132472, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
  cat("TXT files:", length(txt_files), "\n")

  if (length(txt_files) > 0) {
    sample1 <- read.table(txt_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1)
    cat("Columns:", colnames(sample1), "\n")

    if (ncol(sample1) >= 2) {
      igf1r_rows <- grep("IGF1R", sample1[[1]], ignore.case = TRUE)
      cat("IGF1R probes:", length(igf1r_rows), "\n")
      if (length(igf1r_rows) > 0) {
        print(sample1[igf1r_rows, 1:min(5, ncol(sample1))])
      }
    }
  }
}


################################################################################
# FINAL SUMMARY
################################################################################

cat("\n\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║              CROSS-VERIFICATION COMPLETE                           ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\nCheck console output above for VERDICT lines.\n")
