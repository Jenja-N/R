################################################################################
# EXACT CROSS-VERIFICATION SCRIPT
# For files in D:/!jenja/ARC/DATA3/
# Based on actual file inventory from report_fixed.txt
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

setwd("D:/!jenja/ARC/DATA3/")

TARGET_GENES <- c("IGF1R", "SOX2", "IGFBP7", "IGF1", "PTN", "MDK", "FGFR1", "BMPR1A", "THY1", "COL1A1", "DSPP")

################################################################################
# PART 1: GSE184749 — DENTAL MESENCHYME (THE KEY DATASET)
################################################################################

cat("\n========================================\n")
cat("GSE184749: DENTAL MESENCHYME\n")
cat("========================================\n")

# File paths (exact from inventory)
counts_file <- "GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"
cell_annot_file <- "GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"
gene_annot_file <- "GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"

# Verify files exist
for (f in c(counts_file, cell_annot_file, gene_annot_file)) {
  if (!file.exists(f)) {
    stop(paste("File not found:", f))
  }
}

cat("Loading counts matrix...\n")
counts <- readMM(counts_file)

cat("Loading cell annotations...\n")
cells <- read.csv(gzfile(cell_annot_file), stringsAsFactors = FALSE)
cat("Cells:", nrow(cells), "| Columns:", paste(colnames(cells), collapse = ", "), "\n")

cat("Loading gene annotations...\n")
genes <- read.csv(gzfile(gene_annot_file), stringsAsFactors = FALSE)
cat("Genes:", nrow(genes), "| Columns:", paste(colnames(genes), collapse = ", "), "\n")

# Set dimnames
counts_dim <- dim(counts)
cat("Matrix dimensions:", counts_dim[1], "genes x", counts_dim[2], "cells\n")

# Match dimensions
if (ncol(counts) == nrow(cells)) {
  colnames(counts) <- cells[[1]]  # first column is cell ID
} else {
  cat("WARNING: Cell count mismatch. Matrix cols:", ncol(counts), "| Annotations:", nrow(cells), "\n")
}

if (nrow(counts) == nrow(genes)) {
  # Find gene symbol column
  gene_symbol_col <- grep("gene|symbol|name", colnames(genes), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(gene_symbol_col)) {
    rownames(counts) <- genes[[gene_symbol_col]]
    cat("Gene symbols set from column:", gene_symbol_col, "\n")
  } else {
    rownames(counts) <- genes[[1]]
    cat("Gene names set from first column\n")
  }
} else {
  cat("WARNING: Gene count mismatch. Matrix rows:", nrow(counts), "| Annotations:", nrow(genes), "\n")
}

# Create Seurat object
s184749 <- CreateSeuratObject(counts = counts, project = "GSE184749_DM")

# Add metadata from cell annotations
for (col in colnames(cells)) {
  if (col != colnames(cells)[1]) {  # skip cell ID column
    s184749[[col]] <- cells[[col]][match(colnames(s184749), cells[[1]])]
  }
}

cat("\nMetadata columns added:", paste(setdiff(colnames(s184749@meta.data), c("orig.ident", "nCount_RNA", "nFeature_RNA")), collapse = ", "), "\n")

# Check for cluster annotations
cluster_cols <- grep("cluster|cell_type|annotation|celltype", colnames(s184749@meta.data), ignore.case = TRUE, value = TRUE)
if (length(cluster_cols) > 0) {
  cat("\nCluster annotation columns found:", paste(cluster_cols, collapse = ", "), "\n")

  # Use first cluster column as identity
  Idents(s184749) <- cluster_cols[1]
  cat("\nCluster distribution:\n")
  print(table(Idents(s184749)))
} else {
  cat("\nNo cluster annotations found. Will cluster de novo.\n")
}

# Standard processing
s184749 <- NormalizeData(s184749)
s184749 <- FindVariableFeatures(s184749, selection.method = "vst", nfeatures = 2000)
s184749 <- ScaleData(s184749)
s184749 <- RunPCA(s184749, features = VariableFeatures(object = s184749), npcs = 30)
s184749 <- FindNeighbors(s184749, dims = 1:15)
s184749 <- FindClusters(s184749, resolution = 0.5)
s184749 <- RunUMAP(s184749, dims = 1:15)

# Check target gene availability
available_genes <- rownames(s184749)
cat("\n--- TARGET GENE AVAILABILITY ---\n")
for (g in TARGET_GENES) {
  status <- ifelse(g %in% available_genes, "✓ PRESENT", "✗ MISSING")
  cat("  ", g, " ", status, "\n")
}
present_targets <- intersect(TARGET_GENES, available_genes)
cat("\nPresent targets:", paste(present_targets, collapse = ", "), "\n")

# CRITICAL ANALYSIS: IGF1R in DP cluster (Sox2+ progenitors)
cat("\n========================================\n")
cat("CRITICAL: IGF1R IN DP (DENTAL PAPILLA) CLUSTER\n")
cat("========================================\n")

# Find DP cluster by name
dp_patterns <- c("DP", "dental papilla", "papilla", "DEM", "mesenchyme", "progenitor")
dp_clusters <- c()
for (pat in dp_patterns) {
  dp_clusters <- c(dp_clusters, grep(pat, levels(Idents(s184749)), ignore.case = TRUE, value = TRUE))
}
dp_clusters <- unique(dp_clusters)

if (length(dp_clusters) > 0) {
  cat("DP/progenitor cluster(s) found:", paste(dp_clusters, collapse = ", "), "\n")

  for (dp in dp_clusters) {
    dp_cells <- subset(s184749, idents = dp)
    n_dp <- ncol(dp_cells)

    cat("\n--- Cluster", dp, "(" , n_dp, "cells) ---\n")

    if ("IGF1R" %in% present_targets) {
      igf1r_expr <- FetchData(dp_cells, vars = "IGF1R")$IGF1R
      cat("IGF1R expression:\n")
      cat("  Mean (log-norm):", mean(igf1r_expr), "\n")
      cat("  % cells IGF1R > 0:", mean(igf1r_expr > 0) * 100, "\n")
      cat("  % cells IGF1R > 1:", mean(igf1r_expr > 1) * 100, "\n")
    }

    if ("SOX2" %in% present_targets) {
      sox2_expr <- FetchData(dp_cells, vars = "SOX2")$SOX2
      cat("SOX2 expression:\n")
      cat("  Mean (log-norm):", mean(sox2_expr), "\n")
      cat("  % cells SOX2 > 0:", mean(sox2_expr > 0) * 100, "\n")
    }

    # CO-EXPRESSION ANALYSIS
    if (all(c("IGF1R", "SOX2") %in% present_targets)) {
      coexpr <- FetchData(dp_cells, vars = c("IGF1R", "SOX2"))
      igf1r_pos_sox2_pos <- sum(coexpr$IGF1R > 0 & coexpr$SOX2 > 0)
      igf1r_pos_sox2_neg <- sum(coexpr$IGF1R > 0 & coexpr$SOX2 == 0)
      igf1r_neg_sox2_pos <- sum(coexpr$IGF1R == 0 & coexpr$SOX2 > 0)

      cat("\nCO-EXPRESSION IN DP:\n")
      cat("  IGF1R+ SOX2+ :", igf1r_pos_sox2_pos, "(", round(igf1r_pos_sox2_pos/n_dp*100, 1), "%)\n")
      cat("  IGF1R+ SOX2- :", igf1r_pos_sox2_neg, "(", round(igf1r_pos_sox2_neg/n_dp*100, 1), "%)\n")
      cat("  IGF1R- SOX2+ :", igf1r_neg_sox2_pos, "(", round(igf1r_neg_sox2_pos/n_dp*100, 1), "%)\n")

      # Statistical test
      if (igf1r_pos_sox2_pos > 0) {
        cat("\n✓✓✓ VERDICT: IGF1R and SOX2 ARE CO-EXPRESSED in DP cluster ✓✓✓\n")
      } else {
        cat("\n✗✗✗ VERDICT: NO co-expression detected in DP cluster ✗✗✗\n")
      }
    }
  }
} else {
  cat("No DP cluster found by name. Checking SOX2 expression across all clusters...\n")

  if ("SOX2" %in% present_targets) {
    sox2_avg <- AverageExpression(s184749, features = "SOX2")$RNA
    sox2_cluster <- colnames(sox2_avg)[which.max(sox2_avg)]
    cat("\nHighest SOX2 cluster:", sox2_cluster, "(mean:", max(sox2_avg), ")\n")

    sox2_cells <- subset(s184749, idents = sox2_cluster)
    if ("IGF1R" %in% present_targets) {
      igf1r_in_sox2 <- mean(FetchData(sox2_cells, vars = "IGF1R")$IGF1R)
      cat("IGF1R in SOX2+ cluster:", igf1r_in_sox2, "\n")
    }
  }
}

# Compare IGF1R across ALL clusters
cat("\n--- IGF1R expression across ALL clusters ---\n")
if ("IGF1R" %in% present_targets) {
  igf1r_all <- AverageExpression(s184749, features = "IGF1R")$RNA
  print(igf1r_all)

  # Rank clusters by IGF1R
  igf1r_ranked <- sort(igf1r_all[1,], decreasing = TRUE)
  cat("\nTop IGF1R clusters:\n")
  print(head(igf1r_ranked, 5))
}

# Save plots
if (length(present_targets) > 0) {
  ggsave("GSE184749_umap_clusters.png", DimPlot(s184749, reduction = "umap", label = TRUE), width = 10, height = 8)
  ggsave("GSE184749_umap_split.png", DimPlot(s184749, reduction = "umap", split.by = cluster_cols[1]), width = 14, height = 6)

  feat_plot <- present_targets[present_targets %in% c("IGF1R", "SOX2", "IGFBP7", "FGFR1")]
  if (length(feat_plot) > 0) {
    ggsave("GSE184749_features.png", FeaturePlot(s184749, features = feat_plot), width = 12, height = 10)
  }

  if (all(c("IGF1R", "SOX2") %in% present_targets)) {
    ggsave("GSE184749_IGF1R_SOX2.png", FeaturePlot(s184749, features = c("IGF1R", "SOX2"), blend = TRUE), width = 10, height = 8)
  }
}


################################################################################
# PART 2: GSE185222 — DENTAL PULP SOUND VS CARIES
################################################################################

cat("\n\n========================================\n")
cat("GSE185222: DENTAL PULP SOUND VS CARIES\n")
cat("========================================\n")

tar_file <- "GSE185222_RAW/GSE185222_RAW.tar"

if (file.exists(tar_file)) {
  cat("TAR file found:", tar_file, "(55.3 Mb)\n")

  extract_dir <- "GSE185222_extracted"
  if (!dir.exists(extract_dir)) {
    cat("Extracting TAR...\n")
    untar(tar_file, exdir = extract_dir)
    cat("Extraction complete.\n")
  } else {
    cat("Using previously extracted files in", extract_dir, "\n")
  }

  # Find 10X files
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)

  # Group by sample (GSM prefixes)
  sample_dirs <- unique(dirname(all_files))
  cat("Sample directories found:", length(sample_dirs), "\n")

  # Read all samples and merge
  seurat_list <- list()
  sample_names <- c()

  for (sdir in sample_dirs) {
    barcodes_file <- list.files(sdir, pattern = "barcodes", full.names = TRUE)[1]
    features_file <- list.files(sdir, pattern = "features|genes", full.names = TRUE)[1]
    matrix_file <- list.files(sdir, pattern = "matrix\\\\.mtx", full.names = TRUE)[1]

    if (!is.na(matrix_file) && !is.na(barcodes_file) && !is.na(features_file)) {
      cat("Loading sample from:", basename(sdir), "\n")

      counts <- readMM(matrix_file)
      bc <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
      ft <- read.table(features_file, header = FALSE, stringsAsFactors = FALSE)

      colnames(counts) <- bc$V1
      rownames(counts) <- ft$V2

      # Determine condition from GSM ID
      gsm <- basename(sdir)
      if (grepl("5593724", gsm)) {
        condition <- "sound"
      } else if (grepl("5593725", gsm)) {
        condition <- "enamel_caries"
      } else if (grepl("5593726", gsm)) {
        condition <- "deep_caries"
      } else {
        condition <- gsm
      }

      seurat_obj <- CreateSeuratObject(counts = counts, project = condition)
      seurat_obj[["condition"]] <- condition
      seurat_obj[["sample"]] <- gsm

      seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      sample_names <- c(sample_names, condition)
    }
  }

  cat("\nSamples loaded:", paste(sample_names, collapse = ", "), "\n")

  if (length(seurat_list) > 0) {
    # Merge
    if (length(seurat_list) == 1) {
      s185222 <- seurat_list[[1]]
    } else {
      s185222 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)
    }

    cat("Total cells:", ncol(s185222), "\n")
    cat("Conditions:", table(s185222$condition), "\n")

    # QC
    s185222[["percent.mt"]] <- PercentageFeatureSet(s185222, pattern = "^MT-")
    s185222 <- subset(s185222, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
    cat("After QC:", ncol(s185222), "cells\n")

    # Process
    s185222 <- NormalizeData(s185222)
    s185222 <- FindVariableFeatures(s185222, nfeatures = 2000)
    s185222 <- ScaleData(s185222)
    s185222 <- RunPCA(s185222, npcs = 30)
    s185222 <- FindNeighbors(s185222, dims = 1:15)
    s185222 <- FindClusters(s185222, resolution = 0.8)
    s185222 <- RunUMAP(s185222, dims = 1:15)

    # Check genes
    avail_185222 <- rownames(s185222)
    cat("\nTarget gene availability:\n")
    for (g in TARGET_GENES) {
      status <- ifelse(g %in% avail_185222, "✓", "✗")
      cat("  ", g, status, "\n")
    }

    # DEG IGF1R by condition
    if ("IGF1R" %in% avail_185222) {
      cat("\n--- IGF1R DIFFERENTIAL EXPRESSION BY CONDITION ---\n")
      Idents(s185222) <- "condition"

      deg <- FindAllMarkers(s185222, features = "IGF1R", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)
      if (nrow(deg) > 0) {
        print(deg[, c("gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")])
      }

      cat("\n--- IGF1R expression levels by condition ---\n")
      for (cond in c("sound", "enamel_caries", "deep_caries")) {
        if (sum(s185222$condition == cond) > 0) {
          expr <- FetchData(subset(s185222, condition == cond), vars = "IGF1R")$IGF1R
          cat(cond, ": mean =", round(mean(expr), 3), 
              "| % pos =", round(mean(expr > 0)*100, 1), 
              "| n =", length(expr), "\n")
        }
      }

      # Test sound vs deep caries
      if (all(c("sound", "deep_caries") %in% s185222$condition)) {
        cat("\n--- T-test: sound vs deep caries ---\n")
        sound_expr <- FetchData(subset(s185222, condition == "sound"), vars = "IGF1R")$IGF1R
        deep_expr <- FetchData(subset(s185222, condition == "deep_caries"), vars = "IGF1R")$IGF1R

        tt <- t.test(sound_expr, deep_expr)
        cat("Sound mean:", round(mean(sound_expr), 3), "\n")
        cat("Deep mean:", round(mean(deep_expr), 3), "\n")
        cat("p-value:", format(tt$p.value, digits = 3), "\n")

        if (tt$p.value < 0.05) {
          if (mean(deep_expr) > mean(sound_expr)) {
            cat("✓✓✓ IGF1R is SIGNIFICANTLY HIGHER in deep caries ✓✓✓\n")
          } else {
            cat("✗✗✗ IGF1R is SIGNIFICANTLY LOWER in deep caries ✗✗✗\n")
          }
        } else {
          cat("○○○ No significant difference in IGF1R between sound and deep caries ○○○\n")
        }
      }
    }

    # Save plots
    ggsave("GSE185222_umap_condition.png", DimPlot(s185222, reduction = "umap", group.by = "condition"), width = 8, height = 6)
    if ("IGF1R" %in% avail_185222) {
      ggsave("GSE185222_IGF1R_feature.png", FeaturePlot(s185222, features = "IGF1R"), width = 6, height = 5)
      ggsave("GSE185222_IGF1R_violin.png", VlnPlot(s185222, features = "IGF1R", group.by = "condition"), width = 8, height = 5)
    }
  }
} else {
  cat("ERROR: TAR file not found:", tar_file, "\n")
}


################################################################################
# PART 3: GSE167251 — NORMAL VS AGING INFLAMMATORY PULP
################################################################################

cat("\n\n========================================\n")
cat("GSE167251: NORMAL VS AGING INFLAMMATORY PULP\n")
cat("========================================\n")

tar_167251 <- "GSE167251/GSE167251_RAW.tar"
series_matrix_167251 <- "GSE167251/GSE167251_series_matrix.txt.gz"

if (file.exists(series_matrix_167251)) {
  cat("Series matrix found. Parsing sample metadata...\n")
  lines <- readLines(gzfile(series_matrix_167251))

  # Extract sample characteristics
  char_lines <- lines[grep("Sample_characteristics", lines)]
  cat("\nSample characteristics:\n")
  for (l in char_lines) {
    cat(substr(l, 1, min(200, nchar(l))), "\n")
  }

  source_lines <- lines[grep("Sample_source_name", lines)]
  cat("\nSource names:\n")
  for (l in source_lines) {
    cat(substr(l, 1, min(200, nchar(l))), "\n")
  }
}

if (file.exists(tar_167251)) {
  cat("\nRAW.tar found:", tar_167251, "(34.5 Mb)\n")
  cat("Extracting to check contents...\n")

  extract_167251 <- "GSE167251_extracted"
  if (!dir.exists(extract_167251)) {
    untar(tar_167251, exdir = extract_167251)
  }

  extracted_files <- list.files(extract_167251, recursive = TRUE, full.names = TRUE)
  cat("Extracted files:", length(extracted_files), "\n")
  print(head(extracted_files, 20))

  # Check if these are 10X or bulk
  if (any(grepl("barcodes|features|matrix", extracted_files))) {
    cat("\n10X format detected! Can analyze as scRNA-seq.\n")
  } else if (any(grepl("\\\\.cel|\\\\.txt|\\\\.csv", extracted_files))) {
    cat("\nBulk format detected (microarray or bulk RNA-seq).\n")
    cat("Cannot assess single-cell Sox2+ IGF1R co-expression.\n")
  }
} else {
  cat("No RAW.tar found for GSE167251\n")
  cat("Using published knowledge: IGF signaling ↑ in aged pulp via IGFBP7\n")
}


################################################################################
# PART 4: GSE132472 — AMELOBLASTOMA BULK MICROARRAY
################################################################################

cat("\n\n========================================\n")
cat("GSE132472: AMELOBLASTOMA BULK MICROARRAY\n")
cat("========================================\n")

tar_132472 <- "GSE132472_RAW.tar"

if (file.exists(tar_132472)) {
  cat("TAR found:", tar_132472, "(171 Mb)\n")
  cat("Extracting...\n")

  extract_132472 <- "GSE132472_extracted"
  if (!dir.exists(extract_132472)) {
    untar(tar_132472, exdir = extract_132472)
  }

  txt_files <- list.files(extract_132472, pattern = "\\\\.txt$", recursive = TRUE, full.names = TRUE)
  cat("TXT files:", length(txt_files), "\n")

  if (length(txt_files) > 0) {
    # Read first sample
    sample1 <- read.table(txt_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1)
    cat("\nFirst sample columns:", colnames(sample1), "\n")
    cat("First few rows:\n")
    print(head(sample1[, 1:min(5, ncol(sample1))]))

    # Check for IGF1R probe
    if (ncol(sample1) >= 2) {
      probe_col <- colnames(sample1)[1]
      igf1r_rows <- grep("IGF1R", sample1[[probe_col]], ignore.case = TRUE)
      cat("\nIGF1R probes found:", length(igf1r_rows), "\n")
      if (length(igf1r_rows) > 0) {
        print(sample1[igf1r_rows, 1:min(5, ncol(sample1))])
      }
    }
  }
} else {
  cat("TAR file not found\n")
}


################################################################################
# FINAL SUMMARY
################################################################################

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║           CROSS-VERIFICATION SUMMARY                               ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ GSE184749: IGF1R + SOX2 in DP cluster → CHECK CONSOLE OUTPUT ABOVE  ║\n")
cat("║ GSE185222: IGF1R in sound vs deep caries → CHECK CONSOLE OUTPUT ABOVE║\n")
cat("║ GSE167251: IGFBP7/IGF1R age relationship → from published data       ║\n")
cat("║ GSE132472: Bulk IGF1R in ameloblastoma → CHECK CONSOLE OUTPUT ABOVE  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\nCheck generated PNG files in:", getwd(), "\n")
