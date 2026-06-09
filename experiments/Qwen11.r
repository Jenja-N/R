#==============================================================================
# ГДЕ РЕАЛЬНО DSPP+ КЛЕТКИ? Поиск истинного пути к одонтобластам
# Только зубы. Только данные.
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")
library(Matrix)

# Загрузка данных (если не загружены)
if(!exists("fetal_log") || !exists("adult_log") || !exists("cell_annot")) {
  cat("Загружаю данные...\n")
  
  cell_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"))
  gene_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"))
  mtx <- readMM(gzfile("GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"))
  rownames(mtx) <- gene_annot$gene_short_name
  colnames(mtx) <- cell_annot$cell_id
  mtx <- mtx[!duplicated(rownames(mtx)), ]
  
  fetal_lib_size <- Matrix::colSums(mtx)
  fetal_cpm <- t(t(mtx) / fetal_lib_size) * 1e6
  fetal_log <- as.matrix(log1p(fetal_cpm))
  
  adult_data <- readRDS("GSE185222_processed.RDS")
  adult_log <- adult_data$expr
  adult_conditions <- gsub("_.*", "", adult_data$conditions)
  
  # Общие гены
  shared_genes <- intersect(rownames(fetal_log), rownames(adult_log))
  fetal_log <- fetal_log[shared_genes, ]
  adult_log <- adult_log[shared_genes, ]
}

#==============================================================================
# ЧАСТЬ 1: Ко-экспрессия DLX5/MSX1/DSPP по всем fetal кластерам
#==============================================================================
cat("\n=== ЧАСТЬ 1: DLX5/MSX1/DSPP ПО ВСЕМ КЛАСТЕРАМ ===\n")

clusters <- c("DEM", "OB", "POB", "SOB", "SOBP")
coexpr_results <- data.frame()

for(cl in clusters) {
  cells <- which(cell_annot$cluster == cl)
  n_cells <- length(cells)
  
  dlx5_pos <- which(fetal_log["DLX5", cells] > 0)
  msx1_pos <- which(fetal_log["MSX1", cells] > 0)
  dspp_pos <- which(fetal_log["DSPP", cells] > 0)
  
  dlx5_msx1 <- length(intersect(dlx5_pos, msx1_pos))
  dlx5_dspp <- length(intersect(dlx5_pos, dspp_pos))
  msx1_dspp <- length(intersect(msx1_pos, dspp_pos))
  triple <- length(Reduce(intersect, list(dlx5_pos, msx1_pos, dspp_pos)))
  
  coexpr_results <- rbind(coexpr_results, data.frame(
    cluster = cl,
    n_cells = n_cells,
    DLX5 = length(dlx5_pos),
    MSX1 = length(msx1_pos),
    DSPP = length(dspp_pos),
    DLX5_MSX1 = dlx5_msx1,
    DLX5_DSPP = dlx5_dspp,
    MSX1_DSPP = msx1_dspp,
    TRIPLE = triple
  ))
}

print(coexpr_results, row.names = FALSE)

#==============================================================================
# ЧАСТЬ 2: Профиль DSPP+ клеток в каждом кластере
#==============================================================================
cat("\n=== ЧАСТЬ 2: ПРОФИЛЬ DSPP+ КЛЕТОК ===\n")

# Какие гены высокие в DSPP+ клетках POB?
pob_cells <- which(cell_annot$cluster == "POB")
pob_dspp_cells <- pob_cells[fetal_log["DSPP", pob_cells] > 0]
pob_non_dspp_cells <- pob_cells[fetal_log["DSPP", pob_cells] == 0]

cat("DSPP+ в POB:", length(pob_dspp_cells), "клеток\n")
cat("DSPP- в POB:", length(pob_non_dspp_cells), "клеток\n")

# Топ гены в DSPP+ vs DSPP-
genes_to_check <- c("DLX5", "DLX6", "MSX1", "MSX2", "PTN", "MDK", "BMP4",
                    "WNT10A", "AXIN2", "LEF1", "IGF1R", "RUNX2", "SP7",
                    "COL1A1", "COL1A2", "DMP1", "FOS", "JUN", "JUNB",
                    "NFKB1", "SOSTDC1", "SOX2", "BMPR1A", "BMPR1B", "BMPR2",
                    "TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP7")

dspp_profile <- data.frame()
for(g in genes_to_check) {
  if(!g %in% shared_genes) next
  
  dspp_vals <- fetal_log[g, pob_dspp_cells]
  non_dspp_vals <- fetal_log[g, pob_non_dspp_cells]
  
  dspp_profile <- rbind(dspp_profile, data.frame(
    gene = g,
    DSPP_pos_pct = round(100 * sum(dspp_vals > 0) / length(dspp_vals), 1),
    DSPP_pos_mean = round(mean(dspp_vals), 3),
    DSPP_neg_pct = round(100 * sum(non_dspp_vals > 0) / length(non_dspp_vals), 1),
    DSPP_neg_mean = round(mean(non_dspp_vals), 3),
    fold_change = round(mean(dspp_vals) / (mean(non_dspp_vals) + 0.001), 2)
  ))
}

dspp_profile <- dspp_profile[order(dspp_profile$fold_change, decreasing = TRUE), ]
cat("\nТоп гены в DSPP+ vs DSPP- (в POB):\n")
print(dspp_profile, row.names = FALSE)

write.csv(dspp_profile, "dspp_positive_profile.csv", row.names = FALSE)

#==============================================================================
# ЧАСТЬ 3: DSPP+ в OB (зрелые одонтобласты) — их профиль
#==============================================================================
cat("\n=== ЧАСТЬ 3: DSPP+ В OB (зрелые одонтобласты) ===\n")

ob_cells <- which(cell_annot$cluster == "OB")
ob_dspp_cells <- ob_cells[fetal_log["DSPP", ob_cells] > 0]
ob_non_dspp_cells <- ob_cells[fetal_log["DSPP", ob_cells] == 0]

cat("DSPP+ в OB:", length(ob_dspp_cells), "из", length(ob_cells), "\n")

ob_profile <- data.frame()
for(g in genes_to_check) {
  if(!g %in% shared_genes) next
  
  dspp_vals <- fetal_log[g, ob_dspp_cells]
  non_dspp_vals <- fetal_log[g, ob_non_dspp_cells]
  
  ob_profile <- rbind(ob_profile, data.frame(
    gene = g,
    OB_DSPP_pos_pct = round(100 * sum(dspp_vals > 0) / length(dspp_vals), 1),
    OB_DSPP_pos_mean = round(mean(dspp_vals), 3),
    OB_DSPP_neg_pct = round(100 * sum(non_dspp_vals > 0) / length(non_dspp_vals), 1),
    OB_DSPP_neg_mean = round(mean(non_dspp_vals), 3),
    fold_change = round(mean(dspp_vals) / (mean(non_dspp_vals) + 0.001), 2)
  ))
}

ob_profile <- ob_profile[order(ob_profile$fold_change, decreasing = TRUE), ]
cat("\nТоп гены в DSPP+ vs DSPP- (в OB):\n")
print(ob_profile, row.names = FALSE)

write.csv(ob_profile, "ob_dspp_positive_profile.csv", row.names = FALSE)

#==============================================================================
# ЧАСТЬ 4: 3 клетки DSPP+ в adult deep — где они?
#==============================================================================
cat("\n=== ЧАСТЬ 4: DSPP+ КЛЕТКИ В ADULT ===\n")

adult_dspp_cells <- which(adult_log["DSPP", ] > 0)
cat("Всего DSPP+ в adult:", length(adult_dspp_cells), "\n")

if(length(adult_dspp_cells) > 0) {
  cat("Условия этих клеток:\n")
  print(table(adult_conditions[adult_dspp_cells]))
  
  # Их профиль экспрессии
  cat("\nПрофиль DSPP+ клеток в adult:\n")
  adult_dspp_profile <- data.frame()
  for(g in genes_to_check) {
    if(!g %in% shared_genes) next
    vals <- adult_log[g, adult_dspp_cells]
    adult_dspp_profile <- rbind(adult_dspp_profile, data.frame(
      gene = g,
      pct_pos = round(100 * sum(vals > 0) / length(vals), 1),
      mean_expr = round(mean(vals), 3)
    ))
  }
  adult_dspp_profile <- adult_dspp_profile[order(adult_dspp_profile$mean_expr, decreasing = TRUE), ]
  print(adult_dspp_profile, row.names = FALSE)
}

#==============================================================================
# ЧАСТЬ 5: Ключевой вопрос — что объединяет DSPP+ клетки в fetal?
#==============================================================================
cat("\n=== ЧАСТЬ 5: ОБЪЕДИНЁННЫЙ ПРОФИЛЬ ВСЕХ DSPP+ КЛЕТОК В FETAL ===\n")

all_fetal_dspp <- which(fetal_log["DSPP", ] > 0)
all_fetal_non_dspp <- which(fetal_log["DSPP", ] == 0)
cat("Всего DSPP+ в fetal:", length(all_fetal_dspp), "\n")
cat("Распределение по кластерам:\n")
print(table(cell_annot$cluster[all_fetal_dspp]))

# Топ гены, специфичные для DSPP+ во всём fetal
fetal_dspp_signature <- data.frame()
for(g in genes_to_check) {
  if(!g %in% shared_genes) next
  pos_vals <- fetal_log[g, all_fetal_dspp]
  neg_vals <- fetal_log[g, all_fetal_non_dspp]
  fetal_dspp_signature <- rbind(fetal_dspp_signature, data.frame(
    gene = g,
    DSPP_pos_pct = round(100 * sum(pos_vals > 0) / length(pos_vals), 1),
    DSPP_pos_mean = round(mean(pos_vals), 3),
    DSPP_neg_pct = round(100 * sum(neg_vals > 0) / length(neg_vals), 1),
    DSPP_neg_mean = round(mean(neg_vals), 3),
    fold_change = round(mean(pos_vals) / (mean(neg_vals) + 0.001), 2)
  ))
}

fetal_dspp_signature <- fetal_dspp_signature[order(fetal_dspp_signature$fold_change, decreasing = TRUE), ]
cat("\nСигнатура DSPP+ клеток во всём fetal (по fold change):\n")
print(fetal_dspp_signature, row.names = FALSE)

write.csv(fetal_dspp_signature, "fetal_dspp_signature.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Сохранено:\n")
cat("  dspp_positive_profile.csv (POB: DSPP+ vs DSPP-)\n")
cat("  ob_dspp_positive_profile.csv (OB: DSPP+ vs DSPP-)\n")
cat("  fetal_dspp_signature.csv (все DSPP+ во fetal)\n")