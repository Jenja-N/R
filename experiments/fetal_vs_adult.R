#==============================================================================
# FETAL vs ADULT TOOTH: ПРЯМОЕ СРАВНЕНИЕ
# GSE184749 (fetal 9-22w) vs GSE185222 (adult sound/enamel/deep caries)
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")

#==============================================================================
# ШАГ 0: Исследование структуры fetal данных
#==============================================================================
cat("=== ШАГ 0: СТРУКТУРА FETAL ДАННЫХ ===\n")

fetal_dir <- "GSE184749"
fetal_files <- list.files(fetal_dir, full.names = TRUE)
print(basename(fetal_files))

# Читаем аннотации (CSV, маленькие)
cell_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"))
gene_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"))

cat("\nКлеточные аннотации:\n")
cat("  Колонок:", ncol(cell_annot), "\n")
cat("  Строк:", nrow(cell_annot), "\n")
cat("  Имена:", paste(colnames(cell_annot), collapse=", "), "\n")
cat("  Первые строки:\n")
print(head(cell_annot, 3))

cat("\nГенные аннотации:\n")
cat("  Колонок:", ncol(gene_annot), "\n")
cat("  Строк:", nrow(gene_annot), "\n")
cat("  Имена:", paste(colnames(gene_annot), collapse=", "), "\n")
cat("  Первые строки:\n")
print(head(gene_annot, 3))

# Читаем MTX матрицу
cat("\nЗагружаю MTX матрицу (это займёт 1-2 минуты)...\n")
if(!require("Matrix")) install.packages("Matrix")
library(Matrix)

mtx <- readMM(gzfile("GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"))
cat("Размер MTX:", dim(mtx), "\n")

# Определяем ориентацию (гены x клетки или наоборот)
if(nrow(mtx) == nrow(gene_annot) & ncol(mtx) == nrow(cell_annot)) {
  cat("Ориентация: гены x клетки (стандартная)\n")
  rownames(mtx) <- gene_annot[,1]  # первая колонка = gene_id/symbol
  colnames(mtx) <- cell_annot[,1]  # первая колонка = cell_id
} else if(ncol(mtx) == nrow(gene_annot) & nrow(mtx) == nrow(cell_annot)) {
  cat("Ориентация: клетки x гены (транспонированная)\n")
  mtx <- t(mtx)
  rownames(mtx) <- gene_annot[,1]
  colnames(mtx) <- cell_annot[,1]
} else {
  stop("Размерности MTX не совпадают с аннотациями!")
}

cat("Финальная матрица:", nrow(mtx), "генов x", ncol(mtx), "клеток\n")

# Добавляем аннотации клеток
cell_annot$cell_id_check <- cell_annot[,1]
cell_annot <- cell_annot[match(colnames(mtx), cell_annot$cell_id_check), ]

#==============================================================================
# ШАГ 1: Загрузка adult данных
#==============================================================================
cat("\n=== ШАГ 1: ADULT ДАННЫЕ ===\n")
adult_data <- readRDS("GSE185222_processed.RDS")
adult_expr <- adult_data$expr
adult_cond <- adult_data$conditions

cat("Adult:", nrow(adult_expr), "генов x", ncol(adult_expr), "клеток\n")

#==============================================================================
# ШАГ 2: Нормализация fetal (log1p CPM)
#==============================================================================
cat("\n=== ШАГ 2: НОРМАЛИЗАЦИЯ ===\n")

# Fetal: CPM + log1p
fetal_lib_size <- Matrix::colSums(mtx)
fetal_cpm <- t(t(mtx) / fetal_lib_size) * 1e6
fetal_log <- log1p(fetal_cpm)
# Преобразуем в обычную матрицу для удобства
fetal_log <- as.matrix(fetal_log)

# Adult: уже log-normalized (из Seurat pipeline), оставляем как есть
adult_log <- adult_expr

cat("Fetal log range:", range(fetal_log), "\n")
cat("Adult log range:", range(adult_log), "\n")

#==============================================================================
# ШАГ 3: Общие гены
#==============================================================================
cat("\n=== ШАГ 3: ОБЩИЕ ГЕНЫ ===\n")
shared_genes <- intersect(rownames(fetal_log), rownames(adult_log))
cat("Общих генов:", length(shared_genes), "\n")

fetal_shared <- fetal_log[shared_genes, ]
adult_shared <- adult_log[shared_genes, ]

#==============================================================================
# ШАГ 4: Сравнение ключевых генов
#==============================================================================
cat("\n=== ШАГ 4: СРАВНЕНИЕ КЛЮЧЕВЫХ ГЕНОВ ===\n")

# Fetal кластеры
fetal_clusters <- cell_annot$cluster
cat("Fetal кластеры:", paste(unique(fetal_clusters), collapse=", "), "\n")

# Adult условия
adult_conditions <- gsub("_.*", "", adult_cond)
cat("Adult условия:", paste(unique(adult_conditions), collapse=", "), "\n")

# Ключевые гены
key_genes <- c("SOX2", "DSPP", "DMP1", "DLX5", "DLX6", "MSX1", "MSX2",
               "PTN", "MDK", "BMP4", "WNT10A", "AXIN2", "LEF1",
               "IGF1R", "RUNX2", "SP7", "COL1A1", "COL1A2",
               "FOS", "JUN", "JUNB", "NFKB1", "RELA",
               "SOSTDC1", "BMPR1A", "BMPR1B", "BMPR2")

# Функция: % позитивных клеток и средняя экспрессия
calc_stats <- function(expr_matrix, groups) {
  results <- data.frame()
  for(g in key_genes) {
    if(!g %in% rownames(expr_matrix)) {
      row <- data.frame(gene = g, group = unique(groups), 
                        pct_pos = NA, mean_expr = NA)
      results <- rbind(results, row)
      next
    }
    for(grp in unique(groups)) {
      cells <- which(groups == grp)
      vals <- expr_matrix[g, cells]
      pct_pos <- 100 * sum(vals > 0) / length(vals)
      mean_expr <- mean(vals)
      results <- rbind(results, data.frame(gene = g, group = grp, 
                                           pct_pos = round(pct_pos, 2), 
                                           mean_expr = round(mean_expr, 3)))
    }
  }
  results
}

fetal_stats <- calc_stats(fetal_shared, fetal_clusters)
adult_stats <- calc_stats(adult_shared, adult_conditions)

# Сводная таблица: широкая форма
library(reshape2)
fetal_wide <- dcast(fetal_stats, gene ~ group, value.var = "pct_pos")
adult_wide <- dcast(adult_stats, gene ~ group, value.var = "pct_pos")

cat("\n=== FETAL: % позитивных клеток по кластерам ===\n")
print(fetal_wide)

cat("\n=== ADULT: % позитивных клеток по условиям ===\n")
print(adult_wide)

# Сохраняем
write.csv(fetal_wide, "fetal_pct_pos_by_cluster.csv", row.names = FALSE)
write.csv(adult_wide, "adult_pct_pos_by_condition.csv", row.names = FALSE)

#==============================================================================
# ШАГ 5: Прямое сравнение: fetal OB vs adult deep_caries
#==============================================================================
cat("\n=== ШАГ 5: FETAL OB vs ADULT DEEP_CARIES ===\n")

# Fetal OB — зрелые одонтобласты (где DSPP = 34%)
# Adult deep_caries — где программа ломается

fetal_ob_cells <- which(fetal_clusters == "OB")
adult_deep_cells <- which(adult_conditions == "deep_caries")

cat("Fetal OB клеток:", length(fetal_ob_cells), "\n")
cat("Adult deep_caries клеток:", length(adult_deep_cells), "\n")

comparison <- data.frame()
for(g in key_genes) {
  if(!g %in% shared_genes) next
  
  fetal_vals <- fetal_shared[g, fetal_ob_cells]
  adult_vals <- adult_shared[g, adult_deep_cells]
  
  fetal_pct <- 100 * sum(fetal_vals > 0) / length(fetal_vals)
  adult_pct <- 100 * sum(adult_vals > 0) / length(adult_vals)
  fetal_mean <- mean(fetal_vals)
  adult_mean <- mean(adult_vals)
  
  comparison <- rbind(comparison, data.frame(
    gene = g,
    fetal_OB_pct = round(fetal_pct, 2),
    fetal_OB_mean = round(fetal_mean, 3),
    adult_deep_pct = round(adult_pct, 2),
    adult_deep_mean = round(adult_mean, 3),
    delta_pct = round(fetal_pct - adult_pct, 2),
    delta_mean = round(fetal_mean - adult_mean, 3)
  ))
}

# Сортируем по наибольшей потере в adult
comparison <- comparison[order(comparison$delta_pct, decreasing = TRUE), ]

cat("\n=== FETAL OB vs ADULT DEEP_CARIES (отсортировано по потере) ===\n")
print(comparison, row.names = FALSE)

write.csv(comparison, "fetal_OB_vs_adult_deep.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Сохранено:\n")
cat("  fetal_pct_pos_by_cluster.csv\n")
cat("  adult_pct_pos_by_condition.csv\n")
cat("  fetal_OB_vs_adult_deep.csv\n")