#==============================================================================
# FETAL vs ADULT TOOTH: ПРЯМОЕ СРАВНЕНИЕ (ИСПРАВЛЕННАЯ ВЕРСИЯ)
# GSE184749 (fetal) vs GSE185222 (adult)
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")
if(!require("Matrix")) install.packages("Matrix")
if(!require("reshape2")) install.packages("reshape2")
library(Matrix)
library(reshape2)

#==============================================================================
# ШАГ 0: Загрузка и исправление имен генов в FETAL данных
#==============================================================================
cat("=== ШАГ 0: ЗАГРУЗКА FETAL ДАННЫХ ===\n")

cell_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_cell_annotations.csv.gz"))
gene_annot <- read.csv(gzfile("GSE184749/GSE184749_dental_mesenchyme_gene_annotations.csv.gz"))

cat("Загружаю MTX матрицу...\n")
mtx <- readMM(gzfile("GSE184749/GSE184749_dental_mesenchyme_counts.mtx.gz"))

# КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ 1: Используем gene_short_name вместо Ensembl ID
rownames(mtx) <- gene_annot$gene_short_name
colnames(mtx) <- cell_annot$cell_id

# Убираем дубликаты символов генов (оставляем первое вхождение)
mtx <- mtx[!duplicated(rownames(mtx)), ]

cat("Fetal матрица:", nrow(mtx), "генов x", ncol(mtx), "клеток\n")

#==============================================================================
# ШАГ 1: Загрузка ADULT данных
#==============================================================================
cat("\n=== ШАГ 1: ЗАГРУЗКА ADULT ДАННЫХ ===\n")
adult_data <- readRDS("GSE185222_processed.RDS")
adult_expr <- adult_data$expr
adult_cond <- adult_data$conditions

# Приводим названия условий к короткому виду (как в предыдущем выводе)
adult_conditions <- gsub("_.*", "", adult_cond)
cat("Adult условия:", paste(unique(adult_conditions), collapse=", "), "\n")

#==============================================================================
# ШАГ 2: Нахождение ОБЩИХ генов (до конвертации в dense, чтобы сэкономить память)
#==============================================================================
cat("\n=== ШАГ 2: ОБЩИЕ ГЕНЫ ===\n")
shared_genes <- intersect(rownames(mtx), rownames(adult_expr))
cat("Общих генов найдено:", length(shared_genes), "\n")

if(length(shared_genes) < 1000) {
  stop("Слишком мало общих генов! Проверьте аннотации.")
}

# Оставляем только общие гены в обеих матрицах
fetal_shared_sparse <- mtx[shared_genes, ]
adult_shared <- adult_expr[shared_genes, ]

#==============================================================================
# ШАГ 3: Нормализация (теперь безопасно, так как матрица меньше)
#==============================================================================
cat("\n=== ШАГ 3: НОРМАЛИЗАЦИЯ ===\n")

# Fetal: CPM + log1p (конвертируем в dense только для общих генов)
cat("Нормализация fetal данных (может занять 30-60 сек)...\n")
fetal_lib_size <- Matrix::colSums(fetal_shared_sparse)
fetal_cpm <- t(t(fetal_shared_sparse) / fetal_lib_size) * 1e6
fetal_log <- as.matrix(log1p(fetal_cpm))

# Adult: уже log-normalized, просто берем нужные гены
adult_log <- adult_shared

cat("Fetal log range:", round(range(fetal_log), 2), "\n")
cat("Adult log range:", round(range(adult_log), 2), "\n")

#==============================================================================
# ШАГ 4: Сравнение ключевых генов
#==============================================================================
cat("\n=== ШАГ 4: СРАВНЕНИЕ КЛЮЧЕВЫХ ГЕНОВ ===\n")

fetal_clusters <- cell_annot$cluster

# Функция расчета статистики
calc_stats <- function(expr_matrix, groups) {
  results <- data.frame()
  for(g in key_genes) {
    if(!g %in% rownames(expr_matrix)) {
      for(grp in unique(groups)) {
        results <- rbind(results, data.frame(gene = g, group = grp, pct_pos = NA, mean_expr = NA))
      }
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

key_genes <- c("SOX2", "DSPP", "DMP1", "DLX5", "DLX6", "MSX1", "MSX2",
               "PTN", "MDK", "BMP4", "WNT10A", "AXIN2", "LEF1",
               "IGF1R", "RUNX2", "SP7", "COL1A1", "COL1A2",
               "FOS", "JUN", "JUNB", "NFKB1", "RELA", "SOSTDC1", 
               "BMPR1A", "BMPR1B", "BMPR2")

fetal_stats <- calc_stats(fetal_log, fetal_clusters)
adult_stats <- calc_stats(adult_log, adult_conditions)

fetal_wide <- dcast(fetal_stats, gene ~ group, value.var = "pct_pos")
adult_wide <- dcast(adult_stats, gene ~ group, value.var = "pct_pos")

cat("\n=== FETAL: % позитивных клеток по кластерам ===\n")
print(fetal_wide)

cat("\n=== ADULT: % позитивных клеток по условиям ===\n")
print(adult_wide)

write.csv(fetal_wide, "fetal_pct_pos_by_cluster.csv", row.names = FALSE)
write.csv(adult_wide, "adult_pct_pos_by_condition.csv", row.names = FALSE)

#==============================================================================
# ШАГ 5: Прямое сравнение: fetal OB vs adult deep
#==============================================================================
cat("\n=== ШАГ 5: FETAL OB vs ADULT DEEP ===\n")

fetal_ob_cells <- which(fetal_clusters == "OB")
# КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ 2: ищем "deep", а не "deep_caries"
adult_deep_cells <- which(adult_conditions == "deep") 

cat("Fetal OB клеток:", length(fetal_ob_cells), "\n")
cat("Adult deep клеток:", length(adult_deep_cells), "\n")

comparison <- data.frame()
for(g in key_genes) {
  if(!g %in% shared_genes) next
  
  fetal_vals <- fetal_log[g, fetal_ob_cells]
  adult_vals <- adult_log[g, adult_deep_cells]
  
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

# Сортируем по наибольшей потере в adult (delta_pct)
comparison <- comparison[order(comparison$delta_pct, decreasing = TRUE), ]

cat("\n=== FETAL OB vs ADULT DEEP (отсортировано по потере экспрессии) ===\n")
print(comparison, row.names = FALSE)

write.csv(comparison, "fetal_OB_vs_adult_deep.csv", row.names = FALSE)

cat("\n=== АНАЛИЗ ЗАВЕРШЕН УСПЕШНО ===\n")
cat("Проверьте файлы:\n")
cat("1. fetal_pct_pos_by_cluster.csv\n")
cat("2. adult_pct_pos_by_condition.csv\n")
cat("3. fetal_OB_vs_adult_deep.csv (ГЛАВНЫЙ РЕЗУЛЬТАТ)\n")