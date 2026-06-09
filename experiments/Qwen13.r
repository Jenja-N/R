#==============================================================================
# ПРОВЕРКА ГИПОТЕЗЫ: TGFB/BMP vs СТРЕСС (Fetal OB vs Adult deep)
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")
if(!require("Matrix")) install.packages("Matrix")
library(Matrix)

# 1. Загрузка данных (если еще не в памяти)
if(!exists("fetal_log") || !exists("adult_log") || !exists("cell_annot")) {
  cat("Загрузка данных...\n")
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
  
  shared_genes <- intersect(rownames(fetal_log), rownames(adult_log))
  fetal_log <- fetal_log[shared_genes, ]
  adult_log <- adult_log[shared_genes, ]
}

# 2. Определение целевых генов
target_genes <- c(
  "TGFB1", "TGFB2", "TGFB3",  # TGF-beta pathway
  "BMP2", "BMP4", "BMP7",     # BMP pathway
  "IGF1R", "RUNX2",           # Differentiation competence
  "FOS", "JUN", "JUNB"        # Stress/AP-1 (negative regulators of differentiation)
)

# Фильтруем только те, что есть в данных
available_genes <- intersect(target_genes, shared_genes)
cat("Анализируемые гены:", paste(available_genes, collapse=", "), "\n\n")

# 3. Выбор клеток
fetal_ob_cells <- which(cell_annot$cluster == "OB")
adult_deep_cells <- which(adult_conditions == "deep")

cat("Fetal OB клеток:", length(fetal_ob_cells), "\n")
cat("Adult deep клеток:", length(adult_deep_cells), "\n\n")

# 4. Расчет статистики
results <- data.frame()

for(g in available_genes) {
  # Fetal OB
  fetal_vals <- fetal_log[g, fetal_ob_cells]
  fetal_pct <- 100 * sum(fetal_vals > 0) / length(fetal_vals)
  fetal_mean <- mean(fetal_vals)
  
  # Adult deep
  adult_vals <- adult_log[g, adult_deep_cells]
  adult_pct <- 100 * sum(adult_vals > 0) / length(adult_vals)
  adult_mean <- mean(adult_vals)
  
  results <- rbind(results, data.frame(
    gene = g,
    fetal_OB_pct = round(fetal_pct, 2),
    fetal_OB_mean = round(fetal_mean, 3),
    adult_deep_pct = round(adult_pct, 2),
    adult_deep_mean = round(adult_mean, 3),
    delta_pct = round(fetal_pct - adult_pct, 2),
    delta_mean = round(fetal_mean - adult_mean, 3)
  ))
}

# 5. Сортировка: сначала гены, которые ТЕРЯЮТСЯ в adult (положительный delta)
# затем гены, которые ПРИОБРЕТАЮТСЯ (отрицательный delta, т.е. стресс)
results <- results[order(results$delta_pct, decreasing = TRUE), ]

cat("=== РЕЗУЛЬТАТЫ СРАВНЕНИЯ (Fetal OB vs Adult deep) ===\n")
print(results, row.names = FALSE)

# 6. Сохранение
write.csv(results, "tgfb_bmp_stress_comparison.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Результат сохранен в: tgfb_bmp_stress_comparison.csv\n")
cat("\nИнтерпретация:\n")
cat("- Положительный delta_pct: гены, необходимые для дифференцировки (теряются в adult).\n")
cat("- Отрицательный delta_pct: гены стресса, блокирующие дифференцировку (доминируют в adult).\n")