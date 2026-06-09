#==============================================================================
# АНАЛИЗ EMBRYO: Проверка гипотезы TGFB2/IGF1R vs стресс
#==============================================================================

library(Seurat)
library(dplyr)

setwd("D:/!jenja/ARC/DATA3/embryonic tooth germ scRNA and ST")
load("embryo.RData")

# 1. Смотрим типы клеток
cat("=== ТИПЫ КЛЕТОК В EMBRYO ===\n")
print(table(embryo$celltype))

# 2. Извлекаем экспрессию с типами клеток
target_genes <- c("DSPP", "DMP1", "TGFB2", "TGFB3", "BMP2", "BMP4", "BMP7", 
                  "IGF1R", "RUNX2", "FOS", "JUN", "JUNB", "DLX5", "MSX1", "MSX2")
available_genes <- intersect(target_genes, rownames(embryo))

embryo_data <- FetchData(embryo, vars = c("celltype", "orig.ident", available_genes))

# 3. Считаем % позитивных клеток по типам
cat("\n=== % ПОЗИТИВНЫХ КЛЕТОК ПО ТИПАМ ===\n")
results <- data.frame()

for(ct in unique(embryo_data$celltype)) {
  cells <- which(embryo_data$celltype == ct)
  n_cells <- length(cells)
  
  for(g in available_genes) {
    vals <- embryo_data[cells, g]
    pct_pos <- 100 * sum(vals > 0) / n_cells
    mean_expr <- mean(vals)
    
    results <- rbind(results, data.frame(
      celltype = ct,
      gene = g,
      n_cells = n_cells,
      pct_pos = round(pct_pos, 2),
      mean_expr = round(mean_expr, 3)
    ))
  }
}

# Показываем сводку по ключевым генам
cat("\n=== КЛЮЧЕВЫЕ ГЕНЫ ПО ТИПАМ КЛЕТОК ===\n")
key_genes <- c("TGFB2", "IGF1R", "RUNX2", "FOS", "JUN", "JUNB", "DSPP", "BMP2", "BMP4")

for(g in key_genes) {
  cat(sprintf("\n%s:\n", g))
  subset <- results[results$gene == g, c("celltype", "pct_pos", "mean_expr")]
  subset <- subset[order(subset$pct_pos, decreasing = TRUE), ]
  print(subset, row.names = FALSE)
}

# 4. Сравнение: эмбрион vs adult deep_caries
cat("\n=== СРАВНЕНИЕ: EMBRYO vs ADULT DEEP ===\n")

# Загружаем adult данные
setwd("D:/!jenja/ARC/DATA3")
adult_data <- readRDS("GSE185222_processed.RDS")
adult_log <- adult_data$expr
adult_conditions <- gsub("_.*", "", adult_data$conditions)
adult_deep_cells <- which(adult_conditions == "deep")

# Считаем для adult deep
adult_results <- data.frame()
for(g in key_genes) {
  if(g %in% rownames(adult_log)) {
    vals <- adult_log[g, adult_deep_cells]
    pct_pos <- 100 * sum(vals > 0) / length(vals)
    mean_expr <- mean(vals)
    adult_results <- rbind(adult_results, data.frame(
      celltype = "adult_deep",
      gene = g,
      n_cells = length(adult_deep_cells),
      pct_pos = round(pct_pos, 2),
      mean_expr = round(mean_expr, 3)
    ))
  }
}

# Сравниваем: берем все клетки эмбриона vs adult deep
cat("\n=== ВСЕ КЛЕТКИ EMBRYO vs ADULT DEEP ===\n")
comparison <- data.frame()

for(g in key_genes) {
  # Embryo: все клетки
  embryo_vals <- embryo_data[[g]]
  embryo_pct <- 100 * sum(embryo_vals > 0) / length(embryo_vals)
  embryo_mean <- mean(embryo_vals)
  
  # Adult deep
  if(g %in% rownames(adult_log)) {
    adult_vals <- adult_log[g, adult_deep_cells]
    adult_pct <- 100 * sum(adult_vals > 0) / length(adult_vals)
    adult_mean <- mean(adult_vals)
    
    comparison <- rbind(comparison, data.frame(
      gene = g,
      embryo_pct = round(embryo_pct, 2),
      embryo_mean = round(embryo_mean, 3),
      adult_deep_pct = round(adult_pct, 2),
      adult_deep_mean = round(adult_mean, 3),
      delta_pct = round(embryo_pct - adult_pct, 2)
    ))
  }
}

comparison <- comparison[order(comparison$delta_pct, decreasing = TRUE), ]
print(comparison, row.names = FALSE)

write.csv(comparison, "embryo_vs_adult_deep_comparison.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Сохранено: embryo_vs_adult_deep_comparison.csv\n")