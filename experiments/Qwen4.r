library(dplyr)
setwd("D:/!jenja/ARC/DATA3")

# Загружаем полные данные (не top-5000, а все 20293 гена)
data <- readRDS("GSE185222_processed.RDS")
expr_full <- data$expr  # 20293 генов x 6582 клеток

# Проверяем ключевые гены
genes_to_check <- c("SOX2", "DSPP", "DMP1", "SOSTDC1", "USAG1", 
                    "BMP4", "WNT10A", "AXIN2", "LEF1")

cat("=== ЭКСПРЕССИЯ КЛЮЧЕВЫХ ГЕНОВ ===\n\n")

for(gene in genes_to_check) {
  if(gene %in% rownames(expr_full)) {
    expr_gene <- expr_full[gene, ]
    n_cells <- sum(expr_gene > 0)
    pct_cells <- 100 * n_cells / length(expr_gene)
    mean_expr <- mean(expr_gene)
    
    cat(sprintf("%-10s: %5d клеток (%.2f%%), mean=%.3f\n", 
                gene, n_cells, pct_cells, mean_expr))
  } else {
    cat(sprintf("%-10s: НЕТ В ДАТАСЕТЕ\n", gene))
  }
}

# Проверяем по условиям
cat("\n=== ЭКСПРЕССИЯ ПО УСЛОВИЯМ ===\n")
conditions <- data$conditions

for(gene in c("SOX2", "DSPP", "SOSTDC1")) {
  if(gene %in% rownames(expr_full)) {
    cat(sprintf("\n%s:\n", gene))
    for(cond in c("sound", "enamel_caries", "deep_caries")) {
      cells_cond <- which(conditions == paste0(cond, "_"))
      if(length(cells_cond) > 0) {
        expr_cond <- expr_full[gene, cells_cond]
        n_pos <- sum(expr_cond > 0)
        pct_pos <- 100 * n_pos / length(expr_cond)
        cat(sprintf("  %-15s: %4d/%4d клеток (%.2f%%)\n", 
                    cond, n_pos, length(expr_cond), pct_pos))
      }
    }
  }
}