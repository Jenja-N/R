library(dplyr)
setwd("D:/!jenja/ARC/DATA3")

# Загружаем данные
data <- readRDS("GSE185222_processed.RDS")
expr_full <- data$expr
conditions <- data$conditions

# Проверяем по условиям
cat("=== ЭКСПРЕССИЯ ПО УСЛОВИЯМ ===\n")

genes_to_check <- c("SOX2", "DSPP", "DMP1", "SOSTDC1", "BMP4", "WNT10A")

for(gene in genes_to_check) {
  if(gene %in% rownames(expr_full)) {
    cat(sprintf("\n%s:\n", gene))
    
    for(cond in c("sound", "enamel_caries", "deep_caries")) {
      # Правильный синтаксис: ищем колонки, начинающиеся с cond_
      cells_cond <- grep(paste0("^", cond, "_"), colnames(expr_full))
      
      if(length(cells_cond) > 0) {
        expr_cond <- expr_full[gene, cells_cond]
        n_pos <- sum(expr_cond > 0)
        pct_pos <- 100 * n_pos / length(expr_cond)
        mean_expr <- mean(expr_cond)
        
        cat(sprintf("  %-15s: %4d/%4d клеток (%.2f%%), mean=%.3f\n", 
                    cond, n_pos, length(expr_cond), pct_pos, mean_expr))
      }
    }
  }
}

# Дополнительно: проверяем, в каких клетках экспрессируется SOX2
cat("\n=== SOX2+ КЛЕТКИ ПО УСЛОВИЯМ ===\n")
sox2_expr <- expr_full["SOX2", ]
sox2_pos <- which(sox2_expr > 0)

if(length(sox2_pos) > 0) {
  sox2_conditions <- conditions[sox2_pos]
  sox2_cond_clean <- gsub("_.*", "", sox2_conditions)
  
  cat("SOX2+ клетки по условиям:\n")
  print(table(sox2_cond_clean))
} else {
  cat("SOX2+ клеток не найдено\n")
}