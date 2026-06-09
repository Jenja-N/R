library(GENIE3)
library(dplyr)
setwd("D:/!jenja/ARC/DATA3")

# Загружаем уже обработанные данные
data <- readRDS("GSE185222_processed.RDS")
expr_top <- data$expr_top
conditions <- data$conditions

# TF список
key_tfs <- c(
  "DLX5", "DLX6", "DLX3", "MSX1", "MSX2", "TCF7", "RUNX2", "SOX9", "NFIC",
  "TWIST1", "TWIST2", "GLI3", "FOS", "JUN", "JUNB", "FOSB", "EGR1", "NFKB1",
  "NFKB2", "RELA", "REL", "STAT1", "STAT3", "NFKBIA", "MYC", "SMAD5", "SMAD9", "KLF4"
)

conditions_list <- c("sound", "enamel_caries", "deep_caries")
genie3_results <- list()

for(cond in conditions_list) {
  cat("\n=== Обработка условия:", cond, "===\n")
  
  cond_cols <- grep(paste0("^", cond, "_"), colnames(expr_top))
  if(length(cond_cols) == 0) {
    cat("  WARNING: Нет клеток для", cond, "\n")
    next
  }
  
  sub_expr <- expr_top[, cond_cols]
  cat("  Клеток:", ncol(sub_expr), "\n")
  
  # Убираем zero-variance
  sub_var <- apply(sub_expr, 1, var, na.rm = TRUE)
  sub_expr <- sub_expr[!is.na(sub_var) & sub_var > 0, ]
  
  tfs_sub <- intersect(key_tfs, rownames(sub_expr))
  cat("  TF доступно:", length(tfs_sub), "\n")
  
  if(length(tfs_sub) == 0) {
    cat("  WARNING: Нет TF, пропускаю\n")
    next
  }
  
  cat("  Запускаю GENIE3 (nTrees=500)...\n")
  set.seed(42)
  
  # ИСПРАВЛЕНО: убран nThreads
  genie3_out <- GENIE3(sub_expr, regulators = tfs_sub, nTrees = 500)
  
  link_list <- getLinkList(genie3_out, threshold = 0.005)
  link_list$condition <- cond
  
  cat("  Найдено связей:", nrow(link_list), "\n")
  genie3_results[[cond]] <- link_list
  
  saveRDS(genie3_results, "genie3_by_condition.RDS")
  cat("  Промежуточный результат сохранен\n")
}

cat("\n=== GENIE3 ЗАВЕРШЕН ===\n")