#==============================================================================
# КРИТИЧЕСКИЙ ТЕСТ: WNT-сигналинг (AXIN2/LEF1) в DSPP+ клетках adult
# Есть ли WNT-путь в клетках, пытающихся экспрессировать DSPP?
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")
library(Matrix)

# Загрузка данных (если не загружены)
if(!exists("fetal_log") || !exists("adult_log")) {
  cat("Загружаю данные...\n")
  source("QWEN9.R")
}

#==============================================================================
# ЧАСТЬ 1: Ко-экспрессия AXIN2/LEF1/DSPP в adult enamel_caries
#==============================================================================
cat("\n=== ЧАСТЬ 1: ADULT ENAMEL_CARIES ===\n")

adult_enamel_cells <- which(adult_conditions == "enamel")
cat("Клеток в enamel_caries:", length(adult_enamel_cells), "\n")

axin2_pos <- which(adult_log["AXIN2", adult_enamel_cells] > 0)
lef1_pos <- which(adult_log["LEF1", adult_enamel_cells] > 0)
dspp_pos <- which(adult_log["DSPP", adult_enamel_cells] > 0)

cat("AXIN2+ клеток:", length(axin2_pos), "\n")
cat("LEF1+ клеток:", length(lef1_pos), "\n")
cat("DSPP+ клеток:", length(dspp_pos), "\n")

axin2_lef1 <- intersect(axin2_pos, lef1_pos)
axin2_dspp <- intersect(axin2_pos, dspp_pos)
lef1_dspp <- intersect(lef1_pos, dspp_pos)
triple <- Reduce(intersect, list(axin2_pos, lef1_pos, dspp_pos))

cat("\nAXIN2+LEF1+:", length(axin2_lef1), "\n")
cat("AXIN2+DSPP+:", length(axin2_dspp), "\n")
cat("LEF1+DSPP+:", length(lef1_dspp), "\n")
cat("AXIN2+LEF1+DSPP+:", length(triple), "\n")

# Профиль DSPP+ клеток в enamel_caries
if(length(dspp_pos) > 0) {
  cat("\n=== ПРОФИЛЬ DSPP+ КЛЕТОК В ENAMEL ===\n")
  dspp_cells_enamel <- adult_enamel_cells[dspp_pos]
  
  wnt_genes <- c("AXIN2", "LEF1", "WNT10A", "WNT5A", "CTNNB1", "TCF7", "TCF7L1", "TCF7L2")
  available_wnt <- intersect(wnt_genes, shared_genes)
  
  cat("WNT-гены в DSPP+ клетках enamel:\n")
  for(g in available_wnt) {
    vals <- adult_log[g, dspp_cells_enamel]
    pct_pos <- 100 * sum(vals > 0) / length(vals)
    mean_expr <- mean(vals)
    cat(sprintf("  %-10s: %5.1f%% pos, mean=%.3f\n", g, pct_pos, mean_expr))
  }
}

#==============================================================================
# ЧАСТЬ 2: Ко-экспрессия AXIN2/LEF1/DSPP в adult deep_caries
#==============================================================================
cat("\n=== ЧАСТЬ 2: ADULT DEEP_CARIES ===\n")

adult_deep_cells <- which(adult_conditions == "deep")
cat("Клеток в deep_caries:", length(adult_deep_cells), "\n")

axin2_pos_deep <- which(adult_log["AXIN2", adult_deep_cells] > 0)
lef1_pos_deep <- which(adult_log["LEF1", adult_deep_cells] > 0)
dspp_pos_deep <- which(adult_log["DSPP", adult_deep_cells] > 0)

cat("AXIN2+ клеток:", length(axin2_pos_deep), "\n")
cat("LEF1+ клеток:", length(lef1_pos_deep), "\n")
cat("DSPP+ клеток:", length(dspp_pos_deep), "\n")

axin2_lef1_deep <- intersect(axin2_pos_deep, lef1_pos_deep)
axin2_dspp_deep <- intersect(axin2_pos_deep, dspp_pos_deep)
lef1_dspp_deep <- intersect(lef1_pos_deep, dspp_pos_deep)
triple_deep <- Reduce(intersect, list(axin2_pos_deep, lef1_pos_deep, dspp_pos_deep))

cat("\nAXIN2+LEF1+:", length(axin2_lef1_deep), "\n")
cat("AXIN2+DSPP+:", length(axin2_dspp_deep), "\n")
cat("LEF1+DSPP+:", length(lef1_dspp_deep), "\n")
cat("AXIN2+LEF1+DSPP+:", length(triple_deep), "\n")

#==============================================================================
# ЧАСТЬ 3: Сравнение с fetal POB и OB
#==============================================================================
cat("\n=== ЧАСТЬ 3: FETAL POB vs OB ===\n")

# POB
pob_cells <- which(cell_annot$cluster == "POB")
pob_axin2 <- which(fetal_log["AXIN2", pob_cells] > 0)
pob_lef1 <- which(fetal_log["LEF1", pob_cells] > 0)
pob_dspp <- which(fetal_log["DSPP", pob_cells] > 0)

pob_triple <- Reduce(intersect, list(pob_axin2, pob_lef1, pob_dspp))

cat("Fetal POB:\n")
cat("  AXIN2+:", length(pob_axin2), "\n")
cat("  LEF1+:", length(pob_lef1), "\n")
cat("  DSPP+:", length(pob_dspp), "\n")
cat("  AXIN2+LEF1+DSPP+:", length(pob_triple), "\n")

# OB
ob_cells <- which(cell_annot$cluster == "OB")
ob_axin2 <- which(fetal_log["AXIN2", ob_cells] > 0)
ob_lef1 <- which(fetal_log["LEF1", ob_cells] > 0)
ob_dspp <- which(fetal_log["DSPP", ob_cells] > 0)

ob_triple <- Reduce(intersect, list(ob_axin2, ob_lef1, ob_dspp))

cat("\nFetal OB:\n")
cat("  AXIN2+:", length(ob_axin2), "\n")
cat("  LEF1+:", length(ob_lef1), "\n")
cat("  DSPP+:", length(ob_dspp), "\n")
cat("  AXIN2+LEF1+DSPP+:", length(ob_triple), "\n")

#==============================================================================
# ЧАСТЬ 4: Сводная таблица
#==============================================================================
cat("\n=== СВОДНАЯ ТАБЛИЦА ===\n")

summary_table <- data.frame(
  condition = c("Fetal POB", "Fetal OB", "Adult enamel", "Adult deep"),
  n_cells = c(length(pob_cells), length(ob_cells), length(adult_enamel_cells), length(adult_deep_cells)),
  AXIN2_pos = c(length(pob_axin2), length(ob_axin2), length(axin2_pos), length(axin2_pos_deep)),
  LEF1_pos = c(length(pob_lef1), length(ob_lef1), length(lef1_pos), length(lef1_pos_deep)),
  DSPP_pos = c(length(pob_dspp), length(ob_dspp), length(dspp_pos), length(dspp_pos_deep)),
  TRIPLE = c(length(pob_triple), length(ob_triple), length(triple), length(triple_deep))
)

summary_table$AXIN2_pct <- round(100 * summary_table$AXIN2_pos / summary_table$n_cells, 2)
summary_table$LEF1_pct <- round(100 * summary_table$LEF1_pos / summary_table$n_cells, 2)
summary_table$DSPP_pct <- round(100 * summary_table$DSPP_pos / summary_table$n_cells, 2)

print(summary_table, row.names = FALSE)

write.csv(summary_table, "wnt_signaling_summary.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Сохранено: wnt_signaling_summary.csv\n")