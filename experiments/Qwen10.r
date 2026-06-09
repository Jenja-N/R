#==============================================================================
# FETAL POB vs ADULT DEEP: Почему высокие DLX5/MSX1 не включают DSPP?
#==============================================================================

setwd("D:/!jenja/ARC/DATA3")
library(Matrix)

# Загружаем данные (если еще не загружены)
if(!exists("fetal_log") || !exists("adult_log")) {
  source("QWEN9.R")  # или fetal_vs_adult_FIXED.R
}

# Fetal POB (pre-odontoblasts) - клетки на пути к дифференцировке
fetal_pob_cells <- which(cell_annot$cluster == "POB")
cat("Fetal POB клеток:", length(fetal_pob_cells), "\n")

# Adult deep_caries
adult_deep_cells <- which(adult_conditions == "deep")
cat("Adult deep клеток:", length(adult_deep_cells), "\n")

# Ключевые гены для проверки гипотезы
key_genes <- c("DSPP", "DMP1", "DLX5", "DLX6", "MSX1", "MSX2",
               "PTN", "MDK", "BMP4", "WNT10A", "IGF1R", "RUNX2",
               "COL1A1", "COL1A2", "FOS", "JUN", "JUNB", "NFKB1",
               "SOSTDC1", "SOX2", "AXIN2", "LEF1")

# Сравниваем POB vs deep
comparison <- data.frame()
for(g in key_genes) {
  if(!g %in% shared_genes) next
  
  fetal_vals <- fetal_log[g, fetal_pob_cells]
  adult_vals <- adult_log[g, adult_deep_cells]
  
  fetal_pct <- 100 * sum(fetal_vals > 0) / length(fetal_vals)
  adult_pct <- 100 * sum(adult_vals > 0) / length(adult_vals)
  fetal_mean <- mean(fetal_vals)
  adult_mean <- mean(adult_vals)
  
  comparison <- rbind(comparison, data.frame(
    gene = g,
    fetal_POB_pct = round(fetal_pct, 2),
    fetal_POB_mean = round(fetal_mean, 3),
    adult_deep_pct = round(adult_pct, 2),
    adult_deep_mean = round(adult_mean, 3),
    delta_pct = round(fetal_pct - adult_pct, 2),
    delta_mean = round(fetal_mean - adult_mean, 3)
  ))
}

# Сортируем по delta_pct
comparison <- comparison[order(comparison$delta_pct, decreasing = TRUE), ]

cat("\n=== FETAL POB vs ADULT DEEP ===\n")
print(comparison, row.names = FALSE)

write.csv(comparison, "fetal_POB_vs_adult_deep.csv", row.names = FALSE)

# Ко-экспрессия: DLX5+MSX1+DSPP+ в POB
cat("\n=== КО-ЭКСПРЕССИЯ В FETAL POB ===\n")
pob_dlx5 <- which(fetal_log["DLX5", fetal_pob_cells] > 0)
pob_msx1 <- which(fetal_log["MSX1", fetal_pob_cells] > 0)
pob_dspp <- which(fetal_log["DSPP", fetal_pob_cells] > 0)

cat("DLX5+ клеток в POB:", length(pob_dlx5), "\n")
cat("MSX1+ клеток в POB:", length(pob_msx1), "\n")
cat("DSPP+ клеток в POB:", length(pob_dspp), "\n")

# Пересечения
dlx5_msx1 <- intersect(pob_dlx5, pob_msx1)
dlx5_dspp <- intersect(pob_dlx5, pob_dspp)
msx1_dspp <- intersect(pob_msx1, pob_dspp)
triple_pos <- intersect(dlx5_msx1, pob_dspp)

cat("\nDLX5+MSX1+:", length(dlx5_msx1), "\n")
cat("DLX5+DSPP+:", length(dlx5_dspp), "\n")
cat("MSX1+DSPP+:", length(msx1_dspp), "\n")
cat("DLX5+MSX1+DSPP+:", length(triple_pos), "\n")

cat("\n=== ГОТОВО ===\n")
cat("Сохранено: fetal_POB_vs_adult_deep.csv\n")