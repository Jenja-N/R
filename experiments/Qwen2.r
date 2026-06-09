#==============================================================================
# ПРОДОЛЖЕНИЕ АНАЛИЗА: STEP 2-7
# Использует сохраненные данные из GSE185222_processed.RDS
#==============================================================================

library(GENIE3)
library(dplyr)

setwd("D:/!jenja/ARC/DATA3")

# Загружаем уже обработанные данные
cat("=== ЗАГРУЗКА СОХРАНЕННЫХ ДАННЫХ ===\n")
data <- readRDS("GSE185222_processed.RDS")
expr_top <- data$expr_top
conditions <- data$conditions
cat("Загружено:", nrow(expr_top), "генов x", ncol(expr_top), "клеток\n")

#==============================================================================
# STEP 2: СПИСОК TRANSCRIPTION FACTORS
#==============================================================================
cat("\n=== STEP 2: TF CURATION ===\n")

key_tfs <- c(
  "DLX5", "DLX6", "DLX3", "MSX1", "MSX2", "PAX9", "AXIN2", 
  "LEF1", "TCF7", "RUNX2", "SP7", "SOX9", "SOX2", "BARX1", 
  "BARX2", "PITX2", "NFIC", "ETV1", "ETV4", "ETV5", "TWIST1",
  "TWIST2", "HAND1", "HAND2", "GLI1", "GLI2", "GLI3",
  "FOS", "JUN", "JUNB", "FOSB", "EGR1", "NFKB1", "NFKB2", 
  "RELA", "RELB", "REL", "STAT1", "STAT3", "NFKBIA",
  "MYC", "MYCN", "SMAD1", "SMAD5", "SMAD9", "KLF4", "KLF5",
  "SP1", "GATA1", "GATA2", "GATA3"
)

tfs_in_data <- intersect(key_tfs, rownames(expr_top))
cat("TF в анализе:", length(tfs_in_data), "\n")
cat("TF:", paste(tfs_in_data, collapse = ", "), "\n")

if(length(tfs_in_data) == 0) {
  stop("Не найдено ни одного TF в top вариабельных генах!")
}

#==============================================================================
# STEP 3: GENIE3 ПО УСЛОВИЯМ
#==============================================================================
cat("\n=== STEP 3: GENIE3 STRATIFIED RUNS ===\n")

conditions_list <- c("sound", "enamel_caries", "deep_caries")
genie3_results <- list()

for(cond in conditions_list) {
  cat("\nОбработка условия:", cond, "\n")
  
  # Фильтруем клетки
  cond_cols <- grep(paste0("^", cond, "_"), colnames(expr_top))
  if(length(cond_cols) == 0) {
    cat("  WARNING: Нет клеток для", cond, "\n")
    next
  }
  
  sub_expr <- expr_top[, cond_cols]
  cat("  Клеток:", ncol(sub_expr), "\n")
  
  # Убираем zero-variance гены в подвыборке
  sub_var <- apply(sub_expr, 1, var, na.rm = TRUE)
  sub_expr <- sub_expr[!is.na(sub_var) & sub_var > 0, ]
  
  # TF для этой подвыборки
  tfs_sub <- intersect(tfs_in_data, rownames(sub_expr))
  cat("  TF доступно:", length(tfs_sub), "\n")
  
  if(length(tfs_sub) == 0) {
    cat("  WARNING: Нет TF, пропускаю\n")
    next
  }
  
  cat("  Запускаю GENIE3 (nTrees=500, nThreads=4)...\n")
  cat("  (Это может занять 20-40 минут для каждого условия)\n")
  set.seed(42)
  
  # GENIE3
  genie3_out <- GENIE3(sub_expr, regulators = tfs_sub, 
                       nTrees = 500, nThreads = 4)
  
  # Извлекаем связи
  link_list <- getLinkList(genie3_out, threshold = 0.005)
  link_list$condition <- cond
  
  cat("  Найдено связей:", nrow(link_list), "\n")
  genie3_results[[cond]] <- link_list
  
  # Сохраняем промежуточные результаты
  saveRDS(genie3_results, "genie3_by_condition.RDS")
}

cat("\n=== GENIE3 ЗАВЕРШЕН ===\n")

#==============================================================================
# STEP 4: СРАВНЕНИЕ REGULONS
#==============================================================================
cat("\n=== STEP 4: СРАВНЕНИЕ REGULONS ===\n")

target_genes <- c("PTN", "COL1A1", "COL1A2", "MDK", "SOX2", "DSPP")
key_tfs <- c("DLX5", "DLX6", "MSX1", "MSX2", "SOX9", "FOS", "JUN", 
             "NFKB1", "RELA", "RUNX2", "SP7")

regulon_matrix <- data.frame()

for(tf in key_tfs) {
  for(tg in target_genes) {
    row <- data.frame(TF = tf, Target = tg)
    for(cond in conditions_list) {
      links <- genie3_results[[cond]]
      if(!is.null(links)) {
        w <- links$weight[links$regulatoryGene == tf & links$targetGene == tg]
        row[[cond]] <- ifelse(length(w) > 0, w[1], 0)
      } else {
        row[[cond]] <- NA
      }
    }
    regulon_matrix <- rbind(regulon_matrix, row)
  }
}

cat("\n=== МАТРИЦА ВЕСОВ REGULONS ===\n")
print(regulon_matrix, row.names = FALSE)
write.csv(regulon_matrix, "regulon_comparison.csv", row.names = FALSE)

#==============================================================================
# STEP 5: ТЕСТ ГИПОТЕЗЫ NF-κB
#==============================================================================
cat("\n=== STEP 5: ТЕСТ ГИПОТЕЗЫ INFLAMMAGING ===\n")

cat("Проверяем: NFKB1 подавляет репаративные гены?\n\n")

nfkb_test <- regulon_matrix %>%
  filter(TF %in% c("NFKB1", "RELA")) %>%
  select(TF, Target, sound, enamel_caries, deep_caries)

print(nfkb_test)

cat("\n=== ИНТЕРПРЕТАЦИЯ ===\n")
nfkb1_data <- nfkb_test %>% filter(TF == "NFKB1")
if(nrow(nfkb1_data) > 0) {
  for(i in 1:nrow(nfkb1_data)) {
    tg <- nfkb1_data$Target[i]
    s <- as.numeric(nfkb1_data$sound[i])
    e <- as.numeric(nfkb1_data$enamel_caries[i])
    d <- as.numeric(nfkb1_data$deep_caries[i])
    
    if(!is.na(d) && !is.na(s)) {
      if(d > s && d > e) {
        cat(sprintf("✓ NFKB1 -> %s: РОСТ при кариесе (%.3f -> %.3f -> %.3f)\n", 
                    tg, s, e, d))
      } else if(s > d && e > d) {
        cat(sprintf("✗ NFKB1 -> %s: ПАДЕНИЕ при кариесе (%.3f -> %.3f -> %.3f)\n", 
                    tg, s, e, d))
      } else {
        cat(sprintf("~ NFKB1 -> %s: Нет четкого тренда (%.3f -> %.3f -> %.3f)\n", 
                    tg, s, e, d))
      }
    }
  }
}

#==============================================================================
# STEP 6: MASTER REGULATORS
#==============================================================================
cat("\n=== STEP 6: MASTER REGULATORS ===\n")

module1_genes <- c("PTN", "COL1A1", "COL1A2", "VIM", "CLU", "FOS", "JUN", 
                   "JUNB", "HSP90AA1", "HSPA1A")
module1_in_data <- intersect(module1_genes, rownames(expr_top))

mr_results <- data.frame()

for(cond in conditions_list) {
  links <- genie3_results[[cond]]
  if(is.null(links)) next
  
  for(tf in key_tfs) {
    tf_targets <- links$targetGene[links$regulatoryGene == tf]
    overlap <- length(intersect(tf_targets, module1_in_data))
    total <- length(unique(links$targetGene))
    module1_total <- length(module1_in_data)
    
    enrichment <- ifelse(length(tf_targets) > 0 & total > 0,
                         (overlap / length(tf_targets)) / (module1_total / total),
                         0)
    
    mr_results <- rbind(mr_results, data.frame(
      Condition = cond,
      TF = tf,
      Total_targets = length(tf_targets),
      Module1_overlap = overlap,
      Enrichment = round(enrichment, 3)
    ))
  }
}

cat("\n=== TOP 5 MASTER REGULATORS ПО УСЛОВИЯМ ===\n")
mr_results %>%
  group_by(Condition) %>%
  arrange(desc(Enrichment)) %>%
  slice_head(n = 5) %>%
  print()

write.csv(mr_results, "master_regulators.csv", row.names = FALSE)

#==============================================================================
# STEP 7: ФИНАЛЬНАЯ СИНТЕЗА
#==============================================================================
cat("\n=== STEP 7: ФИНАЛЬНАЯ СИНТЕЗА ===\n")

cat("
=====================================================
РЕЗУЛЬТАТЫ GENIE3 АНАЛИЗА:
=====================================================

1. Проверьте файл regulon_comparison.csv:
   - Если DLX5/MSX1 -> PTN/COL1A1 падает от sound к deep_caries:
     => Репаративная ось разрушается при кариесе
   - Если NFKB1 -> PTN/COL1A1 растет от sound к deep_caries:
     => ПРЯМОЕ ДОКАЗАТЕЛЬСТВО inflammaging

2. Проверьте файл master_regulators.csv:
   - Какие TF наиболее специфично регулируют Module 1 (PTN/COL1A1)?
   - Меняется ли их активность при кариесе?

3. Фармакологические мишени:
   - Активация DLX5/MSX1 (CHIR99021, BMP7)
   - Блокада NF-kB (BAY 11-7082, Celastrol)
   - Прямая доставка PTN

=====================================================
АНАЛИЗ ЗАВЕРШЕН
Файлы: regulon_comparison.csv, master_regulators.csv
=====================================================
")

cat("\nГотово! Проверьте CSV файлы для детального анализа.\n")