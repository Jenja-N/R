library(dplyr)
setwd("D:/!jenja/ARC/DATA3")

# Загружаем результаты GENIE3 (должен существовать после Step 3)
genie3_results <- readRDS("genie3_by_condition.RDS")

target_genes <- c("PTN", "COL1A1", "COL1A2", "MDK", "SOX2", "DSPP")
key_tfs <- c("DLX5", "DLX6", "MSX1", "MSX2", "SOX9", "FOS", "JUN", 
             "NFKB1", "RELA", "RUNX2", "SP7")

conditions_list <- c("sound", "enamel_caries", "deep_caries")

regulon_matrix <- data.frame()
for(tf in key_tfs) {
  for(tg in target_genes) {
    row <- data.frame(TF = tf, Target = tg, stringsAsFactors = FALSE)
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

cat("\n=== NF-kB SUPPRESSION TEST ===\n")
nfkb_test <- regulon_matrix %>%
  filter(TF %in% c("NFKB1", "RELA")) %>%
  select(TF, Target, sound, enamel_caries, deep_caries)
print(nfkb_test)