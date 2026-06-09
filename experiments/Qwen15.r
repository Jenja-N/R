#==============================================================================
# ЗАГРУЗКА И АНАЛИЗ EMBRYO.RDATA (Human embryonic tooth germ)
#==============================================================================

# 1. Установка недостающих пакетов (выполнится только если их нет)
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  cat("Установка SeuratObject...\n")
  install.packages("SeuratObject")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  cat("Установка Seurat...\n")
  install.packages("Seurat")
}

library(Seurat)
library(dplyr)

setwd("D:/!jenja/ARC/DATA3/embryonic tooth germ scRNA and ST")

# 2. Загрузка объекта
cat("Загрузка embryo.RData...\n")
load("embryo.RData")

# 3. Базовая инспекция объекта
cat("\n=== ИНСПЕКЦИЯ ОБЪЕКТА EMBRYO ===\n")
print(embryo)

# Посмотрим, какие колонки есть в метаданных, чтобы найти кластеры/типы клеток
cat("\nКолонки метаданных (первые 10):\n")
print(head(colnames(embryo@meta.data), 10))

# Автоматический поиск колонки с кластерами или типами клеток
cluster_col <- grep("cluster|cell_type|ident|annotation", colnames(embryo@meta.data), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(cluster_col)) {
  cat(sprintf("\nРаспределение клеток по '%s':\n", cluster_col))
  print(table(embryo@meta.data[[cluster_col]]))
} else {
  cat("\nНе удалось автоматически найти колонку с кластерами. Проверьте colnames(embryo@meta.data)\n")
  cluster_col <- colnames(embryo@meta.data)[1] # fallback на первую колонку
}

# 4. Извлечение экспрессии ключевых генов для сравнения с adult deep
# Фокус на генах, которые мы выявили как критичные: TGFB2, IGF1R, RUNX2 (позитивные) 
# и FOS/JUN/JUNB (стресс, негативные), а также маркеры одонтогенеза.
target_genes <- c("DSPP", "DMP1", "TGFB2", "TGFB3", "BMP2", "BMP4", "BMP7", 
                  "IGF1R", "RUNX2", "FOS", "JUN", "JUNB", "DLX5", "MSX1", "MSX2")

# Фильтруем только те гены, которые реально есть в этом объекте
available_genes <- intersect(target_genes, rownames(embryo))
cat("\nДоступные целевые гены в объекте embryo:", paste(available_genes, collapse = ", "), "\n")

# Извлекаем данные (FetchData безопасно достает экспрессию и присоединяет метаданные)
cat("Извлечение данных экспрессии (это может занять 10-20 секунд)...\n")
embryo_expr <- FetchData(embryo, vars = c(cluster_col, available_genes))

# Сохраняем результат для кросс-сравнения
write.csv(embryo_expr, "../embryo_key_genes_expression.csv", row.names = FALSE)

cat("\n=== ГОТОВО ===\n")
cat("Данные сохранены в: D:/!jenja/ARC/DATA3/embryo_key_genes_expression.csv\n")
cat("Теперь мы можем напрямую сравнить эмбриональные паттерны этих генов с adult deep_caries.\n")