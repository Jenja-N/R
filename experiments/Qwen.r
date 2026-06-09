#==============================================================================
# STEP 0: РАСПАКОВКА И ЗАГРУЗКА ДАННЫХ (ИСПРАВЛЕНО)
#==============================================================================
cat("=== STEP 0: РАСПАКОВКА ДАННЫХ ===\n")

# Распаковка (если еще не сделано)
extract_dir <- "GSE185222_extracted"
if(!dir.exists(extract_dir)) dir.create(extract_dir)
tar_file <- "GSE185222_RAW/GSE185222_RAW.tar"
if(file.exists(tar_file)) {
  untar(tar_file, exdir = extract_dir)
}

# Поиск CSV файлов
csv_files <- list.files(extract_dir, pattern = "\\.csv\\.gz$", 
                        recursive = TRUE, full.names = TRUE)
cat("Найдено CSV файлов:", length(csv_files), "\n")

# Функция определения условия
get_condition <- function(filename) {
  if(grepl("DTP_01", filename)) return("sound")
  if(grepl("DTP_02", filename)) return("enamel_caries")
  if(grepl("DTP_04|DTP_05", filename)) return("deep_caries")
  return("unknown")
}

# Загружаем каждый файл отдельно, сохраняя gene_id как обычную колонку
sample_list <- list()
condition_annotations <- c()

for(f in csv_files) {
  cat("  Загружаю:", basename(f), "\n")
  cond <- get_condition(basename(f))
  
  # Читаем БЕЗ row.names, чтобы gene_id был обычной колонкой
  df <- read.csv(gzfile(f), stringsAsFactors = FALSE)
  
  # Первая колонка всегда gene_id (проверяем)
  first_col <- colnames(df)[1]
  if(first_col != "gene_id" && first_col != "Gene" && first_col != "gene") {
    # Переименовываем первую колонку в gene_id
    colnames(df)[1] <- "gene_id"
  } else {
    colnames(df)[1] <- "gene_id"
  }
  
  # Добавляем префикс условия ко всем остальным колонкам (клеткам)
  other_cols <- colnames(df)[-1]
  colnames(df)[-1] <- paste0(cond, "_", other_cols)
  
  sample_list[[basename(f)]] <- df
  condition_annotations <- c(condition_annotations, rep(cond, ncol(df) - 1))
}

# Объединяем через merge по gene_id (без проблем с Row.names)
cat("Объединяю матрицы...\n")
expr_df <- sample_list[[1]]
for(i in 2:length(sample_list)) {
  expr_df <- merge(expr_df, sample_list[[i]], by = "gene_id", all = TRUE)
}

# Убираем дубликаты генов (если есть)
expr_df <- expr_df[!duplicated(expr_df$gene_id), ]

# Переводим в матрицу
rownames(expr_df) <- expr_df$gene_id
expr_df$gene_id <- NULL

# Принудительно преобразуем в numeric (убираем любые строки)
expr_matrix <- as.matrix(expr_df)
mode(expr_matrix) <- "numeric"

# Заполняем NA нулями
expr_matrix[is.na(expr_matrix)] <- 0

cat("Итоговая матрица:", nrow(expr_matrix), "генов x", ncol(expr_matrix), "клеток\n")

# Проверка размерности (должно быть ~20293 x ~6582)
if(nrow(expr_matrix) > 30000) {
  cat("ВНИМАНИЕ: Очень много генов. Фильтруем по дисперсии...\n")
}

#==============================================================================
# STEP 1: ФИЛЬТРАЦИЯ ГЕНОВ (ИСПРАВЛЕНО)
#==============================================================================
cat("\n=== STEP 1: ФИЛЬТРАЦИЯ ГЕНОВ ===\n")

# Считаем дисперсию с na.rm=TRUE
gene_var <- apply(expr_matrix, 1, var, na.rm = TRUE)

# Убираем NA и гены с нулевой дисперсией
valid_genes <- !is.na(gene_var) & gene_var > 0
expr_filtered <- expr_matrix[valid_genes, ]
cat("Генов после фильтрации:", nrow(expr_filtered), "\n")

# Top 5000 вариабельных генов (с защитой от NA)
n_top <- min(5000, sum(valid_genes))
top_var_idx <- order(gene_var[valid_genes], decreasing = TRUE)[1:n_top]
expr_top <- expr_filtered[top_var_idx, ]
cat("Top вариабельных генов:", nrow(expr_top), "\n")

# Проверка целевых генов
targets <- c("PTN", "MDK", "COL1A1", "COL1A2", "DLX5", "MSX1", "NFKB1", "SOX2")
for(t in targets) {
  present <- t %in% rownames(expr_top)
  var_t <- ifelse(present, var(expr_top[t, ], na.rm = TRUE), NA)
  cat(sprintf("  %s: in_top=%s, var=%.4f\n", t, present, ifelse(is.na(var_t), 0, var_t)))
}

# Сохраняем для повторных запусков
saveRDS(list(expr = expr_matrix, conditions = condition_annotations, 
             expr_top = expr_top), "GSE185222_processed.RDS")
cat("Данные сохранены в GSE185222_processed.RDS для повторных запусков\n")