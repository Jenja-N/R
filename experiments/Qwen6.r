setwd("D:/!jenja/ARC/DATA3")

# 1. Убедимся, что данные загружены (если сессия была перезагружена)
if(!exists("expr_full") || !exists("conditions")) {
  data <- readRDS("GSE185222_processed.RDS")
  expr_full <- data$expr
  conditions <- data$conditions
}

# 2. Находим клетки, где ОБА гена экспрессируются (> 0)
# expr_full["SOSTDC1", ] > 0 возвращает логический вектор
# expr_full["BMP4", ] > 0 возвращает логический вектор
# Оператор & находит их пересечение
coexpressing_cells <- colnames(expr_full)[expr_full["SOSTDC1", ] > 0 & expr_full["BMP4", ] > 0]

# 3. Выводим результаты
cat("=== АНАЛИЗ КО-ЭКСПРЕССИИ SOSTDC1 И BMP4 ===\n")
cat("Всего клеток с ко-экспрессией:", length(coexpressing_cells), "\n\n")

if(length(coexpressing_cells) > 0) {
  cat("Имена клеток:\n")
  print(coexpressing_cells)
  
  cat("\nРаспределение по условиям:\n")
  # Извлекаем условия для этих конкретных клеток
  cell_conditions <- conditions[coexpressing_cells]
  # Очищаем названия условий (убираем "_GSM...")
  clean_conditions <- gsub("_.*", "", cell_conditions)
  print(table(clean_conditions))
} else {
  cat("ПЕРЕСЕЧЕНИЕ РАВНО 0.\n")
  cat("Это означает, что SOSTDC1 и BMP4 экспрессируются в РАЗНЫХ клетках.\n")
  cat("SOSTDC1 действует как паракринный (секретируемый) ингибитор BMP, а не аутокринный.\n")
}