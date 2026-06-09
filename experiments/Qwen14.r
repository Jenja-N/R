setwd("D:/!jenja/ARC/DATA3/embryonic tooth germ scRNA and ST")

cat("=== ЗАГРУЗКА embryo.RData ===\n")
load("embryo.RData")

cat("\nОбъекты в workspace:\n")
print(ls())

# Проверяем тип объекта
for(obj_name in ls()) {
  obj <- get(obj_name)
  cat(sprintf("\n%s: класс = %s\n", obj_name, class(obj)[1]))
  
  if(inherits(obj, "Seurat")) {
    cat("  Это Seurat объект\n")
    cat("  Клеток:", ncol(obj), "\n")
    cat("  Генов:", nrow(obj), "\n")
    if("ident" %in% names(obj@meta.data)) {
      cat("  Кластеров:", length(unique(obj@meta.data$ident)), "\n")
      print(table(obj@meta.data$ident))
    }
    if("orig.ident" %in% names(obj@meta.data)) {
      cat("  Образцов:", length(unique(obj@meta.data$orig.ident)), "\n")
      print(table(obj@meta.data$orig.ident))
    }
  }
}