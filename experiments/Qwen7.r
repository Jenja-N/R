# Проверить экспрессию BMP-рецепторов
receptors <- c("BMPR1A", "BMPR1B", "BMPR2")
for(rec in receptors) {
  if(rec %in% rownames(expr_full)) {
    n_pos <- sum(expr_full[rec, ] > 0)
    pct_pos <- 100 * n_pos / ncol(expr_full)
    cat(sprintf("%s: %d клеток (%.2f%%)\n", rec, n_pos, pct_pos))
  }
}