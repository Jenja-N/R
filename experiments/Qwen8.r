# Средняя экспрессия по образцам
sostdc1_mean <- tapply(expr_full["SOSTDC1", ], conditions, mean)
bmp4_mean <- tapply(expr_full["BMP4", ], conditions, mean)
msx1_mean <- tapply(expr_full["MSX1", ], conditions, mean)

cat("SOSTDC1:", sostdc1_mean, "\n")
cat("BMP4:", bmp4_mean, "\n")
cat("MSX1:", msx1_mean, "\n")