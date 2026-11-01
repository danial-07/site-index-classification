# =============================================================================
# 06_comparacao_metodos.R   (extra script, for the consolidated report)
#
# Reads the five curvas_*.csv files and the five indice_local_parcelas_*.csv
# files produced by scripts 01-05, draws ONE comparison plot of the five
# families of curves, and tabulates (a) the shape test (is the curve family
# anamorphic?) and (b) the S class given to each plot by each method.
#
# HOW TO USE: run scripts 01-05 first, then edit the paths and run this one.
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                          # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"              # <-- change
BASE_DIR     <- "C:/path/to/your/desafio-07-classificacao-sitio"      # <-- the repository folder
OUT_DIR      <- file.path(BASE_DIR, "06_comparacao_metodos/outputs")
# =============================================================================

source(FUNCOES_FILE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
d <- ler_sitio(SITIO_FILE)

metodos <- list(
  curva_guia            = "01_curva_guia",
  diferenca_algebrica   = "02_diferenca_algebrica",
  atribuicao_preliminar = "03_atribuicao_preliminar",
  hammer                = "04_hammer",
  predicao_parametros   = "05_predicao_parametros")
rotulos <- c("Guide curve", "Algebraic difference", "Preliminary assignment",
             "Hammer", "Parameter prediction")
cores <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
ltys  <- c(1, 2, 4, 3, 5)

curvas <- lapply(names(metodos), function(m)
  read.csv(file.path(BASE_DIR, metodos[[m]], "outputs", paste0("curvas_", m, ".csv"))))
names(curvas) <- names(metodos)

png(file.path(OUT_DIR, "comparacao_metodos.png"), width = 1500, height = 1000, res = 150)
par(mar = c(4.5, 4.5, 3, 1))
plot(d$idade, d$Hd, pch = 16, col = "grey80", xlim = c(2.3, 6.8), ylim = c(10, 28),
     xlab = "Age (years)", ylab = "Dominant height Hd (m)",
     main = paste0("Site index curves, five methods (S = ",
                   paste(CLASSES_S, collapse = ", "), " m at Ii = ", IDADE_INDICE, " yr)"))
for (i in seq_along(curvas)) {
  for (s in CLASSES_S) {
    g <- curvas[[i]][curvas[[i]]$S == s, ]
    lines(g$idade, g$Hd_estimada_m, col = cores[i], lty = ltys[i], lwd = 2)
  }
}
legend("topleft", legend = rotulos, col = cores, lty = ltys, lwd = 2, bty = "n", cex = 0.9)
dev.off()

# (a) shape test: Hd(S_max)/Hd(S_min) at a young and an old age
forma <- do.call(rbind, lapply(seq_along(curvas), function(i) {
  cv <- curvas[[i]]; smin <- min(CLASSES_S); smax <- max(CLASSES_S)
  r <- function(a) cv$Hd_estimada_m[cv$S == smax & cv$idade == a] /
                   cv$Hd_estimada_m[cv$S == smin & cv$idade == a]
  data.frame(metodo = rotulos[i], razao_idade_2.5 = r(2.5), razao_idade_6.5 = r(6.5),
             tipo = if (abs(r(2.5) - r(6.5)) > 0.01) "polymorphic" else "anamorphic")
}))
cat("\n=== Shape test: Hd(S_max)/Hd(S_min) ===\n"); print(forma)
write.csv(forma, file.path(OUT_DIR, "teste_forma_curvas.csv"), row.names = FALSE)

# (b) S class given to each plot by each method
col_S <- c(curva_guia = "classe_S", diferenca_algebrica = "classe_S",
           atribuicao_preliminar = "classe_reclassificada", hammer = "classe_S",
           predicao_parametros = "classe_S_equacao_final")
tab <- NULL
for (m in names(metodos)) {
  x <- read.csv(file.path(BASE_DIR, metodos[[m]], "outputs", paste0("indice_local_parcelas_", m, ".csv")))
  x <- x[, c("parcela", col_S[[m]])]; names(x)[2] <- m
  tab <- if (is.null(tab)) x else merge(tab, x, by = "parcela")
}
cat("\n=== S class per plot and method ===\n"); print(tab)
write.csv(tab, file.path(OUT_DIR, "classe_S_por_metodo.csv"), row.names = FALSE)
cat("\nCOMPARISON FINISHED.\n")
