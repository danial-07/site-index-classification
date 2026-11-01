# =============================================================================
# 02_diferenca_algebrica.R
#
# DATA CHALLENGE 07 -- PART 2: ALGEBRAIC DIFFERENCE METHOD
#
# Idea: take the base model (Chapman-Richards, Hd = b0*(1-exp(-b1*t))^b2),
# write it at an earlier age I1 and a later age I2 of the SAME plot, and
# eliminate the parameter that carries the site quality. What remains is an
# equation that predicts Hd2 from (I1, Hd1, I2) alone -- no S needed as input.
# It is fitted to PAIRS of measurements from the same plot.
#
# Which parameter is eliminated decides the curve family:
#   Variant A: isolate b0 (asymptote varies with site)  -> ANAMORPHIC curves
#       b0 = Hd1 / (1-exp(-b1*I1))^b2
#       Hd2 = Hd1 * [(1-exp(-b1*I2)) / (1-exp(-b1*I1))]^b2        (fit b1, b2)
#   Variant B: isolate b2 (shape varies with site)      -> POLYMORPHIC curves
#       b2 = ln(Hd1/b0) / ln(1-exp(-b1*I1))
#       Hd2 = b0 * (1-exp(-b1*I2))^[ ln(Hd1/b0) / ln(1-exp(-b1*I1)) ]   (fit b0, b1)
# Both are fitted; the one with the smaller residual sum of squares is used.
# Index curves: set I1 = Ii (index age) and Hd1 = S.
#
# HOW TO USE: edit the 3 paths below, run  Rscript 02_diferenca_algebrica.R
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                          # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"              # <-- change
OUT_DIR      <- "C:/path/to/your/02_diferenca_algebrica/outputs"      # <-- change
# =============================================================================

if (!file.exists(SITIO_FILE))   stop("File not found: ", SITIO_FILE)
if (!file.exists(FUNCOES_FILE)) stop("File not found: ", FUNCOES_FILE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
source(FUNCOES_FILE)

d <- ler_sitio(SITIO_FILE)
parcelas <- sort(unique(d$parcela))

# -----------------------------------------------------------------------
# STEP 1 -- build every pair (earlier, later) within each plot
# -----------------------------------------------------------------------
# A plot with m measurements gives m*(m-1)/2 pairs (4 -> 6, 5 -> 10).
# Pairs from one plot share measurements, so they are not independent; this
# does not bias the estimates but makes standard errors optimistic (README).
pares <- do.call(rbind, lapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  g <- g[order(g$idade), ]
  idx <- t(combn(nrow(g), 2))              # all (i < j) index pairs
  data.frame(parcela = pl,
             t1 = g$idade[idx[, 1]], h1 = g$Hd[idx[, 1]],
             t2 = g$idade[idx[, 2]], h2 = g$Hd[idx[, 2]])
}))
cat("Pairs of measurements:", nrow(pares), "\n")

# -----------------------------------------------------------------------
# STEP 2 -- fit both difference equations by nonlinear regression
# -----------------------------------------------------------------------
fit_A <- ajustar_nls(
  h2 ~ h1 * ((1 - exp(-b1 * t2)) / (1 - exp(-b1 * t1)))^b2, data = pares,
  starts = list(list(b1 = 0.15, b2 = 0.85), list(b1 = 0.30, b2 = 1.10),
                list(b1 = 0.05, b2 = 0.50)),
  lower = c(0.001, 0.05), upper = c(5, 10))

# b0 must exceed every Hd1 (otherwise ln(Hd1/b0) is undefined), so its lower
# bound is just above the largest observed Hd1.
lim_b0 <- max(pares$h1) * 1.01
fit_B <- ajustar_nls(
  h2 ~ b0 * (1 - exp(-b1 * t2))^(log(h1 / b0) / log(1 - exp(-b1 * t1))), data = pares,
  starts = list(list(b0 = 34, b1 = 0.15), list(b0 = 50, b1 = 0.10),
                list(b0 = 80, b1 = 0.05)),
  lower = c(lim_b0, 0.001), upper = c(1000, 5))

mA <- metricas_ajuste(pares$h2, predict(fit_A), 2, "A_b0_isolado_anamorfica")
mB <- metricas_ajuste(pares$h2, predict(fit_B), 2, "B_b2_isolado_polimorfica")
cand <- rbind(mA, mB)
cand$parametros <- c(paste(names(coef(fit_A)), round(coef(fit_A), 5), collapse = "; "),
                     paste(names(coef(fit_B)), round(coef(fit_B), 5), collapse = "; "))
cat("\n=== Difference-equation variants (metrics on Hd2) ===\n"); print(cand)

# -----------------------------------------------------------------------
# STEP 3 -- select the variant (same number of parameters -> smallest SSE)
# -----------------------------------------------------------------------
usa_B <- mB$sse < mA$sse
cat("\nSelected variant:", if (usa_B) "B (polymorphic)" else "A (anamorphic)", "\n")
fit_sel <- if (usa_B) fit_B else fit_A
write.csv(cand, file.path(OUT_DIR, "ajuste_diferenca_algebrica.csv"), row.names = FALSE)
grafico_residuos(pares$h2, predict(fit_sel), pares$t2,
                 file.path(OUT_DIR, "residuos_diferenca_algebrica.png"),
                 "Algebraic difference")

# -----------------------------------------------------------------------
# STEP 4 -- site-index curves and S of any observation
# -----------------------------------------------------------------------
if (usa_B) {
  b0 <- unname(coef(fit_B)["b0"]); b1 <- unname(coef(fit_B)["b1"])
  hd_ada <- function(idade, S) b0 * (1 - exp(-b1 * idade))^(log(S / b0) / log(1 - exp(-b1 * IDADE_INDICE)))
  S_ada  <- function(idade, hd) b0 * (1 - exp(-b1 * IDADE_INDICE))^(log(hd / b0) / log(1 - exp(-b1 * idade)))
} else {
  b1 <- unname(coef(fit_A)["b1"]); b2 <- unname(coef(fit_A)["b2"])
  hd_ada <- function(idade, S) S * ((1 - exp(-b1 * idade)) / (1 - exp(-b1 * IDADE_INDICE)))^b2
  S_ada  <- function(idade, hd) hd * ((1 - exp(-b1 * IDADE_INDICE)) / (1 - exp(-b1 * idade)))^b2
}

curvas <- montar_curvas(hd_ada)
write.csv(curvas, file.path(OUT_DIR, "curvas_diferenca_algebrica.csv"), row.names = FALSE)
grafico_feixe(curvas, d, file.path(OUT_DIR, "feixe_diferenca_algebrica.png"),
              paste0("Algebraic difference - site index curves (Ii = ", IDADE_INDICE, " yr)"))

# -----------------------------------------------------------------------
# STEP 5 -- anamorphic or polymorphic?
# -----------------------------------------------------------------------
# Anamorphic curves are proportional: Hd(S=21)/Hd(S=15) is the same at every
# age. If the ratio changes with age the curves have different shapes.
cv <- montar_curvas(hd_ada, classes = c(min(CLASSES_S), max(CLASSES_S)), idades = c(2.5, 6.5))
razao_jovem <- cv$Hd_estimada_m[cv$S == max(CLASSES_S) & cv$idade == 2.5] /
               cv$Hd_estimada_m[cv$S == min(CLASSES_S) & cv$idade == 2.5]
razao_velha <- cv$Hd_estimada_m[cv$S == max(CLASSES_S) & cv$idade == 6.5] /
               cv$Hd_estimada_m[cv$S == min(CLASSES_S) & cv$idade == 6.5]
cat(sprintf("\nHd(S=%d)/Hd(S=%d): %.3f at 2.5 yr, %.3f at 6.5 yr -> %s\n",
            max(CLASSES_S), min(CLASSES_S), razao_jovem, razao_velha,
            if (abs(razao_jovem - razao_velha) > 0.01) "POLYMORPHIC" else "anamorphic"))

# -----------------------------------------------------------------------
# STEP 6 -- classify the 9 plots
# -----------------------------------------------------------------------
d$S_obs <- S_ada(d$idade, d$Hd)
classif <- do.call(rbind, lapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  i <- which.min(abs(g$idade - IDADE_INDICE))
  data.frame(parcela = pl, n_medicoes = nrow(g),
             idade_mais_proxima = g$idade[i],
             S_medicao_mais_proxima = g$S_obs[i],
             S_media_medicoes = mean(g$S_obs),
             classe_S = classe_mais_proxima(mean(g$S_obs)))
}))
cat("\n=== Plot classification (algebraic difference) ===\n"); print(classif)
write.csv(classif, file.path(OUT_DIR, "indice_local_parcelas_diferenca_algebrica.csv"), row.names = FALSE)

cat("\nPART 2 FINISHED.\n")
