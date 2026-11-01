# =============================================================================
# 01_curva_guia.R
#
# DATA CHALLENGE 07 -- PART 1: GUIDE CURVE METHOD
#
# Idea: fit ONE growth curve Hd = f(age) to ALL observations pooled together
# (ignoring which plot each point came from). That average curve is the
# "guide curve". Because the curves it generates are the guide curve scaled up
# or down, they all have the same shape (ANAMORPHIC):
#       Hd(t; S) = S * g(t) / g(Ii)      ->      S = Hd * g(Ii) / g(t)
# where g(t) is the shape part of the fitted model and Ii is the index age.
#
# HOW TO USE:
#   1. Edit SITIO_FILE, FUNCOES_FILE and OUT_DIR below.
#   2. Run: Rscript 01_curva_guia.R        Requires: readxl
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                          # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"              # <-- change
OUT_DIR      <- "C:/path/to/your/01_curva_guia/outputs"               # <-- change
TOLERANCIA_SYX <- 0.02   # another model must beat Chapman-Richards' Syx by >2% to replace it
# =============================================================================

if (!file.exists(SITIO_FILE))   stop("File not found: ", SITIO_FILE)
if (!file.exists(FUNCOES_FILE)) stop("File not found: ", FUNCOES_FILE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
source(FUNCOES_FILE)

d <- ler_sitio(SITIO_FILE)
cat("Observations:", nrow(d), "| plots:", length(unique(d$parcela)), "\n")

# -----------------------------------------------------------------------
# STEP 1 -- explore the data (measurements per plot, ranges, first plot)
# -----------------------------------------------------------------------
cat("\nMeasurements per plot:\n"); print(table(d$parcela))
cat("\nAge range:", round(range(d$idade), 2), "| Hd range:", range(d$Hd), "\n")

parcelas <- sort(unique(d$parcela))
cores <- rainbow(length(parcelas), s = 0.8, v = 0.8)
png(file.path(OUT_DIR, "exploracao_hd_idade.png"), width = 1200, height = 800, res = 140)
plot(d$idade, d$Hd, type = "n", xlab = "Age (years)", ylab = "Hd (m)",
     main = "Dominant height vs age, by plot")
for (i in seq_along(parcelas)) {
  g <- d[d$parcela == parcelas[i], ]
  lines(g$idade, g$Hd, type = "b", pch = 16, col = cores[i])
}
legend("topleft", legend = parcelas, col = cores, pch = 16, ncol = 3, cex = 0.7, bty = "n")
dev.off()

# -----------------------------------------------------------------------
# STEP 2 -- fit the candidate guide-curve models on the pooled data
# -----------------------------------------------------------------------
# (a) Chapman-Richards  Hd = b0*(1-exp(-b1*t))^b2     (nonlinear, nls)
fit_cr <- ajustar_nls(
  Hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = d,
  starts = list(list(b0 = 34, b1 = 0.15, b2 = 0.85),
                list(b0 = 50, b1 = 0.10, b2 = 1.00),
                list(b0 = 25, b1 = 0.40, b2 = 1.50)),
  lower = c(5, 0.001, 0.1), upper = c(500, 5, 10))
p_cr <- coef(fit_cr)
b1_cr <- unname(p_cr["b1"]); b2_cr <- unname(p_cr["b2"])

# (b) Backman  ln(Hd) = b0 + b1*ln(t) + b2*ln(t)^2      (linear after the log)
fit_bk <- lm(log(Hd) ~ log(idade) + I(log(idade)^2), data = d)
p_bk <- unname(coef(fit_bk))

# (c) Schumacher  ln(Hd) = b0 + b1/t                    (linear after the log)
fit_sc <- lm(log(Hd) ~ I(1 / idade), data = d)
p_sc <- unname(coef(fit_sc))

# every metric is computed on the ORIGINAL Hd scale so the three are comparable
pred_cr <- predict(fit_cr)
pred_bk <- exp(predict(fit_bk))
pred_sc <- exp(predict(fit_sc))

cand <- rbind(
  metricas_ajuste(d$Hd, pred_cr, 3, "chapman_richards"),
  metricas_ajuste(d$Hd, pred_bk, 3, "backman"),
  metricas_ajuste(d$Hd, pred_sc, 2, "schumacher"))
cat("\n=== Candidate guide-curve models ===\n"); print(cand)

# -----------------------------------------------------------------------
# STEP 3 -- select the model
# -----------------------------------------------------------------------
# The three fits are statistically almost identical (see Syx). The assignment
# requires ONE base model reused in Parts 1, 2, 3 and 5, so Chapman-Richards
# is kept unless another candidate is clearly better (Syx lower by more than
# TOLERANCIA_SYX). Reasons for the default: its asymptote b0 has a biological
# meaning, and the algebraic-difference method (Part 2) needs it.
syx_cr <- cand$syx[cand$modelo == "chapman_richards"]
outros <- cand[cand$modelo != "chapman_richards", ]
melhor_outro <- outros[which.min(outros$syx), ]
modelo_sel <- if (melhor_outro$syx < (1 - TOLERANCIA_SYX) * syx_cr) melhor_outro$modelo else "chapman_richards"
cat("\nSelected guide-curve model:", modelo_sel, "\n")

# shape function g(t) = Hd(t; S = 1 at Ii) of the selected model.
# Hd(t; S) = S * g(t);  and  S = Hd / g(t).
forma_g <- function(t) {
  if (modelo_sel == "chapman_richards") {
    ((1 - exp(-b1_cr * t)) / (1 - exp(-b1_cr * IDADE_INDICE)))^b2_cr
  } else if (modelo_sel == "backman") {
    exp(p_bk[2] * (log(t) - log(IDADE_INDICE)) +
        p_bk[3] * (log(t)^2 - log(IDADE_INDICE)^2))
  } else {
    exp(p_sc[2] * (1 / t - 1 / IDADE_INDICE))
  }
}
hd_guia <- function(idade, S) S * forma_g(idade)
S_guia  <- function(idade, hd) hd / forma_g(idade)

pred_sel <- switch(modelo_sel, chapman_richards = pred_cr, backman = pred_bk, schumacher = pred_sc)
write.csv(cand, file.path(OUT_DIR, "ajuste_curva_guia.csv"), row.names = FALSE)
grafico_residuos(d$Hd, pred_sel, d$idade, file.path(OUT_DIR, "residuos_curva_guia.png"),
                 "Guide curve")

# -----------------------------------------------------------------------
# STEP 4 -- family of site-index curves
# -----------------------------------------------------------------------
curvas <- montar_curvas(hd_guia)
write.csv(curvas, file.path(OUT_DIR, "curvas_curva_guia.csv"), row.names = FALSE)
grafico_feixe(curvas, d, file.path(OUT_DIR, "feixe_curva_guia.png"),
              paste0("Guide curve - site index curves (Ii = ", IDADE_INDICE, " yr)"))

# -----------------------------------------------------------------------
# STEP 5 -- classify the 9 plots
# -----------------------------------------------------------------------
# Every measurement gives its own S (Hd projected to the index age). Both are
# reported: (1) the measurement closest to Ii, and (2) the mean over ALL the
# plot's measurements. The mean is used for the class because it uses all the
# information and is less sensitive to a single noisy measurement.
d$S_obs <- S_guia(d$idade, d$Hd)
classif <- do.call(rbind, lapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  i <- which.min(abs(g$idade - IDADE_INDICE))
  data.frame(parcela = pl, n_medicoes = nrow(g),
             idade_mais_proxima = g$idade[i],
             S_medicao_mais_proxima = g$S_obs[i],
             S_media_medicoes = mean(g$S_obs),
             classe_S = classe_mais_proxima(mean(g$S_obs)))
}))
cat("\n=== Plot classification (guide curve) ===\n"); print(classif)
write.csv(classif, file.path(OUT_DIR, "indice_local_parcelas_curva_guia.csv"), row.names = FALSE)

cat("\nPART 1 FINISHED.\n")
