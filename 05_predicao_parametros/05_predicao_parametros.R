# =============================================================================
# 05_predicao_parametros.R
#
# DATA CHALLENGE 07 -- PART 5: PARAMETER PREDICTION
#
# Idea: (1) fit the growth model to EACH plot separately; (2) give each plot a
# site index S = its own fitted curve at the index age; (3) regress the model
# parameters on S (b0 = a0 + a1*S, b2 = c0 + c1*S); (4) substitute those
# relations back into the model -> ONE equation Hd = f(age, S).
#
# IMPORTANT PRACTICAL DECISION (documented in README and report):
# each plot has only 4-5 measurements, and Chapman-Richards has 3 parameters,
# so a free 3-parameter fit per plot is unstable (1-2 degrees of freedom).
# The rate b1 is therefore FIXED at the value from the pooled fit (same as the
# guide curve) and b0 and b2 are fitted per plot. The free 3-parameter fits are
# still computed and saved as columns *_livre, only to show the instability.
#
# HOW TO USE: edit the 3 paths below, run  Rscript 05_predicao_parametros.R
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                          # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"              # <-- change
OUT_DIR      <- "C:/path/to/your/05_predicao_parametros/outputs"      # <-- change
# =============================================================================

if (!file.exists(SITIO_FILE))   stop("File not found: ", SITIO_FILE)
if (!file.exists(FUNCOES_FILE)) stop("File not found: ", FUNCOES_FILE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
source(FUNCOES_FILE)

d <- ler_sitio(SITIO_FILE)
parcelas <- sort(unique(d$parcela))

# -----------------------------------------------------------------------
# STEP 1 -- pooled fit (only to obtain the common rate b1)
# -----------------------------------------------------------------------
fit_pool <- ajustar_nls(
  Hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = d,
  starts = list(list(b0 = 34, b1 = 0.15, b2 = 0.85), list(b0 = 50, b1 = 0.10, b2 = 1.00),
                list(b0 = 25, b1 = 0.40, b2 = 1.50)),
  lower = c(5, 0.001, 0.1), upper = c(500, 5, 10))
B1_FIXO <- unname(coef(fit_pool)["b1"])
cat("Pooled fit:", paste(names(coef(fit_pool)), round(coef(fit_pool), 4), collapse = "  "), "\n")
cat("b1 fixed at", round(B1_FIXO, 5), "for the per-plot fits\n\n")

# -----------------------------------------------------------------------
# STEP 2 -- per-plot fits (b0, b2 free; b1 fixed) + diagnostic free fits
# -----------------------------------------------------------------------
por_parcela <- do.call(rbind, lapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  fit <- ajustar_nls(
    Hd ~ b0 * (1 - exp(-B1_FIXO * idade))^b2, data = g,
    starts = list(list(b0 = 34, b2 = 0.85), list(b0 = 30, b2 = 0.80), list(b0 = 40, b2 = 1.00)),
    lower = c(5, 0.1), upper = c(500, 10))
  b0 <- unname(coef(fit)["b0"]); b2 <- unname(coef(fit)["b2"])
  # diagnostic: all three parameters free (may fail -> NA)
  livre <- tryCatch(ajustar_nls(
    Hd ~ b0 * (1 - exp(-b1 * idade))^b2, data = g,
    starts = list(list(b0 = 34, b1 = 0.15, b2 = 0.85), list(b0 = 60, b1 = 0.08, b2 = 1.00),
                  list(b0 = 25, b1 = 0.30, b2 = 1.20), list(b0 = 20, b1 = 0.50, b2 = 1.80)),
    lower = c(5, 0.001, 0.1), upper = c(500, 5, 10)),
    error = function(e) NULL)
  cl <- if (is.null(livre)) c(NA, NA, NA) else unname(coef(livre)[c("b0", "b1", "b2")])
  data.frame(parcela = pl, n_medicoes = nrow(g), b0 = b0, b1 = B1_FIXO, b2 = b2,
             S_estimado = cr(IDADE_INDICE, b0, B1_FIXO, b2), sse = deviance(fit),
             b0_livre = cl[1], b1_livre = cl[2], b2_livre = cl[3],
             sse_livre = if (is.null(livre)) NA else deviance(livre))
}))
cat("=== Parameters per plot ===\n"); print(por_parcela)
write.csv(por_parcela, file.path(OUT_DIR, "parametros_por_parcela.csv"), row.names = FALSE)

# -----------------------------------------------------------------------
# STEP 3 -- relate the parameters to S (simple linear regressions)
# -----------------------------------------------------------------------
reg_b0 <- lm(b0 ~ S_estimado, data = por_parcela)
reg_b2 <- lm(b2 ~ S_estimado, data = por_parcela)
rel <- data.frame(
  parametro = c("b0", "b2"),
  intercepto = unname(c(coef(reg_b0)[1], coef(reg_b2)[1])),
  inclinacao = unname(c(coef(reg_b0)[2], coef(reg_b2)[2])),
  r2 = c(summary(reg_b0)$r.squared, summary(reg_b2)$r.squared),
  p_valor_inclinacao = c(summary(reg_b0)$coefficients[2, 4], summary(reg_b2)$coefficients[2, 4]),
  row.names = NULL)
cat("\n=== Parameter ~ S relations ===\n"); print(rel)
write.csv(rel, file.path(OUT_DIR, "relacoes_parametros_S.csv"), row.names = FALSE)

# -----------------------------------------------------------------------
# STEP 4 -- the single equation Hd = f(age, S)
# -----------------------------------------------------------------------
c0 <- unname(coef(reg_b0)[1]); c1 <- unname(coef(reg_b0)[2])   # b0 = c0 + c1*S
e0 <- unname(coef(reg_b2)[1]); e1 <- unname(coef(reg_b2)[2])   # b2 = e0 + e1*S
hd_pp <- function(idade, S) (c0 + c1 * S) * (1 - exp(-B1_FIXO * idade))^(e0 + e1 * S)
cat(sprintf("\nFinal equation: Hd = (%.4f + %.4f*S) * (1 - exp(-%.5f*t))^(%.4f + %.4f*S)\n",
            c0, c1, B1_FIXO, e0, e1))

d$S <- por_parcela$S_estimado[match(d$parcela, por_parcela$parcela)]
pred <- hd_pp(d$idade, d$S)
met <- metricas_ajuste(d$Hd, pred, 4, "predicao_de_parametros")
cat("\nFit of the final equation on the 43 observations (plot S taken from step 2):\n"); print(met)
write.csv(met, file.path(OUT_DIR, "ajuste_predicao_parametros.csv"), row.names = FALSE)
grafico_residuos(d$Hd, pred, d$idade,
                 file.path(OUT_DIR, "residuos_predicao_parametros.png"), "Parameter prediction")

# -----------------------------------------------------------------------
# STEP 5 -- curves and plot classification
# -----------------------------------------------------------------------
curvas <- montar_curvas(hd_pp)
write.csv(curvas, file.path(OUT_DIR, "curvas_predicao_parametros.csv"), row.names = FALSE)
grafico_feixe(curvas, d, file.path(OUT_DIR, "feixe_predicao_parametros.png"),
              paste0("Parameter prediction - site index curves (Ii = ", IDADE_INDICE, " yr)"))
cat("\nRealised Hd at the index age for each class (should be close to S):\n")
print(data.frame(S = CLASSES_S, Hd_na_idade_indice = hd_pp(IDADE_INDICE, CLASSES_S)))

S_eq <- sapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  estimar_S_parcela(g$idade, g$Hd, hd_pp)
})
S_proj <- por_parcela$S_estimado[match(parcelas, por_parcela$parcela)]
classif <- data.frame(parcela = parcelas,
                      S_projecao_curva_parcela = S_proj, classe_S = classe_mais_proxima(S_proj),
                      S_equacao_final = S_eq, classe_S_equacao_final = classe_mais_proxima(S_eq))
cat("\n=== Plot classification (parameter prediction) ===\n"); print(classif)
write.csv(classif, file.path(OUT_DIR, "indice_local_parcelas_predicao_parametros.csv"), row.names = FALSE)

cat("\nPART 5 FINISHED.\n")
