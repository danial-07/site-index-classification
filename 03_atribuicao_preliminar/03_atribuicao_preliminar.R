# =============================================================================
# 03_atribuicao_preliminar.R
#
# DATA CHALLENGE 07 -- PART 3: PRELIMINARY SITE INDEX ASSIGNMENT
#
# Idea: (1) give each plot a preliminary site index S (here: from the guide
# curve of script 01); (2) treat S as a KNOWN covariate and fit
# Hd = f(age, S) to all observations; (3) generate the curves from that
# model; (4) re-classify the plots with the final model and compare.
#
# Input: Sitio.xls + indice_local_parcelas_curva_guia.csv (output of script 01)
#
# HOW TO USE: edit the 4 paths below, run  Rscript 03_atribuicao_preliminar.R
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                                    # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"                        # <-- change
INDICE_GUIA_FILE <- "C:/path/to/your/01_curva_guia/outputs/indice_local_parcelas_curva_guia.csv"  # <-- from script 01
OUT_DIR      <- "C:/path/to/your/03_atribuicao_preliminar/outputs"             # <-- change
# =============================================================================

for (f in c(SITIO_FILE, FUNCOES_FILE, INDICE_GUIA_FILE))
  if (!file.exists(f)) stop("File not found: ", f,
                            "\n(The guide-curve file is created by running script 01 first.)")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
source(FUNCOES_FILE)

d <- ler_sitio(SITIO_FILE)
parcelas <- sort(unique(d$parcela))

# -----------------------------------------------------------------------
# STEP 1 -- preliminary S per plot (from script 01), attached to every row
# -----------------------------------------------------------------------
prelim <- read.csv(INDICE_GUIA_FILE, stringsAsFactors = FALSE)
prelim$parcela <- as.character(prelim$parcela)
d$S <- prelim$S_media_medicoes[match(d$parcela, prelim$parcela)]
if (any(is.na(d$S))) stop("Some plots have no preliminary S - check script 01's output.")
cat("Preliminary S per plot:\n"); print(prelim[, c("parcela", "S_media_medicoes", "classe_S")])

# -----------------------------------------------------------------------
# STEP 2 -- three ways of putting S into the Chapman-Richards model
# -----------------------------------------------------------------------
#  M1: asymptote linear in S          b0 = a0 + a1*S
#  M2: asymptote as a power of S      b0 = a0 * S^a1
#  M3: asymptote AND shape linear in S   b0 = a0 + a1*S ; b2 = c0 + c1*S  (polymorphic)
fit_M1 <- ajustar_nls(
  Hd ~ (a0 + a1 * S) * (1 - exp(-b1 * idade))^b2, data = d,
  starts = list(list(a0 = 10, a1 = 1, b1 = 0.15, b2 = 0.85),
                list(a0 = 5, a1 = 1.5, b1 = 0.30, b2 = 1.00),
                list(a0 = 20, a1 = 0.5, b1 = 0.15, b2 = 0.90)),
  lower = c(-100, -10, 0.001, 0.1), upper = c(500, 50, 5, 10))
fit_M2 <- ajustar_nls(
  Hd ~ a0 * S^a1 * (1 - exp(-b1 * idade))^b2, data = d,
  starts = list(list(a0 = 3, a1 = 0.8, b1 = 0.15, b2 = 0.85),
                list(a0 = 2, a1 = 1.0, b1 = 0.30, b2 = 1.00),
                list(a0 = 5, a1 = 0.5, b1 = 0.15, b2 = 0.90)),
  lower = c(0.001, 0.01, 0.001, 0.1), upper = c(500, 5, 5, 10))
fit_M3 <- ajustar_nls(
  Hd ~ (a0 + a1 * S) * (1 - exp(-b1 * idade))^(c0 + c1 * S), data = d,
  starts = list(list(a0 = 10, a1 = 1, b1 = 0.15, c0 = 0.85, c1 = 0),
                list(a0 = 5, a1 = 1.5, b1 = 0.30, c0 = 1.00, c1 = 0),
                list(a0 = 20, a1 = 0.5, b1 = 0.15, c0 = 0.90, c1 = -0.01)),
  lower = c(-100, -10, 0.001, -5, -1), upper = c(500, 50, 5, 10, 1))

fits <- list(M1_b0_linear_em_S = fit_M1, M2_b0_potencia_de_S = fit_M2,
             M3_b0_e_b2_lineares_em_S = fit_M3)
cand <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  m <- metricas_ajuste(d$Hd, predict(f), length(coef(f)), nm)
  m$aic <- AIC(f)
  m$parametros <- paste(names(coef(f)), round(coef(f), 5), collapse = "; ")
  m
}))
cat("\n=== Candidate models (metrics on Hd) ===\n"); print(cand[, c("modelo", "r2", "r2_ajustado", "syx", "syx_pct", "aic")])

# -----------------------------------------------------------------------
# STEP 3 -- select by AIC (lower = better trade-off between fit and size)
# -----------------------------------------------------------------------
sel <- cand$modelo[which.min(cand$aic)]
fit_sel <- fits[[sel]]
cat("\nSelected model (lowest AIC):", sel, "\n"); print(summary(fit_sel)$coefficients)
write.csv(cand, file.path(OUT_DIR, "ajuste_atribuicao_preliminar.csv"), row.names = FALSE)
grafico_residuos(d$Hd, predict(fit_sel), d$idade,
                 file.path(OUT_DIR, "residuos_atribuicao_preliminar.png"),
                 "Preliminary assignment")

# model as a function of (idade, S) for any S (continuous)
cf <- as.list(coef(fit_sel))
hd_prelim <- function(idade, S) {
  if (sel == "M1_b0_linear_em_S") {
    (cf$a0 + cf$a1 * S) * (1 - exp(-cf$b1 * idade))^cf$b2
  } else if (sel == "M2_b0_potencia_de_S") {
    cf$a0 * S^cf$a1 * (1 - exp(-cf$b1 * idade))^cf$b2
  } else {
    (cf$a0 + cf$a1 * S) * (1 - exp(-cf$b1 * idade))^(cf$c0 + cf$c1 * S)
  }
}

# -----------------------------------------------------------------------
# STEP 4 -- family of curves
# -----------------------------------------------------------------------
curvas <- montar_curvas(hd_prelim)
write.csv(curvas, file.path(OUT_DIR, "curvas_atribuicao_preliminar.csv"), row.names = FALSE)
grafico_feixe(curvas, d, file.path(OUT_DIR, "feixe_atribuicao_preliminar.png"),
              paste0("Preliminary assignment - site index curves (Ii = ", IDADE_INDICE, " yr)"))
cat("\nRealised Hd at the index age for each class (should be close to S):\n")
print(data.frame(S = CLASSES_S, Hd_na_idade_indice = hd_prelim(IDADE_INDICE, CLASSES_S)))

# -----------------------------------------------------------------------
# STEP 5 -- re-classify the plots with the final model and compare
# -----------------------------------------------------------------------
# For each plot: the S whose curve fits that plot's own points best.
S_novo <- sapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  estimar_S_parcela(g$idade, g$Hd, hd_prelim)
})
S_pre <- prelim$S_media_medicoes[match(parcelas, prelim$parcela)]
comp <- data.frame(parcela = parcelas,
                   S_preliminar = S_pre, classe_preliminar = classe_mais_proxima(S_pre),
                   S_reclassificado = S_novo, classe_reclassificada = classe_mais_proxima(S_novo))
comp$mudou_classe <- comp$classe_preliminar != comp$classe_reclassificada
comp$diferenca_S <- comp$S_reclassificado - comp$S_preliminar
cat("\n=== Preliminary vs re-classified ===\n"); print(comp)
cat("\nPlots that changed class:", sum(comp$mudou_classe), "of", nrow(comp), "\n")
write.csv(comp, file.path(OUT_DIR, "indice_local_parcelas_atribuicao_preliminar.csv"), row.names = FALSE)

cat("\nPART 3 FINISHED.\n")
