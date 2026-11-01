# =============================================================================
# 04_hammer.R
#
# DATA CHALLENGE 07 -- PART 4: HAMMER'S METHOD
#
# Idea: do NOT follow plots through time. Instead (1) cut the data into AGE
# CLASSES; (2) inside each age class rank the plots by Hd and split them into
# PRODUCTIVITY CLASSES (here: 3 equal-sized rank groups, 1 = lowest, 3 =
# highest); (3) take the mean age and mean Hd of each group as a
# representative point; (4) fit one growth model to those points with the
# productivity class as a covariate; (5) turn the productivity class into the
# site index S by evaluating the curve at the index age.
#
# HOW TO USE: edit the 3 paths below, run  Rscript 04_hammer.R
# =============================================================================

# ============================ CONFIG (edit these) ============================
SITIO_FILE   <- "C:/path/to/your/Sitio.xls"                          # <-- change
FUNCOES_FILE <- "C:/path/to/your/00_funcoes_desafio07.R"              # <-- change
OUT_DIR      <- "C:/path/to/your/04_hammer/outputs"                   # <-- change
N_CLASSES_PROD <- 3        # productivity classes inside each age class
# Age classes [lower, upper). The last class includes its upper limit.
# They follow the five measurement rounds seen in the data (~2.6, 3.6, 4.6,
# 5.7 and 6.3 years), so each plot appears at most once per class.
CLASSES_IDADE <- rbind(c(2.5, 3.0), c(3.5, 4.0), c(4.5, 5.0), c(5.5, 6.0), c(6.0, 6.5))
# =============================================================================

if (!file.exists(SITIO_FILE))   stop("File not found: ", SITIO_FILE)
if (!file.exists(FUNCOES_FILE)) stop("File not found: ", FUNCOES_FILE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
source(FUNCOES_FILE)

d <- ler_sitio(SITIO_FILE)
parcelas <- sort(unique(d$parcela))

# -----------------------------------------------------------------------
# STEP 1 -- age class of every observation
# -----------------------------------------------------------------------
classe_idade <- function(a) {
  nc <- nrow(CLASSES_IDADE)
  for (i in seq_len(nc)) {
    dentro <- a >= CLASSES_IDADE[i, 1] && (a < CLASSES_IDADE[i, 2] ||
              (i == nc && a == CLASSES_IDADE[i, 2]))
    if (dentro) return(i)
  }
  NA_integer_
}
d$cls_idade <- sapply(d$idade, classe_idade)
if (any(is.na(d$cls_idade))) stop("Some ages fall outside CLASSES_IDADE.")
cat("Observations per age class:\n"); print(table(d$cls_idade))

# -----------------------------------------------------------------------
# STEP 2 -- productivity classes inside each age class, representative points
# -----------------------------------------------------------------------
# Plots are ordered by Hd (ties broken by plot name); the position in that
# order decides the class: class = ceiling(position * N / number_of_plots).
reps <- do.call(rbind, lapply(sort(unique(d$cls_idade)), function(ci) {
  g <- d[d$cls_idade == ci, ]
  g <- g[order(g$Hd, g$parcela), ]
  n <- nrow(g)
  g$classe_prod <- ceiling(seq_len(n) * N_CLASSES_PROD / n)
  do.call(rbind, lapply(sort(unique(g$classe_prod)), function(z) {
    gz <- g[g$classe_prod == z, ]
    data.frame(cls_idade = ci, classe_prod = z, n_parcelas = nrow(gz),
               idade = mean(gz$idade), Hd = mean(gz$Hd))
  }))
}))
cat("\n=== Representative points (age class x productivity class) ===\n"); print(reps)
write.csv(reps, file.path(OUT_DIR, "pontos_representativos_hammer.csv"), row.names = FALSE)

# -----------------------------------------------------------------------
# STEP 3 -- one model with the productivity class z as covariate
# -----------------------------------------------------------------------
#   Hd = (a0 + a1*z) * (1 - exp(-b1*t))^(b2 + b3*z)
# (asymptote and shape both vary with z, so the curves may differ in shape).
fit_h <- ajustar_nls(
  Hd ~ (a0 + a1 * classe_prod) * (1 - exp(-b1 * idade))^(b2 + b3 * classe_prod), data = reps,
  starts = list(list(a0 = 10, a1 = 5, b1 = 0.15, b2 = 0.85, b3 = 0),
                list(a0 = 20, a1 = 3, b1 = 0.30, b2 = 1.00, b3 = 0),
                list(a0 = 30, a1 = 2, b1 = 0.10, b2 = 0.80, b3 = -0.05)),
  lower = c(-100, -50, 0.001, -5, -2), upper = c(500, 50, 5, 10, 2))
cf <- as.list(coef(fit_h))
cat("\n=== Hammer model ===\n"); print(summary(fit_h)$coefficients)
met <- metricas_ajuste(reps$Hd, predict(fit_h), 5, "hammer_CR_com_classe_produtividade")
print(met)
write.csv(met, file.path(OUT_DIR, "ajuste_hammer.csv"), row.names = FALSE)
grafico_residuos(reps$Hd, predict(fit_h), reps$idade,
                 file.path(OUT_DIR, "residuos_hammer.png"), "Hammer (representative points)")

hd_z <- function(idade, z) (cf$a0 + cf$a1 * z) * (1 - exp(-cf$b1 * idade))^(cf$b2 + cf$b3 * z)

# -----------------------------------------------------------------------
# STEP 4 -- connect the productivity class to the site index S
# -----------------------------------------------------------------------
# S = Hd of the class-z curve at the index age. To draw the curve for a given S
# (e.g. 17 m) we need the (continuous) z that gives exactly that Hd at Ii.
z_de_S <- function(S) optimize(function(z) (hd_z(IDADE_INDICE, z) - S)^2,
                               interval = c(-3, 8), tol = 1e-10)$minimum
zS <- sapply(CLASSES_S, z_de_S)
cat("\nz needed for each S class:\n"); print(data.frame(S = CLASSES_S, z = zS))

curvas <- do.call(rbind, lapply(seq_along(CLASSES_S), function(i)
  data.frame(idade = IDADES_CURVA, S = CLASSES_S[i], Hd_estimada_m = hd_z(IDADES_CURVA, zS[i]))))
write.csv(curvas, file.path(OUT_DIR, "curvas_hammer.csv"), row.names = FALSE)
grafico_feixe(curvas, d, file.path(OUT_DIR, "feixe_hammer.png"),
              paste0("Hammer - site index curves (Ii = ", IDADE_INDICE, " yr)"))

# -----------------------------------------------------------------------
# STEP 5 -- classify the 9 plots (z that best fits each plot's own points)
# -----------------------------------------------------------------------
res <- do.call(rbind, lapply(parcelas, function(pl) {
  g <- d[d$parcela == pl, ]
  z <- optimize(function(z) sum((g$Hd - hd_z(g$idade, z))^2), interval = c(-3, 8), tol = 1e-10)$minimum
  S <- hd_z(IDADE_INDICE, z)
  data.frame(parcela = pl, n_medicoes = nrow(g), indice_classe_prod_z = z,
             S_estimado = S, classe_S = classe_mais_proxima(S))
}))
cat("\n=== Plot classification (Hammer) ===\n"); print(res)
write.csv(res, file.path(OUT_DIR, "indice_local_parcelas_hammer.csv"), row.names = FALSE)

cat("\nPART 4 FINISHED.\n")
