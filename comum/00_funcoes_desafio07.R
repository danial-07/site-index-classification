# =============================================================================
# 00_funcoes_desafio07.R
#
# DATA CHALLENGE 07 -- SHARED SETTINGS AND FUNCTIONS
# Site productive capacity classification (5 methods).
#
# Sourced (never run directly) by every method script. It holds the pieces
# that MUST be identical across the five methods so they can be compared:
#   - the index age (Ii) and the site-index classes (S)
#   - the base growth model (Chapman-Richards)
#   - the data reader, fit metrics, plotting helpers
# =============================================================================

# -----------------------------------------------------------------------
# SETTINGS SHARED BY ALL FIVE METHODS (documented in the README)
# -----------------------------------------------------------------------
IDADE_INDICE <- 4                      # index age Ii (years)
CLASSES_S    <- c(15, 17, 19, 21)      # site-index classes: Hd (m) at age Ii
IDADES_CURVA <- seq(2.5, 6.5, by = 0.25)  # ages used to draw the curves

# -----------------------------------------------------------------------
# BASE GROWTH MODEL: Chapman-Richards   Hd = b0 * (1 - exp(-b1*t))^b2
#   b0 = asymptote, b1 = growth rate, b2 = shape
# -----------------------------------------------------------------------
cr <- function(t, b0, b1, b2) b0 * (1 - exp(-b1 * t))^b2

# -----------------------------------------------------------------------
# DATA READER
# -----------------------------------------------------------------------
# Reads Sitio.xls and standardises the column names BY POSITION (the file
# has accents and a trailing space in "Volume ", which are fragile to type).
# A plot ID is built as "stand-plot" because plot numbers repeat across stands.
ler_sitio <- function(caminho) {
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  d <- as.data.frame(readxl::read_excel(caminho))
  names(d)[1:7] <- c("projeto", "talhao", "parcela_num", "idade",
                     "area_basal", "Hd", "volume")
  d$parcela <- paste(d$talhao, d$parcela_num, sep = "-")
  d <- d[order(d$parcela, d$idade), ]
  rownames(d) <- NULL
  d
}

# -----------------------------------------------------------------------
# FIT METRICS (always computed on the ORIGINAL Hd scale, in metres)
# -----------------------------------------------------------------------
# k = number of estimated parameters.
metricas_ajuste <- function(obs, pred, k, nome) {
  n    <- length(obs)
  sse  <- sum((obs - pred)^2)
  sst  <- sum((obs - mean(obs))^2)
  r2   <- 1 - sse / sst
  data.frame(
    modelo = nome, n = n, k = k, sse = sse,
    r2 = r2,
    r2_ajustado = 1 - (1 - r2) * (n - 1) / (n - k - 1),
    syx = sqrt(sse / (n - k)),
    syx_pct = sqrt(sse / (n - k)) / mean(obs) * 100,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------
# NONLINEAR FIT WITH SEVERAL STARTING POINTS
# -----------------------------------------------------------------------
# nls() can fail or land in a poor spot when started badly. This helper tries
# every start in `starts` (a list of named lists) with the bounded "port"
# algorithm and keeps the fit with the smallest residual sum of squares.
ajustar_nls <- function(formula, data, starts, lower, upper) {
  melhor <- NULL
  for (st in starts) {
    fit <- tryCatch(
      nls(formula, data = data, start = st, algorithm = "port",
          lower = lower, upper = upper,
          control = nls.control(maxiter = 500, warnOnly = TRUE)),
      error = function(e) NULL
    )
    if (!is.null(fit) && (is.null(melhor) || deviance(fit) < deviance(melhor))) {
      melhor <- fit
    }
  }
  if (is.null(melhor)) stop("nls failed from every starting point.")
  melhor
}

# -----------------------------------------------------------------------
# SITE INDEX CLASS HELPERS
# -----------------------------------------------------------------------
classe_mais_proxima <- function(S, classes = CLASSES_S) {
  sapply(S, function(s) classes[which.min(abs(classes - s))])
}

# Builds the long table "age x S class -> estimated Hd" from any function
# f(idade, S) (this is the format of every curvas_*.csv output).
montar_curvas <- function(f, classes = CLASSES_S, idades = IDADES_CURVA) {
  g <- expand.grid(idade = idades, S = classes)
  g$Hd_estimada_m <- f(g$idade, g$S)
  g[order(g$S, g$idade), c("idade", "S", "Hd_estimada_m")]
}

# Finds, for ONE plot, the continuous S whose curve f(idade, S) best fits that
# plot's own observations (least squares). Used wherever a plot must be
# classified from a model in which S is an explicit input.
estimar_S_parcela <- function(idade, hd, f, intervalo = c(8, 35)) {
  optimize(function(S) sum((hd - f(idade, S))^2),
           interval = intervalo, tol = 1e-8)$minimum
}

# -----------------------------------------------------------------------
# PLOT HELPERS
# -----------------------------------------------------------------------
# Residual diagnostics (required: "analise grafica de residuos"):
# residuals vs fitted values, and residuals vs age.
grafico_residuos <- function(obs, pred, idade, arquivo, titulo) {
  res <- obs - pred
  png(arquivo, width = 1500, height = 700, res = 140)
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
  plot(pred, res, pch = 16, col = "steelblue", xlab = "Fitted Hd (m)",
       ylab = "Residual (m)", main = paste(titulo, "- residuals vs fitted"))
  abline(h = 0, col = "red", lty = 2)
  plot(idade, res, pch = 16, col = "darkorange", xlab = "Age (years)",
       ylab = "Residual (m)", main = paste(titulo, "- residuals vs age"))
  abline(h = 0, col = "red", lty = 2)
  dev.off()
}

# Family of site-index curves over the observed points of the 9 plots
# (required deliverable (iv) of every script).
grafico_feixe <- function(curvas, d, arquivo, titulo) {
  parcelas <- sort(unique(d$parcela))
  cores <- rainbow(length(parcelas), s = 0.8, v = 0.8)
  png(arquivo, width = 1300, height = 900, res = 140)
  par(mar = c(4.5, 4.5, 3, 1))
  plot(d$idade, d$Hd, type = "n", xlim = c(2.3, 7), ylim = c(10, 28),
       xlab = "Age (years)", ylab = "Dominant height Hd (m)", main = titulo)
  for (i in seq_along(parcelas)) {
    g <- d[d$parcela == parcelas[i], ]
    points(g$idade, g$Hd, pch = 16, col = cores[i], cex = 1.1)
  }
  for (s in sort(unique(curvas$S))) {
    g <- curvas[curvas$S == s, ]
    lines(g$idade, g$Hd_estimada_m, lwd = 2)
    text(max(g$idade) + 0.05, tail(g$Hd_estimada_m, 1), paste0("S=", s),
         pos = 4, cex = 0.8)
  }
  abline(v = IDADE_INDICE, lty = 3, col = "grey50")
  legend("topleft", legend = paste("plot", parcelas), col = cores, pch = 16,
         ncol = 3, cex = 0.7, bty = "n")
  dev.off()
}

# Gets a plot ready: classe_S column from the continuous S estimate.
tabela_parcelas <- function(parcela, S_estimado) {
  data.frame(parcela = parcela, S_estimado = S_estimado,
             classe_S = classe_mais_proxima(S_estimado),
             stringsAsFactors = FALSE)
}
