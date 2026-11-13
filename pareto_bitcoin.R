# =============================================================================
#  Proyecto de Probabilidad — Segunda entrega (Corte 3): Aplicación en R
#  "Modelación de la concentración de riqueza y ballenas en el mercado de
#   Bitcoin mediante la distribución de Pareto"
#
#  Universidad del Rosario — Programa de Economía — Equipo G
#  Integrantes: Lina Mariana Álvarez Vargas, Tomas Montoya Suarez,
#               Daniel Esteban Moreno Orjuela
#  Docente: Lizeth Alexandra Estévez Joya
#
#  Artículo base:
#    Kondor, D., Pósfai, M., Csabai, I., & Vattay, G. (2014). Do the Rich Get
#    Richer? An Empirical Analysis of the Bitcoin Transaction Network.
#    PLOS ONE, 9(2), e86197. https://doi.org/10.1371/journal.pone.0086197
#
#  Uso de inteligencia artificial: este código se elaboró con apoyo de
#  Claude (Anthropic, modelo Claude Opus 5.5, 13 de noviembre de 2026,
#  https://claude.ai/code). Ver la sección "Uso de IA" del README.md.
#
#  Cómo ejecutarlo (desde la carpeta del proyecto):
#      Rscript pareto_bitcoin.R
#  Solo usa R base (paquetes stats y graphics). Si el paquete 'actuar' está
#  instalado, se usa además para validar nuestras funciones; si no, se omite.
#  Las figuras quedan en ./figuras y las tablas en ./resultados
# =============================================================================

set.seed(2026)                       # semilla: resultados reproducibles
dir.create("figuras",    showWarnings = FALSE)
dir.create("resultados", showWarnings = FALSE)

# Paleta de colores usada en todos los gráficos (la misma del informe:
# gris = alfa 3, naranja = alfa 1.5, azul = alfa 0.984)
col_a <- c("3" = "grey45", "1.5" = "darkorange2", "0.984" = "royalblue3")
alfas <- c(3, 1.5, 0.984)

# Función auxiliar para guardar cada gráfico en PNG
guardar_png <- function(nombre, expr, ancho = 1800, alto = 1200) {
  png(file.path("figuras", nombre), width = ancho, height = alto, res = 200)
  on.exit(dev.off())
  expr
}

linea <- function(titulo) cat("\n", strrep("=", 70), "\n ", titulo, "\n",
                              strrep("=", 70), "\n", sep = "")

# =============================================================================
# 1. IMPLEMENTACIÓN DE LA DISTRIBUCIÓN DE PARETO TIPO I:  X ~ Pareto(xm, alfa)
#    R base no trae la Pareto, así que programamos sus cuatro funciones con la
#    misma convención de R (d = densidad, p = acumulada, q = cuantil,
#    r = números aleatorios).
# =============================================================================

# F.D.P.:  f(x) = alfa * xm^alfa / x^(alfa+1),  x >= xm   (0 en otro caso)
dpareto <- function(x, xm, alfa) {
  ifelse(x >= xm, alfa * xm^alfa / x^(alfa + 1), 0)
}

# F.D.A.:  F(x) = 1 - (xm/x)^alfa,  x >= xm   (0 en otro caso)
ppareto <- function(x, xm, alfa, lower.tail = TRUE) {
  sobrevivencia <- ifelse(x >= xm, (xm / x)^alfa, 1)   # P(X > x)
  if (lower.tail) 1 - sobrevivencia else sobrevivencia
}

# Función cuantil (inversa de F):  F(x) = p  =>  x = xm * (1 - p)^(-1/alfa)
qpareto <- function(p, xm, alfa) xm * (1 - p)^(-1 / alfa)

# Números aleatorios por el MÉTODO DE LA TRANSFORMADA INVERSA:
# si U ~ Uniforme(0,1), entonces X = F^{-1}(U) ~ Pareto(xm, alfa)
rpareto <- function(n, xm, alfa) qpareto(runif(n), xm, alfa)

# Momento de orden k:  E[X^k] = alfa * xm^k / (alfa - k)  si k < alfa;
#                      si k >= alfa la integral diverge (momento infinito)
momento_k <- function(k, xm, alfa) ifelse(k < alfa, alfa * xm^k / (alfa - k), Inf)

# Media:     E[X]   = alfa*xm/(alfa-1)                    (finita si alfa > 1)
media_pareto <- function(xm, alfa) momento_k(1, xm, alfa)

# Varianza:  Var(X) = alfa*xm^2 / ((alfa-1)^2 (alfa-2))   (finita si alfa > 2)
var_pareto <- function(xm, alfa) {
  ifelse(alfa > 2, alfa * xm^2 / ((alfa - 1)^2 * (alfa - 2)), Inf)
}

# Mediana (siempre existe):  xm * 2^(1/alfa)
mediana_pareto <- function(xm, alfa) qpareto(0.5, xm, alfa)

# Gini teórico:  G = 1 / (2*alfa - 1), válido solo si alfa > 1
gini_teorico <- function(alfa) ifelse(alfa > 1, 1 / (2 * alfa - 1), NA)

# =============================================================================
# 2. VERIFICACIONES NUMÉRICAS DE LA TEORÍA DEL INFORME
# =============================================================================
linea("2. Verificaciones numéricas")

# 2.1 La densidad integra 1 (condición de normalización) para cada alfa
for (a in alfas) {
  area <- integrate(dpareto, lower = 1, upper = Inf, xm = 1, alfa = a)$value
  cat(sprintf("Integral de f(x) con xm = 1, alfa = %-5s : %.6f\n", a, area))
}

# 2.2 Media y varianza: fórmula cerrada vs. integración numérica (alfa = 3)
a <- 3; xm <- 1
EX  <- integrate(function(x) x   * dpareto(x, xm, a), xm, Inf)$value
EX2 <- integrate(function(x) x^2 * dpareto(x, xm, a), xm, Inf)$value
cat(sprintf("\nalfa = 3: E[X]  fórmula = %.6f | integral = %.6f\n",
            media_pareto(xm, a), EX))
cat(sprintf("alfa = 3: Var(X) fórmula = %.6f | E[X^2]-E[X]^2 = %.6f\n",
            var_pareto(xm, a), EX2 - EX^2))

# 2.3 Validación cruzada con el paquete 'actuar' (opcional).
#     OJO: actuar::dpareto es la Pareto tipo II (Lomax); la tipo I de nuestro
#     informe se llama actuar::dpareto1(x, shape = alfa, min = xm).
if (requireNamespace("actuar", quietly = TRUE)) {
  xs <- c(1, 1.5, 2, 5, 10, 100)
  dif <- max(abs(dpareto(xs, 1, 1.5) - actuar::dpareto1(xs, shape = 1.5, min = 1)))
  cat(sprintf("\nMáxima diferencia con actuar::dpareto1: %.2e\n", dif))
} else {
  cat("\n(Paquete 'actuar' no instalado: se omite la validación cruzada)\n")
}

# =============================================================================
# 3. GRÁFICOS DE LA F.D.P. Y LA F.D.A.  (xm = 1; alfa = 3, 1.5, 0.984)
#    Reproducen en R las figuras que el informe hizo en GeoGebra.
# =============================================================================
x <- seq(1, 6, length.out = 600)

guardar_png("01_densidad_pareto.png", {
  par(mar = c(5, 5, 4, 2))
  plot(NA, xlim = c(1, 6), ylim = c(0, 3),
       xlab = "x", ylab = "f(x)",
       main = "Función de densidad de Pareto tipo I (xm = 1)")
  for (a in alfas) lines(x, dpareto(x, 1, a), col = col_a[as.character(a)], lwd = 3)
  legend("topright", legend = paste("alfa =", alfas), col = col_a, lwd = 3,
         bty = "n", title = "Parámetro de forma")
  grid()
})

guardar_png("02_acumulada_pareto.png", {
  par(mar = c(5, 5, 4, 2))
  plot(NA, xlim = c(1, 6), ylim = c(0, 1),
       xlab = "x", ylab = "F(x) = P(X <= x)",
       main = "Función de distribución acumulada de Pareto tipo I (xm = 1)")
  for (a in alfas) lines(x, ppareto(x, 1, a), col = col_a[as.character(a)], lwd = 3)
  abline(h = 1, lty = 2)
  legend("bottomright", legend = paste("alfa =", alfas), col = col_a, lwd = 3,
         bty = "n", title = "Parámetro de forma")
  grid()
})

# =============================================================================
# 4. NÚMEROS ALEATORIOS: HISTOGRAMA vs. DENSIDAD TEÓRICA
#    Se simulan n = 10 000 valores con la transformada inversa.
# =============================================================================
n_sim <- 10000

guardar_png("03_histograma_vs_densidad.png", {
  par(mfrow = c(1, 3), mar = c(5, 4.5, 4, 1))
  for (a in alfas) {
    m <- rpareto(n_sim, 1, a)
    # Con cola pesada hay valores enormes: se muestra el tramo [1, 6], pero el
    # histograma se normaliza con TODA la muestra para que sea comparable.
    tope <- 6
    cortes <- seq(1, tope, by = 0.1)
    frec <- table(cut(m[m <= tope], cortes))
    dens_emp <- as.numeric(frec) / (n_sim * 0.1)
    plot(NA, xlim = c(1, tope), ylim = c(0, max(dens_emp, a) * 1.05),
         xlab = "x", ylab = "Densidad",
         main = sprintf("alfa = %s  (n = %d)", a, n_sim))
    rect(head(cortes, -1), 0, tail(cortes, -1), dens_emp,
         col = adjustcolor(col_a[as.character(a)], 0.35), border = "white")
    curve(dpareto(x, 1, a), add = TRUE, lwd = 3, col = col_a[as.character(a)], n = 400)
    legend("topright", c("Simulación", "Densidad teórica"), bty = "n",
           fill = c(adjustcolor(col_a[as.character(a)], 0.35), NA),
           border = c("white", NA), lwd = c(NA, 3),
           col = c(NA, col_a[as.character(a)]), cex = 0.9)
    mtext(sprintf("%.1f%% de la muestra supera x = %d", 100 * mean(m > tope), tope),
          side = 3, line = 0.2, cex = 0.7)
  }
})

# =============================================================================
# 5. MEDIA Y VARIANZA: TEÓRICAS vs. MUESTRALES
#    Aquí se ve la consecuencia de la cola pesada: si alfa <= 1 la media no
#    existe y si alfa <= 2 la varianza no existe, así que los promedios
#    muestrales NO se estabilizan (falla la ley de los grandes números).
# =============================================================================
linea("5. Media y varianza: teoría vs. simulación (xm = 1, n = 100 000)")

n_grande <- 1e5
muestras <- lapply(alfas, function(a) rpareto(n_grande, 1, a))
names(muestras) <- alfas

tabla_mom <- data.frame(
  alfa              = alfas,
  media_teorica     = media_pareto(1, alfas),
  media_muestral    = sapply(muestras, mean),
  varianza_teorica  = var_pareto(1, alfas),
  varianza_muestral = sapply(muestras, var),
  mediana_teorica   = mediana_pareto(1, alfas),
  mediana_muestral  = sapply(muestras, median)
)
print(format(tabla_mom, digits = 5), row.names = FALSE)
write.csv(tabla_mom, "resultados/momentos_teoria_vs_simulacion.csv", row.names = FALSE)

# Media acumulada X̄_n en función de n: converge solo cuando la media existe
guardar_png("04_media_acumulada.png", {
  par(mfrow = c(1, 3), mar = c(5, 4.5, 4, 1))
  for (a in alfas) {
    m  <- muestras[[as.character(a)]]
    ma <- cumsum(m) / seq_along(m)
    plot(ma, type = "l", log = "x", lwd = 2, col = col_a[as.character(a)],
         xlab = "Tamaño de muestra n (escala log)", ylab = "Media muestral acumulada",
         main = sprintf("alfa = %s", a))
    if (is.finite(media_pareto(1, a))) {
      abline(h = media_pareto(1, a), lty = 2, lwd = 2)
      legend("topright", c("Media muestral", sprintf("E[X] = %.3f", media_pareto(1, a))),
             lty = c(1, 2), lwd = 2, col = c(col_a[as.character(a)], "black"), bty = "n")
    } else {
      legend("topleft", c("Media muestral", "E[X] = infinito: no converge"),
             lty = c(1, NA), lwd = 2, col = col_a[as.character(a)], bty = "n")
    }
  }
})

# =============================================================================
# 6. FUNCIÓN GENERADORA DE MOMENTOS Y FUNCIÓN CARACTERÍSTICA
# =============================================================================
linea("6. F.G.M. y función característica")

# 6.1 Para t > 0 la F.G.M. NO existe: la integral de e^{tx} f(x) hasta un
#     límite superior L crece sin cota cuando L aumenta.
a <- 3; t_pos <- 0.5
cat("Integral de e^(t x) f(x) entre 1 y L, con t = 0.5 y alfa = 3:\n")
for (L in c(10, 50, 100, 200, 400)) {
  val <- integrate(function(x) exp(t_pos * x) * dpareto(x, 1, a), 1, L)$value
  cat(sprintf("   L = %4d  ->  %.4e\n", L, val))
}
cat("   => crece sin límite: M_X(t) no existe para t > 0.\n\n")

# 6.2 Para t < 0 sí existe:  M_X(t) = alfa (-xm t)^alfa Γ(-alfa, -xm t).
#     R base no trae la gamma incompleta superior con primer argumento
#     negativo, así que la calculamos integrando su definición.
gamma_inc_sup <- function(s, z) integrate(function(u) u^(s - 1) * exp(-u), z, Inf)$value
fgm_formula  <- function(t, xm, alfa) alfa * (-xm * t)^alfa * gamma_inc_sup(-alfa, -xm * t)
fgm_numerica <- function(t, xm, alfa)
  integrate(function(x) exp(t * x) * dpareto(x, xm, alfa), xm, Inf)$value
for (tt in c(-0.1, -0.5, -1, -2)) {
  cat(sprintf("t = %5.1f | fórmula con Γ = %.6f | integral directa = %.6f | Monte Carlo = %.6f\n",
              tt, fgm_formula(tt, 1, 3), fgm_numerica(tt, 1, 3),
              mean(exp(tt * muestras[["3"]]))))
}

# 6.3 Función característica φ(t) = E[e^{itX}] = E[cos tX] + i E[sin tX].
#     Siempre existe porque |e^{itX}| = 1. La comparamos con su versión
#     empírica (promedio sobre la muestra simulada).
#     El integrando oscila, así que se integra por tramos cortos entre xm y
#     B = 200·xm; integrando por partes, la cola omitida (x > B) es menor que
#     2·alfa·B^(-alfa-1)/t, despreciable para los t graficados.
phi_teorica <- function(t, xm, alfa, B = 200 * xm) {
  if (t == 0) return(complex(real = 1, imaginary = 0))
  nodos <- seq(xm, B, length.out = 401)
  tramo <- function(g) sum(mapply(function(a, b) integrate(g, a, b)$value,
                                  head(nodos, -1), tail(nodos, -1)))
  re <- tramo(function(x) cos(t * x) * dpareto(x, xm, alfa))
  im <- tramo(function(x) sin(t * x) * dpareto(x, xm, alfa))
  complex(real = re, imaginary = im)
}
t_grid <- seq(0, 10, length.out = 101)

guardar_png("05_funcion_caracteristica.png", {
  par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1))
  for (a in c(3, 1.5)) {
    m <- muestras[[as.character(a)]]
    phi_t <- sapply(t_grid, phi_teorica, xm = 1, alfa = a)
    phi_e <- sapply(t_grid, function(t) complex(real = mean(cos(t * m)),
                                                imaginary = mean(sin(t * m))))
    plot(t_grid, Re(phi_t), type = "l", lwd = 3, ylim = c(-1, 1),
         col = col_a[as.character(a)], xlab = "t", ylab = expression(phi(t)),
         main = sprintf("Función característica, alfa = %s", a))
    lines(t_grid, Im(phi_t), lwd = 3, lty = 2, col = col_a[as.character(a)])
    points(t_grid[c(TRUE, FALSE)], Re(phi_e)[c(TRUE, FALSE)], pch = 16, cex = 0.6)
    points(t_grid[c(TRUE, FALSE)], Im(phi_e)[c(TRUE, FALSE)], pch = 1,  cex = 0.6)
    abline(h = 0, col = "grey70")
    legend("topright", bty = "n", cex = 0.8,
           c("Parte real (teórica)", "Parte imaginaria (teórica)",
             "Parte real (simulación)", "Parte imaginaria (simulación)"),
           lty = c(1, 2, NA, NA), lwd = c(3, 3, NA, NA), pch = c(NA, NA, 16, 1),
           col = c(col_a[as.character(a)], col_a[as.character(a)], "black", "black"))
  }
})
cat("\nφ(t) calculada y comparada con la versión empírica (figura 05).\n")

# =============================================================================
# 7. APLICACIÓN: BALLENAS DE BITCOIN (Kondor et al., 2014)
#    Saldos de las direcciones con más de 50 BTC:  X ~ Pareto(xm = 50, alfa)
#    El artículo reporta densidad ∝ x^(-1.984) = x^-(alfa+1)  =>  alfa ≈ 0.984
# =============================================================================
linea("7. Aplicación: ballenas de Bitcoin, X ~ Pareto(xm = 50 BTC, alfa = 0.984)")

xm_btc   <- 50
alfa_btc <- 1.984 - 1          # exponente de la densidad menos 1
n_btc    <- 1e5                # direcciones "ballena" simuladas
saldos   <- rpareto(n_btc, xm_btc, alfa_btc)

# 7.1 Probabilidades e indicadores teóricos del modelo
umbrales <- c(100, 1000, 10000, 100000)
tabla_prob <- data.frame(
  saldo_BTC = umbrales,
  P_mayor_teorica  = ppareto(umbrales, xm_btc, alfa_btc, lower.tail = FALSE),
  P_mayor_simulada = sapply(umbrales, function(u) mean(saldos > u))
)
cat("Probabilidad de que una ballena tenga más de x BTC:\n")
print(format(tabla_prob, digits = 4), row.names = FALSE)
write.csv(tabla_prob, "resultados/probabilidades_ballenas.csv", row.names = FALSE)

cat(sprintf("\nMediana teórica  = %.2f BTC | mediana simulada = %.2f BTC\n",
            mediana_pareto(xm_btc, alfa_btc), median(saldos)))
cat(sprintf("Cuantil 99%% teórico = %.0f BTC | cuantil 99.9%% = %.0f BTC\n",
            qpareto(0.99, xm_btc, alfa_btc), qpareto(0.999, xm_btc, alfa_btc)))
cat(sprintf("E[X] = %s (alfa <= 1) y Var(X) = %s (alfa <= 2)\n",
            media_pareto(xm_btc, alfa_btc), var_pareto(xm_btc, alfa_btc)))
cat(sprintf("Fórmula de Gini 1/(2 alfa - 1) = %.4f > 1  -> no aplica (requiere alfa > 1)\n",
            1 / (2 * alfa_btc - 1)))

# 7.2 La media muestral es inestable: cinco réplicas independientes del mismo
#     modelo dan medias muy distintas, mientras la mediana casi no cambia.
replicas <- t(sapply(1:5, function(i) {
  s <- rpareto(n_btc, xm_btc, alfa_btc)
  c(replica = i, media = mean(s), mediana = median(s), maximo = max(s))
}))
cat("\nCinco réplicas de 100 000 ballenas:\n")
print(format(as.data.frame(replicas), digits = 6, big.mark = ","), row.names = FALSE)
write.csv(replicas, "resultados/replicas_media_vs_mediana.csv", row.names = FALSE)

# 7.3 Desigualdad: Gini muestral, curva de Lorenz, participación del top y
#     un índice tipo Nakamoto (mínimo de direcciones que suman el 50 %).
gini_muestral <- function(v) {
  v <- sort(v); n <- length(v)
  sum((2 * seq_len(n) - n - 1) * v) / (n * sum(v))
}
indicadores_desigualdad <- function(s) {
  desc <- sort(s, decreasing = TRUE)
  top  <- function(p) sum(desc[seq_len(ceiling(p * length(s)))]) / sum(desc)
  c(gini     = gini_muestral(s),
    top_1    = top(0.01),
    top_6.28 = top(0.0628),
    top_20   = top(0.20),
    nakamoto = which(cumsum(desc) / sum(desc) >= 0.5)[1],   # nº de direcciones
    peso_max = desc[1] / sum(desc))                        # peso de la más rica
}

# Con alfa < 1 UNA sola muestra puede estar dominada por su valor máximo, así
# que resumimos 200 réplicas de 100 000 ballenas con su mediana y el rango
# entre los percentiles 5 % y 95 %.
rep_ind <- t(replicate(200, indicadores_desigualdad(rpareto(n_btc, xm_btc, alfa_btc))))
res_ind <- data.frame(
  indicador = c("Gini", "Participación del 1% más rico",
                "Participación del 6.28% más rico", "Participación del 20% más rico",
                "Índice tipo Nakamoto (direcciones para el 50%)",
                "Peso de la dirección más rica"),
  mediana = apply(rep_ind, 2, median),
  p05     = apply(rep_ind, 2, quantile, 0.05),
  p95     = apply(rep_ind, 2, quantile, 0.95)
)
cat("\nIndicadores de desigualdad (200 réplicas de 100 000 ballenas):\n")
print(format(res_ind, digits = 4), row.names = FALSE)
write.csv(res_ind, "resultados/desigualdad_ballenas.csv", row.names = FALSE)

G_sim    <- res_ind$mediana[1]
top_628  <- res_ind$mediana[3]
nakamoto <- res_ind$mediana[5]
cat(sprintf("\nGini mediano = %.4f (artículo, todas las direcciones: 0.985)\n", G_sim))
cat(sprintf("El 6.28%% más rico de las ballenas tiene, en mediana, el %.2f%% de sus BTC\n",
            100 * top_628))
cat(sprintf("En mediana bastan %.0f de %d direcciones (%.3f%%) para reunir el 50%%\n",
            nakamoto, n_btc, 100 * nakamoto / n_btc))

# Con alfa < 1 el Gini muestral sube hacia 1 a medida que crece n: no tiene un
# valor poblacional finito distinto de 1. Se ilustra con 30 réplicas por n.
tamanos <- c(1e2, 1e3, 1e4, 1e5)
gini_n <- sapply(tamanos, function(n) {
  median(replicate(30, gini_muestral(rpareto(n, xm_btc, alfa_btc))))
})
cat("\nGini muestral mediano (30 réplicas) según el tamaño de muestra:\n")
print(data.frame(n = format(tamanos, scientific = FALSE, big.mark = ","),
                 gini = round(gini_n, 4)), row.names = FALSE)

# Curvas de Lorenz de 8 réplicas: todas muy pegadas al eje horizontal (mucha
# desigualdad) y distintas entre sí (la cola pesada hace inestable la muestra).
guardar_png("06_lorenz_ballenas.png", {
  par(mar = c(5, 5, 4, 2))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Proporción acumulada de direcciones (de menor a mayor saldo)",
       ylab = "Proporción acumulada de BTC",
       main = "Curvas de Lorenz: 8 réplicas de 100 000 ballenas (alfa = 0.984)")
  abline(0, 1, lty = 2)
  ginis <- numeric(8)
  for (r in 1:8) {
    asc <- sort(rpareto(n_btc, xm_btc, alfa_btc))
    ginis[r] <- gini_muestral(asc)
    pob <- c(0, seq_len(n_btc) / n_btc)
    riq <- c(0, cumsum(asc) / sum(asc))
    idx <- unique(round(seq(1, length(pob), length.out = 3000)))
    lines(pob[idx], riq[idx], lwd = 2, col = adjustcolor(col_a["0.984"], 0.6))
  }
  legend("topleft", bty = "n",
         c("Igualdad perfecta",
           sprintf("Réplicas (Gini entre %.3f y %.3f)", min(ginis), max(ginis))),
         lty = c(2, 1), lwd = c(1, 2), col = c("black", col_a["0.984"]))
  grid()
})

# 7.4 Firma de ley de potencias en escala log-log
#     - Sobrevivencia P(X > x) = (xm/x)^alfa  -> recta de pendiente -alfa
#     - Densidad con intervalos logarítmicos -> recta de pendiente -(alfa+1),
#       que es el exponente 1.984 que reporta el artículo (Fig. 9).
guardar_png("07_loglog_ballenas.png", ancho = 2200, expr = {
  par(mfrow = c(1, 2), mar = c(5, 4.5, 4, 1))
  # (a) sobrevivencia empírica
  xs <- sort(saldos)
  surv <- 1 - (seq_along(xs) - 1) / length(xs)
  # puntos espaciados en escala log contando desde el más rico, para que la
  # cola (pocas direcciones con saldos enormes) quede bien representada
  idx <- unique(length(xs) + 1 - round(10^seq(0, log10(length(xs)), length.out = 400)))
  plot(xs[idx], surv[idx], log = "xy", pch = 16, cex = 0.6, col = col_a["0.984"],
       xlab = "Saldo x (BTC, escala log)", ylab = "P(X > x) (escala log)",
       main = "Sobrevivencia P(X > x) = 1 - F(x)")
  curve(ppareto(x, xm_btc, alfa_btc, lower.tail = FALSE), add = TRUE, lwd = 2,
        from = xm_btc, to = max(xs))
  legend("bottomleft", bty = "n", cex = 0.85,
         c("Simulación", sprintf("Teórica: pendiente = -%.3f", alfa_btc)),
         pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), col = c(col_a["0.984"], "black"))
  # (b) densidad con intervalos logarítmicos (log-binning)
  cortes <- 10^seq(log10(xm_btc), log10(xm_btc) +
                     0.2 * ceiling(log10(max(saldos) / xm_btc) / 0.2 + 1e-9), by = 0.2)
  cuenta <- hist(saldos, breaks = cortes, plot = FALSE)$counts
  dens   <- cuenta / (n_btc * diff(cortes))
  centro <- sqrt(head(cortes, -1) * tail(cortes, -1))
  ok <- dens > 0
  # La recta se ajusta solo con intervalos que tienen al menos 10 datos: los
  # del extremo de la cola tienen 1 o 2 observaciones y sesgarían la pendiente.
  usar <- cuenta >= 10
  ajuste <- lm(log10(dens[usar]) ~ log10(centro[usar]))
  plot(centro[ok], dens[ok], log = "xy", pch = 16, col = col_a["0.984"],
       xlab = "Saldo x (BTC, escala log)", ylab = "Densidad (escala log)",
       main = "Densidad con intervalos logarítmicos")
  curve(dpareto(x, xm_btc, alfa_btc), add = TRUE, lwd = 2,
        from = xm_btc, to = max(saldos))
  legend("bottomleft", bty = "n", cex = 0.85,
         c("Simulación", sprintf("Pendiente ajustada = %.3f", coef(ajuste)[2]),
           sprintf("Teórica = -(alfa + 1) = -%.3f", alfa_btc + 1)),
         pch = c(16, NA, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
         col = c(col_a["0.984"], NA, "black"))
})

# =============================================================================
# 8. ESTIMACIÓN DE ALFA A PARTIR DE LOS DATOS (réplica ilustrativa)
#    Fingimos no conocer alfa y lo recuperamos de la muestra simulada con:
#    (a) Máxima verosimilitud: alfa_hat = n / sum(log(x_i / xm)),
#        error estándar ≈ alfa_hat / sqrt(n).
#    (b) Mínima distancia de Kolmogorov-Smirnov, la misma idea que usan
#        Kondor et al. (2014): probar varios alfa y quedarse con el que
#        minimiza la máxima diferencia entre la F.D.A. empírica y la teórica.
# =============================================================================
linea("8. Estimación de alfa con los datos simulados")

alfa_mv <- n_btc / sum(log(saldos / xm_btc))
ee_mv   <- alfa_mv / sqrt(n_btc)
cat(sprintf("Máxima verosimilitud: alfa = %.4f  (IC 95%%: %.4f – %.4f)\n",
            alfa_mv, alfa_mv - 1.96 * ee_mv, alfa_mv + 1.96 * ee_mv))

distancia_ks <- function(a, datos, xm) {
  xs <- sort(datos); n <- length(xs)
  Ft <- ppareto(xs, xm, a)
  max(pmax(seq_len(n) / n - Ft, Ft - (seq_len(n) - 1) / n))
}
grilla <- seq(0.80, 1.20, by = 0.002)
ks <- sapply(grilla, distancia_ks, datos = saldos, xm = xm_btc)
alfa_ks <- grilla[which.min(ks)]
cat(sprintf("Mínima distancia KS:  alfa = %.3f  (D = %.4f)\n", alfa_ks, min(ks)))
cat(sprintf("Valor verdadero usado en la simulación: alfa = %.3f\n", alfa_btc))
cat("Prueba KS formal con el alfa verdadero:\n")
# (runif tiene resolución finita y con 100 000 datos puede repetir algún valor;
#  ese empate aislado no afecta la prueba, por eso se silencia la advertencia)
print(suppressWarnings(ks.test(saldos, ppareto, xm = xm_btc, alfa = alfa_btc)))

guardar_png("08_estimacion_ks.png", {
  par(mar = c(5, 5, 4, 2))
  plot(grilla, ks, type = "l", lwd = 3, col = col_a["0.984"],
       xlab = expression(paste("Valor de prueba de ", alpha)),
       ylab = "Distancia de Kolmogorov-Smirnov D",
       main = "Estimación de alfa por mínima distancia KS")
  abline(v = alfa_ks, lty = 2, lwd = 2)
  abline(v = 1, lty = 3, col = "grey40")
  legend("top", bty = "n",
         c(sprintf("Mínimo: alfa = %.3f", alfa_ks), "alfa = 1 (frontera de media finita)"),
         lty = c(2, 3), lwd = c(2, 1), col = c("black", "grey40"))
  grid()
})

# =============================================================================
# 9. RESUMEN FINAL
# =============================================================================
linea("9. Resumen de la aplicación")
resumen <- data.frame(
  indicador = c("xm (BTC)", "alfa del artículo", "alfa por máxima verosimilitud",
                "alfa por mínima KS", "Mediana teórica (BTC)", "Media teórica",
                "Varianza teórica", "Gini muestral mediano (n = 100 000)",
                "Gini del artículo (todas las direcciones)",
                "Participación mediana del 6.28% más rico",
                "Direcciones que concentran el 50% (mediana)"),
  valor = c(xm_btc, alfa_btc, round(alfa_mv, 4), alfa_ks,
            round(mediana_pareto(xm_btc, alfa_btc), 2), "Infinita (alfa <= 1)",
            "Infinita (alfa <= 2)", round(G_sim, 4), 0.985,
            sprintf("%.2f%%", 100 * top_628),
            sprintf("%.0f (%.3f%%)", nakamoto, 100 * nakamoto / n_btc))
)
print(resumen, row.names = FALSE)
write.csv(resumen, "resultados/resumen_aplicacion.csv", row.names = FALSE)
cat("\nFiguras guardadas en ./figuras y tablas en ./resultados\n")
