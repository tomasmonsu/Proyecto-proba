# Proyecto de Probabilidad — Entrega 2 (Corte 3): Aplicación en R

**Modelación de la concentración de riqueza y ballenas en el mercado de Bitcoin mediante la distribución de Pareto**

Universidad del Rosario · Programa de Economía · Equipo G
Integrantes: Lina Mariana Álvarez Vargas, Tomas Montoya Suarez, Daniel Esteban Moreno Orjuela
Docente: Lizeth Alexandra Estévez Joya · Exposición: 16 al 20 de noviembre de 2026

Artículo base: Kondor, D., Pósfai, M., Csabai, I., & Vattay, G. (2014). *Do the Rich Get Richer? An Empirical Analysis of the Bitcoin Transaction Network*. PLOS ONE, 9(2), e86197. https://doi.org/10.1371/journal.pone.0086197

---

## 1. Cómo ejecutar

```r
# En RStudio: abrir la carpeta del proyecto como directorio de trabajo y
source("pareto_bitcoin.R")
```
o desde la terminal: `Rscript pareto_bitcoin.R` (unos 15 segundos).

- Solo usa **R base** (no hay que instalar nada). Si el paquete `actuar` está instalado, el script lo usa además para validar nuestras funciones (`actuar::dpareto1` es la Pareto tipo I; ojo: `actuar::dpareto` es la tipo II).
- Semilla fija (`set.seed(2026)`): los resultados son reproducibles.
- Salidas: gráficos en `figuras/`, tablas en `resultados/` y la salida completa de consola en `resultados/salida_consola.txt`.

## 2. Qué hace el script (cumple la sección 6 de la guía)

| Sección | Contenido | Salida |
|---|---|---|
| 1 | Implementación propia de `dpareto`, `ppareto`, `qpareto`, `rpareto` (transformada inversa), momento de orden k, media, varianza, mediana y Gini | — |
| 2 | Verificación numérica: la densidad integra 1; media y varianza por fórmula = por integración | consola |
| 3 | F.D.P. y F.D.A. con xm = 1 y α = 3; 1.5; 0.984 (las mismas curvas del informe, ahora en R) | `01`, `02` |
| 4 | 10 000 números aleatorios: histograma vs. densidad teórica | `03` |
| 5 | Media y varianza teóricas vs. muestrales; media acumulada (ley de los grandes números) | `04`, CSV |
| 6 | F.G.M.: divergencia para t > 0, fórmula con Γ(−α, −xm·t) para t < 0; función característica teórica vs. empírica | `05` |
| 7 | **Aplicación Bitcoin**: X ~ Pareto(xm = 50 BTC, α = 0.984); probabilidades, cuantiles, inestabilidad de la media, Gini, Lorenz, top 6.28 %, índice tipo Nakamoto, gráficos log-log | `06`, `07`, CSV |
| 8 | **Réplica del método del artículo**: estimar α por máxima verosimilitud y por mínima distancia de Kolmogorov-Smirnov | `08` |
| 9 | Tabla resumen | CSV |

## 3. Resultados e interpretación

### 3.1 La distribución (xm = 1)
- La densidad integra **1.000000** para los tres valores de α, y para α = 3 la fórmula da E[X] = 1.5 y Var(X) = 0.75, iguales a la integración numérica.
- **Media y varianza (n = 100 000):**

| α | E[X] teórica | Media muestral | Var teórica | Var muestral | Mediana teórica | Mediana muestral |
|---|---|---|---|---|---|---|
| 3 | 1.5 | 1.504 | 0.75 | 0.72 | 1.260 | 1.262 |
| 1.5 | 3.0 | 2.982 | ∞ | 400.6 | 1.587 | 1.591 |
| 0.984 | ∞ | 14.60 | ∞ | 357 740 | 2.023 | 2.014 |

- **Figura 04 (idea clave para la exposición):** con α = 3 la media acumulada se estabiliza en 1.5. Con α = 1.5 se acerca a 3 pero con saltos (la varianza infinita la hace lenta). Con α = 0.984 **nunca se estabiliza**: cada ballena nueva la hace saltar. Así se ve en la práctica que la media "es infinita": no hay un valor al que converja. La **mediana**, en cambio, siempre se estima bien.

### 3.2 F.G.M. y función característica
- Para t = 0.5, ∫₁ᴸ e^{tx} f(x) dx pasa de 2.4 (L = 10) a 1.7·10⁷⁷ (L = 400): **la F.G.M. no existe para t > 0**, como decía el informe.
- Para t < 0, la fórmula α(−xm·t)^α Γ(−α, −xm·t) coincide con la integral directa y con Monte Carlo (p. ej. t = −1: 0.258187 / 0.258187 / 0.2577).
- Como la guía pide la **función característica** cuando la F.G.M. no existe: φ(t) = E[e^{itX}] siempre existe porque |e^{itX}| = 1. La figura 05 muestra que la φ(t) teórica (integración numérica) coincide con la empírica de la simulación.

### 3.3 Aplicación: ballenas de Bitcoin (xm = 50 BTC, α = 0.984)
- **Probabilidades:** P(X > 100) = 0.506; P(X > 1 000) = 0.052; P(X > 10 000) = 0.0054; P(X > 100 000) = 0.00056. La simulación reproduce estos valores.
- **Mediana = 101 BTC**, pero el cuantil 99 % es 5 389 BTC y el 99.9 % es 55 944 BTC: la cola es enorme.
- **Media infinita en la práctica:** cinco simulaciones del mismo modelo dan medias de 592, 598, 616, 834 y 1 235 BTC (la mediana siempre ≈ 101 BTC). El promedio depende de cuál fue la ballena más grande que salió.
- **Gini:** la fórmula 1/(2α − 1) = 1.033 > 1, imposible, porque requiere α > 1 (lo que dice el informe). El Gini de las muestras sí se puede calcular: mediana **0.862** con 100 000 ballenas, y **sube con el tamaño de la muestra** (0.67 → 0.77 → 0.83 → 0.86 para n = 10² … 10⁵): con α < 1 no tiene un valor poblacional estable.
- **Concentración (200 réplicas, mediana):** el 1 % más rico de las ballenas tiene el 67.7 % de sus BTC, el 6.28 % más rico tiene el 80.8 % y el 20 % más rico el 88.9 %. En mediana **bastan 94 direcciones de 100 000 (0.09 %) para reunir la mitad de todos los BTC** (índice tipo Nakamoto) y la dirección más rica tiene sola el 11.8 %.
- **Comparación con el artículo (honestidad metodológica):** el 93.72 % y el Gini de 0.985 de Kondor et al. se calculan sobre *todas* las direcciones (incluidas las de saldo pequeño), mientras que nuestro modelo describe solo la cola (> 50 BTC). Por eso nuestros números de concentración *dentro* de las ballenas son algo menores. Aun así, muestran que incluso entre las ballenas la riqueza está extremadamente concentrada.
- **Figura 07 (log-log):** la sobrevivencia P(X > x) es una recta de pendiente −0.984 y la densidad (con intervalos logarítmicos, como la Fig. 9 del artículo) una recta de pendiente ajustada **−1.963** frente a la teórica −(α + 1) = **−1.984**. Esta es la "firma" de una ley de potencias y explica por qué el artículo reporta 1.984 y nosotros usamos α = 0.984.

### 3.4 Estimación de α (réplica del método del artículo)
- **Máxima verosimilitud:** α̂ = n / Σ ln(xᵢ/xm) = **0.9844** (IC 95 %: 0.9783 – 0.9905).
- **Mínima distancia KS** (la idea de Kondor et al.: probar varios α y quedarse con el que minimiza la distancia de Kolmogorov-Smirnov): **α = 0.986** (figura 08).
- Ambos recuperan el valor verdadero 0.984, y la prueba KS no rechaza el modelo (p = 0.71). El intervalo de confianza queda entero por debajo de 1: con datos como los del artículo, la conclusión "media infinita" es estadísticamente firme.

## 4. Guion sugerido para la exposición (10 minutos con preguntas)

Cada integrante debe poder responder sobre todo el proyecto (la docente pregunta a cualquiera). Reparto sugerido:

| Tiempo | Quién | Contenido | Apoyo |
|---|---|---|---|
| 0:00–2:00 | Integrante 1 | Problema: ¿por qué no una Normal? Ballenas, Gini, Nakamoto, efecto Mateo vs. ley de Gibrat; artículo de Kondor et al. | Informe |
| 2:00–4:30 | Integrante 2 | La distribución: f(x), F(x), parámetros xm y α, momento de orden k, media (α > 1), varianza (α > 2), F.G.M. inexistente para t > 0 → función característica. Mostrar secciones 1–2 del código. | Fig. 01, 02, 05 |
| 4:30–7:00 | Integrante 3 | Código en R: transformada inversa, histograma vs. densidad, media acumulada que no converge con α = 0.984. | Fig. 03, 04 |
| 7:00–8:30 | Integrante 1 o 2 | Aplicación: α = 1.984 − 1, probabilidades, mediana vs. media, Lorenz, 94 direcciones con el 50 %, log-log, estimación por MV y KS. | Fig. 06, 07, 08 |
| 8:30–10:00 | Todos | Conclusión y preguntas. Mencionar el uso de IA en el código. | — |

**Conclusión de una frase:** con α ≈ 0.984 < 1 la riqueza de las ballenas de Bitcoin no tiene media ni varianza finitas; el promedio no describe al "usuario típico" (la mediana sí) y unas pocas direcciones concentran la mitad de los bitcoins, lo que es coherente con un mecanismo de "los ricos se hacen más ricos".

**Preguntas probables y respuesta corta:**
- *¿Por qué α = 0.984 y no 1.984?* La densidad de Pareto es ∝ x^−(α+1); el artículo reporta el exponente de la densidad, 1.984 = α + 1.
- *¿Cómo genera R números de Pareto si no viene en R base?* Transformada inversa: si U ~ U(0,1), X = xm·(1 − U)^(−1/α).
- *¿Por qué la media muestral da un número si la media teórica es infinita?* Toda muestra finita tiene promedio finito, pero no converge: crece y salta al aumentar n (figura 04 y las cinco réplicas).
- *¿Por qué no existe la F.G.M.?* e^{tx} crece más rápido que cualquier potencia x^{α+1}, así que la integral diverge para t > 0. La función característica sí existe porque |e^{itx}| = 1.
- *¿Por qué su Gini (0.86) es menor que el del artículo (0.985)?* El del artículo incluye todas las direcciones; nosotros modelamos solo la cola de más de 50 BTC.

## 5. Uso de inteligencia artificial

De acuerdo con la guía, se declara que el código en R (`pareto_bitcoin.R`) y este documento se elaboraron con apoyo de un asistente de IA generativa, a partir de la teoría y la aplicación desarrolladas por el equipo en la primera entrega. El equipo revisó, ejecutó y verificó los resultados. Para el informe/presentación, se sugiere citarlo así (APA 7):

> Anthropic. (2026). *Claude* (versión Claude Opus 5.5) [Modelo de lenguaje de gran tamaño]. Consultado el 13 de noviembre de 2026. https://claude.ai/code

Prompt utilizado (adjuntar como evidencia en el anexo junto con capturas legibles de la conversación): *"Necesito que hagas eso"*, adjuntando la guía del proyecto y el informe corregido de la entrega 1.

## 6. Estructura del repositorio

```
pareto_bitcoin.R      Script principal (comentado paso a paso)
figuras/              8 gráficos PNG para la exposición
resultados/           Tablas CSV y salida completa de consola
README.md             Este documento
```
