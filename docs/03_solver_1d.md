# 03 — Solver 1D: discretización, estabilidad y verificación

> Código: [`src/mcplan/solver1d.py`](../src/mcplan/solver1d.py) ·
> Tests: [`tests/test_solver1d.py`](../tests/test_solver1d.py) ·
> Física: [02_fisica_conduccion.md](02_fisica_conduccion.md)

## 1. Idea: un balance de energía por celda

El enunciado pide una "grilla cartesiana regular con balance térmico por celda". Eso es un
esquema de **volúmenes finitos centrados en celdas**: dividimos el dominio en celdas y a cada
una le aplicamos el balance *acumulación = entra − sale + genera*.

¿Por qué celdas (y no nodos en los bordes, como el código de ejemplo del curso)?

- En 3D cada **vóxel será una celda**: el solver 1D ya tiene la misma lógica que tendrá el 3D.
- **Conserva la energía exactamente**: lo que sale de una celda entra a su vecina.
- Materiales distintos (hormigón–suelo) y moldaje/aislación se tratan con **resistencias en serie**.

## 2. Malla

```
 cara 0      cara 1      cara 2                    cara N
   |  celda 0  |  celda 1  |    ...    | celda N-1  |
  x=0   •x_0        •x_1                    •x_{N-1}  x=L
```

- $N$ celdas de ancho $\Delta x = L/N$; centro de la celda $i$: $x_i = (i + \tfrac12)\Delta x$.
- La incógnita $T_i$ es la temperatura de la celda $i$.
- Trabajamos **por m² de sección** ($A = 1$ m²), por eso las unidades quedan "por m²".
- Detalle útil: con $\Delta x = 5$ cm el centro de la primera celda queda **a 25 mm de la
  superficie**, justo la profundidad del criterio $S_{25}$ del curso.

Código: clase `Malla1D` (`Malla1D.uniforme(...)` o `Malla1D.por_capas(...)`).

## 3. Balance de la celda $i$

$$C_i \frac{dT_i}{dt} = F_{i} - F_{i+1} + Q_i\,\Delta x$$

| Símbolo | Significado | Unidades | Código |
|---|---|---|---|
| $C_i = \rho_i c_{p,i} \Delta x$ | capacidad térmica de la celda | J/(m²K) | `malla.capacidad` |
| $F_j$ | flujo por la cara $j$, **positivo hacia $+x$** | W/m² | arreglo `F` (N+1 caras) |
| $F_i - F_{i+1}$ | lo que entra por la izquierda menos lo que sale por la derecha | W/m² | `F[:-1] - F[1:]` |
| $Q_i \Delta x$ | calor generado en la celda | W/m² | `Q * malla.dx` |

## 4. Flujo entre celdas: conductancias

$$F_{i+1} = G_{i+1}\,(T_i - T_{i+1}), \qquad G = \frac{1}{R}, \qquad R = \frac{\Delta x}{2k_i} + \frac{\Delta x}{2k_{i+1}}$$

El calor va del centro de la celda $i$ a su cara (media celda) y de la cara al centro de la
celda $i+1$ (otra media celda): **dos resistencias en serie**. Código: `malla.conductancia_interior`.

- Si $k$ es igual en ambas: $G = k/\Delta x$ y se recupera la diferencia finita clásica:
  $\dfrac{dT_i}{dt} = \alpha\dfrac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2}$.
- Si los materiales son distintos, la resistencia en serie (media armónica de $k$) garantiza
  que el flujo sea continuo en la interfaz. Un promedio aritmético de $k$ sería incorrecto
  (por ejemplo, entre hormigón y aislación). Lo verifica `test_estacionario_dos_materiales_en_serie`.

## 5. Bordes

La celda del borde está a $\Delta x/2$ de la superficie. Su conductancia hacia "afuera":

| Borde | Conductancia $G_b$ | Flujo que entra |
|---|---|---|
| Dirichlet $T_b$ | $\dfrac{1}{\Delta x/2k} = \dfrac{2k}{\Delta x}$ | $G_b\,(T_b - T_0)$ |
| Robin $(h, T_\infty, R_{extra})$ | $\dfrac{1}{\Delta x/2k + R_{extra} + 1/h}$ | $G_b\,(T_\infty - T_0)$ |
| Neumann $q$ | 0 (no depende de $T$) | $q$ |

Código: `_conductancia_borde()` y `_flujo_entrante()`. Como $F$ es positivo hacia $+x$:
por la izquierda `F[0] = +q_entrante` y por la derecha `F[-1] = -q_entrante`.

## 6. Avance en el tiempo: Euler explícito

$$T_i^{n+1} = T_i^n + \frac{\Delta t}{C_i}\left(F_i^n - F_{i+1}^n + Q_i^n\Delta x\right)$$

**Explícito** = la temperatura nueva se calcula directamente con valores del paso anterior (no
hay que resolver un sistema de ecuaciones). Es simple y fácil de explicar, pero tiene un límite
de estabilidad para $\Delta t$. Precisión: orden 1 en el tiempo y orden 2 en el espacio (§8).

En el código es una sola línea dentro del bucle de `resolver()`:
`T = T + dt * (F[:-1] - F[1:] + Q * malla.dx) / C`.

## 7. Estabilidad

### 7.1 Criterio de coeficientes positivos (el que usa el código)

Reordenando la actualización de la celda $i$ (sin fuente):

$$T_i^{n+1} = \left(1 - \frac{\Delta t\,S_i}{C_i}\right) T_i^n + \sum_{\text{vecinos}} \frac{\Delta t\,G}{C_i}\,T_{\text{vecino}}^n,\qquad S_i = \sum G \text{ de la celda } i$$

Si **todos los coeficientes son ≥ 0**, la nueva temperatura es un **promedio ponderado** de
temperaturas anteriores: nunca puede salirse del rango que ya existía (*principio del máximo
discreto*) y los errores no crecen. Si $1 - \Delta t S_i/C_i < 0$, la celda responde "al revés"
de su valor anterior y aparece un serrucho que crece. Condición:

$$\Delta t \le \min_i \frac{C_i}{S_i}$$

| Celda | $S_i$ | Límite en $Fo = \alpha\Delta t/\Delta x^2$ |
|---|---|---|
| interior (1D) | $2k/\Delta x$ | $Fo \le 1/2$ |
| junto a un borde **Dirichlet** | $k/\Delta x + 2k/\Delta x$ | $Fo \le 1/3$ (más estricto: la superficie está a $\Delta x/2$) |
| junto a un borde Robin | entre los dos anteriores | entre 1/3 y 1/2 |
| interior 2D / 3D | $4k/\Delta x$ / $6k/\Delta x$ | $Fo \le 1/4$ / $Fo \le 1/6$ |

Código: `paso_tiempo_maximo()`. Si no se da `dt`, `resolver()` usa el 90 % del máximo; si se
pide un `dt` mayor, **lanza un error** (salvo `permitir_inestable=True`, solo para demostrarlo).

### 7.2 Análisis de von Neumann (la versión clásica)

Para una malla uniforme infinita, un error con forma $e^{\mathrm{i}\theta j}$ se multiplica en cada
paso por el factor de amplificación

$$g(\theta) = 1 - 4\,Fo\,\sin^2(\theta/2)$$

Se necesita $|g| \le 1$ para todo $\theta$. El peor caso es $\theta = \pi$ (el serrucho
$+,-,+,-$), con $g = 1 - 4Fo$: estable si $Fo \le 1/2$. Con $Fo = 0{,}55$, $g = -1{,}2$: el
serrucho **cambia de signo y crece un 20 % en cada paso**
([fig05](../resultados/figuras/fig05_estabilidad.png)).

### 7.3 Números para el hormigón ($\alpha \approx 9{,}3\cdot10^{-7}$ m²/s)

- 1D, $\Delta x = 5$ cm, borde Dirichlet ($Fo \le 1/3$): $\Delta t_{max} \approx 900$ s = 15 min.
- 3D, $\Delta x = 10$ cm ($Fo \le 1/6$): $\Delta t_{max} \approx 1800$ s = 30 min → 14 días ≈ 700 pasos.

**Estable no significa preciso**: la estabilidad solo garantiza que los errores no crezcan. La
precisión se mide con la verificación.

## 8. Verificación

*Verificar* = demostrar que el programa resuelve bien la ecuación (contra respuestas conocidas).
*Validar* = demostrar que la ecuación representa la realidad (contra mediciones; más adelante,
con Ju y Lei 2019).

| Test (`tests/test_solver1d.py`) | Qué demuestra |
|---|---|
| `test_solucion_exacta_seno_dirichlet` | Coincide con la solución exacta $T = T_b + A\sin(\pi x/L)e^{-\alpha\pi^2t/L^2}$ |
| `test_convergencia_espacial_orden_2` | Duplicar $N$ divide el error por 4 (orden 2) |
| `test_conservacion_energia_bordes_adiabaticos` | Sin flujos ni fuente, la energía total no cambia (error ~10⁻¹⁶) |
| `test_balance_con_flujo_neumann_y_fuente` | Energía ganada = calor entrante por el borde + calor generado |
| `test_fuente_uniforme_adiabatica_aumento_lineal` | Bloque aislado con $Q$: $T = T_0 + Qt/\rho c_p$ |
| `test_robin_tiende_a_temperatura_del_aire` | Con Robin, todo termina a $T_\infty$ |
| `test_estacionario_dirichlet_robin_con_moldaje` | Régimen permanente con resistencias en serie (aire + moldaje) |
| `test_estacionario_dos_materiales_en_serie` | Interfaz entre materiales distintos |
| `test_paso_maximo_*` | Límites $Fo \le 1/2$ (interior) y $Fo \le 1/3$ (borde Dirichlet) |
| `test_esquema_estable_respeta_principio_del_maximo` | Con $Fo \le 1/2$ no aparecen temperaturas fuera del rango inicial |
| `test_esquema_inestable_explota` | Con $Fo = 0{,}6$ la solución explota |

**Estudio de convergencia** (`scripts/04_convergencia_1d.py`,
[fig06](../resultados/figuras/fig06_convergencia.png)):

- Espacio ($Fo = 0{,}25$ fijo, $N = 10 \to 320$): orden observado **2,00**.
- Tiempo ($N = 50$ fijo, $\Delta t = \Delta t_{max} \to \Delta t_{max}/32$): orden observado **1,00**.
  Para aislar el error temporal se compara con una referencia en la **misma malla** (así el
  error espacial es idéntico y se cancela), obtenida por extrapolación de Richardson:
  si $T(\Delta t) \approx T^* + C\,\Delta t$, entonces $2T(h) - T(2h) \approx T^*$.

## 9. Cómo usarlo

```python
import numpy as np
from mcplan.solver1d import Material, Malla1D, Dirichlet, Neumann, Robin, resolver

hormigon = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)
malla = Malla1D.uniforme(L=3.30, n=66, material=hormigon)      # dx = 5 cm

res = resolver(
    malla,
    T_inicial=20.0,                                   # temperatura de colocación [°C]
    borde_izq=Robin(h=10.0, T_inf=15.0, R_extra=0.14),  # superficie con moldaje
    borde_der=Neumann(0.0),                           # base adiabática (provisional)
    t_final=14 * 86400.0,                             # 14 días [s]
    fuente=lambda t, x: 1000.0 * np.exp(-t / 86400),  # Q(t, x) [W/m³] ilustrativa
)
res.x, res.t, res.T      # posiciones, instantes guardados, temperaturas (fila = instante)
res.fourier, res.dt      # número de Fourier y paso usado
```

## 10. Limitaciones actuales (próximas semanas)

- Es 1D: la extensión 2D/3D (semana 2) reutiliza la misma idea de celdas y conductancias.
- La fuente $Q$ no depende de la temperatura; la hidratación real sí (semana 3).
- $R_{extra}$ es constante; el retiro de moldaje requiere que cambie en el tiempo (semana 4).
- No hay activación de vertidas ni suelo explícito (semanas 3–4).
