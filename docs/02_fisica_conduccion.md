# 02 — Física de la conducción de calor en hormigón masivo

> Objetivo: entender **qué ecuación resolvemos y por qué**, antes de ver cómo se discretiza
> ([03_solver_1d.md](03_solver_1d.md)).

## 1. El problema en una frase

El hormigón **genera calor** al hidratarse el cemento. Ese calor viaja por **conducción**
hacia las superficies y desde ahí sale al **aire** (convección) y al **suelo**. Si el calor se
genera más rápido de lo que sale, el núcleo se calienta (riesgo **T_max > 70 °C**) mientras la
superficie se enfría con el ambiente (riesgo **ΔT núcleo–superficie > 20 °C**, que produce
tensiones y fisuras).

## 2. Propiedades térmicas

| Símbolo | Nombre | Unidades | Qué significa físicamente | Hormigón (provisional) |
|---|---|---|---|---|
| $k$ | conductividad térmica | W/(m·K) | Cuánto calor atraviesa 1 m² por cada 1 K/m de gradiente. Mide **qué tan fácil conduce**. | 2,0 |
| $\rho$ | densidad | kg/m³ | Masa por unidad de volumen. | 2400 |
| $c_p$ | calor específico | J/(kg·K) | Energía para subir **1 K a 1 kg**. Mide **cuánto calor almacena**. | 900 |
| $\rho c_p$ | capacidad térmica volumétrica | J/(m³·K) | Energía para subir 1 K a 1 m³. | 2,16·10⁶ |
| $\alpha = \dfrac{k}{\rho c_p}$ | difusividad térmica | m²/s | **Qué tan rápido se propaga** un cambio de temperatura: conduce / almacena. | 9,3·10⁻⁷ |

Para comparar: acero $k \approx 50$, madera contrachapada $k \approx 0{,}13$,
poliestireno expandido $k \approx 0{,}035$ W/(m·K) (valores típicos de referencia).
Un material con $\alpha$ grande "se entera" rápido de lo que pasa en sus bordes.

### 2.1 Longitud de difusión: por qué el hormigón masivo se calienta tanto

En un tiempo $t$, un cambio de temperatura penetra aproximadamente $\ell \approx \sqrt{\alpha t}$:

| tiempo | 1 hora | 1 día | 7 días | 14 días |
|---|---|---|---|---|
| $\sqrt{\alpha t}$ (hormigón) | 6 cm | 28 cm | 75 cm | **1,06 m** |

En los 14 días del análisis el calor "viaja" del orden de **1 m**. El centro de una fundación de
3,30 m de espesor casi no se entera de que existen las superficies: se comporta casi como si
estuviera **aislado (adiabático)** y acumula casi todo el calor que genera.
Al revés, el tiempo para que el calor cruce media fundación (1,65 m) es
$t \approx L^2/\alpha \approx 34$ días.

### 2.2 Número de Biot: ¿basta una sola temperatura?

$$Bi = \frac{h\,L}{k} = \frac{\text{resistencia interna (conducción)}}{\text{resistencia externa (convección)}}$$

Con $h \approx 10$ W/(m²K), $L = 1{,}65$ m y $k = 2$: $Bi \approx 8$.
Si $Bi \ll 0{,}1$ el cuerpo tendría temperatura casi uniforme; con $Bi \gg 1$ hay **grandes
gradientes internos**. Por eso necesitamos un modelo con variación espacial (1D → 3D), y no
un modelo de "una sola temperatura".

## 3. Balance de energía

Para cualquier volumen de control $V$ con superficie $S$:

$$
\underbrace{\frac{d}{dt}\int_V \rho c_p T \, dV}_{\text{acumulación}}
= \underbrace{-\oint_S \mathbf{q}\cdot\mathbf{n}\, dS}_{\text{entra} - \text{sale por la superficie}}
+ \underbrace{\int_V Q\, dV}_{\text{generación}}
$$

- **Ley de Fourier**: $\mathbf{q} = -k\,\nabla T$ [W/m²]. El signo menos dice que el calor fluye de
  caliente a frío.
- Reemplazando y aplicando el teorema de la divergencia se obtiene la **ecuación del calor**:

$$\rho c_p \frac{\partial T}{\partial t} = \nabla\cdot(k\nabla T) + Q$$

- En 1D y con $k$ constante: $\dfrac{\partial T}{\partial t} = \alpha \dfrac{\partial^2 T}{\partial x^2} + \dfrac{Q}{\rho c_p}$.

**Interpretación de $\partial^2T/\partial x^2$ (curvatura):** donde el perfil tiene forma de
"loma" (curvatura negativa) la temperatura baja; donde tiene forma de "valle", sube. La difusión
**suaviza** el perfil (ver [fig01](../resultados/figuras/fig01_dirichlet_vs_exacta.png)).

**Fuente $Q$ [W/m³]:** en nuestro problema es el calor de hidratación (semana 3). Si el
bloque estuviera perfectamente aislado, todo el calor generado se acumularía:
$\Delta T_{ad} = \dfrac{\int Q\,dt}{\rho c_p}$ (aumento adiabático). Por ejemplo, 8,6·10⁷ J/m³
generados equivalen a $\Delta T_{ad} = 40$ °C.

En nuestro código el balance se aplica **a cada celda** (ver [03_solver_1d.md](03_solver_1d.md)):
es exactamente esta ecuación, escrita para volúmenes pequeños.

## 4. Condición inicial

La ecuación es de **primer orden en el tiempo**: necesita conocer $T(x, t=0)$ en todo el dominio.

- En el proyecto: cada vertida empieza a existir cuando se coloca, con su **temperatura de
  colocación** (propiedad `TemperaturaColocacion_C` del IFC). Antes de eso, esas celdas no son
  hormigón activo (activación de vertidas, semana 3).
- El suelo parte con su propia temperatura inicial (semana 4).
- En el código: argumento `T_inicial` de `resolver()` (un número o un valor por celda).

## 5. Condiciones de borde

La ecuación es de **segundo orden en el espacio**: necesita una condición en **cada** borde.

| Tipo | Qué se impone | Ecuación (en el borde) | Ejemplo en el proyecto | Clase en el código |
|---|---|---|---|---|
| **Dirichlet** (1ª especie) | la **temperatura** | $T = T_b(t)$ | contacto con un medio de temperatura conocida; casos de verificación | `Dirichlet(T)` |
| **Neumann** (2ª especie) | el **flujo** de calor | $q_{\text{entrante}} = q(t)$ | $q = 0$: aislación perfecta o **plano de simetría**; radiación solar absorbida | `Neumann(q)` |
| **Robin** (3ª especie, mixta) | relación **flujo–temperatura** | $q_{\text{entrante}} = h\,(T_\infty - T_s)$ | superficie expuesta al aire, con o sin moldaje/aislación | `Robin(h, T_inf, R_extra)` |

Convención de signo del código: **$q$ positivo = calor que entra al hormigón.**

### 5.1 Robin con moldaje o aislación: resistencias en serie

El calor que sale del hormigón atraviesa capas en serie, igual que corriente por resistencias:

$$q = \frac{T_s - T_\infty}{R_{total}}, \qquad R_{total} = \underbrace{\frac{e_{\text{moldaje}}}{k_{\text{moldaje}}} + \frac{e_{\text{aislación}}}{k_{\text{aislación}}}}_{R_{extra}} + \frac{1}{h}\quad [\text{m}^2\text{K/W}]$$

| Capa (valores típicos de referencia) | Resistencia [m²K/W] |
|---|---|
| Película de aire, $h = 10$ W/(m²K) | 0,10 |
| Madera contrachapada 18 mm ($k \approx 0{,}13$) | ≈ 0,14 |
| Poliestireno 5 cm ($k \approx 0{,}035$) | ≈ 1,43 |

- $h$ crece con el **viento** (orden 5–25 W/(m²K)); se definirá con datos meteorológicos en la semana 4.
- **Retirar el moldaje** = cambiar $R_{extra}$ en ese instante: la superficie se enfría de golpe
  (choque térmico). Es una decisión constructiva del plan (`RetiroMoldaje_h`).
- [fig04](../resultados/figuras/fig04_fuente_nucleo_superficie.png) muestra el efecto con una
  fuente ilustrativa: 5 cm de aislación bajan el diferencial de ≈38 °C a ≈12 °C, pero casi no
  cambian $T_{max}$ (el núcleo es casi adiabático, §2.1).

### 5.2 Onda térmica diaria

Si el aire oscila con período $P = 24$ h, la oscilación se amortigua con la profundidad como
$e^{-z/\delta}$, con $\delta = \sqrt{2\alpha/\omega} = \sqrt{\alpha P/\pi} \approx 16$ cm.
El ciclo día-noche afecta solo los primeros centímetros (justo donde se mide la superficie, a
25 mm) y no al núcleo ([fig03](../resultados/figuras/fig03_robin_ambiente_diario.png)).

## 6. Resumen para el oral

1. Ecuación: $\rho c_p\,\partial_t T = \nabla\cdot(k\nabla T) + Q$, que es un **balance de energía**.
2. $k$ conduce, $\rho c_p$ almacena, $\alpha$ = rapidez de propagación ($\sqrt{\alpha t}$ ≈ 1 m en 14 días).
3. Condición inicial = temperatura de colocación; condiciones de borde = qué pasa en las superficies.
4. Dirichlet fija $T$, Neumann fija el flujo, Robin relaciona flujo con $T$ (aire, moldaje, aislación).
