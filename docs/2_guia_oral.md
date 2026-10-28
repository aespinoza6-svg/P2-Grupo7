# Guía para los orales

Una sección por revisión oral. Cada pregunta tiene una respuesta modelo y qué mostrar en el
código. La teoría completa está en [1_teoria.md](1_teoria.md).

## E01 — BIM, difusión y solver 1D (miércoles 21 de octubre)

> En la revisión **se puede preguntar a cualquier integrante** y se pide **navegar por código
> real**. Cada respuesta indica qué mostrar. Las respuestas son guías: hay que poder decirlas
> con palabras propias.
>
> Lecturas previas exigidas: guías del curso G00, G01, G02 y G05.

### Checklist individual (cada integrante, antes del miércoles)

- [ ] Clonar el repo, crear el entorno e instalar (`README.md` → *Instalación*).
- [ ] Correr `pytest` y revisar que todas las pruebas pasen.
- [ ] Correr los scripts 01 a 04 y mirar cada figura en `resultados/figuras/`.
- [ ] Abrir el IFC en Bonsai, cambiar una propiedad, guardar como otra versión y compararla
      con `scripts/02_inspeccionar_ifc.py --comparar`.
- [ ] Cambiar `dt` en un ejemplo para que `Fo > 1/2` y ver el error / la explosión.
- [ ] Poder responder **todas** las preguntas de abajo sin leer.

---

### Bloque A — BIM / IFC

**A1. ¿Qué es BIM?**
Una forma de trabajar con un modelo digital de **información** del proyecto: cada objeto es un
elemento constructivo con identidad, tipo, relaciones y propiedades (no solo un dibujo). En
nuestro caso describe qué se construye y cuándo: la fundación, cada vertida, su material, inicio
planificado, temperatura de colocación, moldaje y aislación.
*Mostrar:* `fundaciones/fundacion.ifc` en Bonsai o la salida de `scripts/02_inspeccionar_ifc.py`.

**A2. ¿Qué es IFC?**
El formato **abierto** e internacional (ISO 16739, buildingSMART) para intercambiar modelos BIM.
Es un archivo de texto con entidades (`IfcFooting`, `IfcSite`…), relaciones y property sets.
Al ser abierto, cualquier programa lo lee: Bonsai, Revit o nuestro Python.
*Mostrar:* abrir el `.ifc` en un editor de texto y ubicar la entidad `IFCFOOTING` de V01.

**A3. ¿Qué rol cumple Blender?**
Es la plataforma de modelado 3D (interfaz gráfica, visualización). Por sí solo no entiende IFC
ni sabe qué es una fundación.

**A4. ¿Qué agrega Bonsai?**
Convierte Blender en una herramienta **BIM nativa IFC**: lo que se edita es directamente el IFC.
Agrega clases IFC, árbol espacial (sitio → edificio → nivel), property sets, materiales y
planificación. Por dentro usa IfcOpenShell, la misma librería que usamos en Python.

**A5. ¿Por qué el solver no depende de Blender?**
Porque el único puente es el **archivo IFC**: (1) cualquiera puede reproducir el cálculo sin
instalar Blender, incluido el cluster HPC sin pantalla; (2) el estado interno de Blender no es
un formato estable ni auditable, el IFC sí; (3) se puede testear con `pytest`; (4) cada vertida
tiene un `GlobalId` permanente y se pueden comparar versiones.
*Mostrar:* `mcplan/__init__.py` (regla: ningún módulo importa Blender) y `bim.leer_vertidas()`.

**A6. ¿Dónde están las propiedades de cada vertida?** En el property set `UANDES_HormigonMasivo` de cada
`IfcFooting` (`VertidaID`, `TipoHormigon`, `Inicio_h`, `InicioPlanificado`, `TemperaturaColocacion_C`,
`Moldaje`, `Aislacion`, `RetiroMoldaje_h`). Las del material en los psets estándar
`Pset_MaterialCommon` y `Pset_MaterialThermal`.
*Mostrar:* `bim.crear_modelo_fundacion()`, paso 5.

**A7. ¿Por qué `UANDES_HormigonMasivo` y no `Pset_Vertida`?** El prefijo `Pset_` está reservado
para los property sets estándar de buildingSMART; el set docente propio usa un nombre distinto.

**A8. ¿Qué es el GlobalId y por qué lo fijamos?** Es la identidad única y permanente de cada
objeto IFC. IfcOpenShell lo genera al azar; nosotros lo calculamos con `uuid5` a partir de una
clave (`"vertida/V01"`) para que regenerar el modelo no cambie las identidades.
*Mostrar:* `bim._guid_estable()` y `test_globalid_estable_al_regenerar`.

**A9. ¿Qué dimensiones tiene la geometría docente?** Planta octagonal con 18 m entre caras,
1,10 m de espesor en el borde, 3,30 m al centro y pedestal cilíndrico Ø7 m × 1 m; por eso la
altura total es 4,30 m. Es una simplificación docente, no una reconstrucción CAD del artículo.
*Mostrar:* `config/caso_base.yaml` y [3_supuestos.md](3_supuestos.md).

---

### Bloque B — Física

**B1. Significado físico de k, ρ, c_p, α**
- $k$ [W/(m·K)]: qué tan fácil **conduce** calor (flujo por unidad de gradiente).
- $\rho$ [kg/m³]: masa por volumen.
- $c_p$ [J/(kg·K)]: energía para subir 1 K a 1 kg; mide cuánto calor **almacena**.
- $\alpha = k/(\rho c_p)$ [m²/s]: **qué tan rápido se propaga** un cambio de temperatura.
  Para el hormigón $\alpha \approx 9{,}3\cdot10^{-7}$ m²/s: en 14 días el calor penetra solo
  $\sqrt{\alpha t} \approx 1$ m, por eso el núcleo de 3,30 m se calienta casi como si estuviera aislado.
*Mostrar:* clase `Material` en `solver1d.py`; [teoría B.2](1_teoria.md).

**B2. Balance de energía**
Acumulación = lo que entra − lo que sale + lo que se genera:
$\frac{d}{dt}\int_V \rho c_p T\,dV = -\oint_S \mathbf{q}\cdot\mathbf{n}\,dS + \int_V Q\,dV$, con
Fourier $\mathbf{q} = -k\nabla T$. Resulta $\rho c_p \partial_t T = \nabla\cdot(k\nabla T) + Q$.
En el código se aplica a cada celda: $C_i\,dT_i/dt = F_i - F_{i+1} + Q_i\Delta x$.
*Mostrar:* línea `T = T + dt * (F[:-1] - F[1:] + Q * malla.dx) / C` en `resolver()`.

**B3. Condición inicial**
Como la ecuación es de primer orden en el tiempo, se necesita $T(x, 0)$ en todo el dominio. En
el proyecto: la temperatura de colocación de cada vertida al momento de colocarla.
*Mostrar:* argumento `T_inicial` de `resolver()`.

**B4. Dirichlet, Neumann y Robin**
- Dirichlet: se impone la **temperatura** en el borde. Ej.: verificación; medio a temperatura fija.
- Neumann: se impone el **flujo**. $q = 0$ = aislado perfecto o plano de simetría; $q > 0$ = radiación.
- Robin: el flujo depende de la diferencia con el ambiente: $q = (T_\infty - T_s)/(1/h + R_{extra})$.
  Superficie al aire; moldaje y aislación son resistencias en serie ($R = e/k$).
*Mostrar:* clases `Dirichlet`, `Neumann`, `Robin` y `_conductancia_borde()`.

**B5. ¿Qué es h y de qué depende?** Coeficiente de convección [W/(m²K)]; crece con el viento
(orden 5–25). Se definirá con datos meteorológicos en la semana 4.

**B6. ¿Qué efecto tiene la aislación?** Aumenta la resistencia: la superficie queda más caliente
y el diferencial baja (fig04: de ≈38 °C a ≈12 °C con 5 cm de poliestireno, fuente ilustrativa).
Puede subir $T_{max}$; por eso hay que mirar ambos criterios.

---

### Bloque C — Discretización y estabilidad

**C1. ¿Cómo discretizaron en 1D?**
Volúmenes finitos centrados en celdas: $N$ celdas de ancho $\Delta x$, una temperatura por
celda. Flujo entre vecinas $F = G(T_i - T_{i+1})$ con $G = 1/(\Delta x/2k_i + \Delta x/2k_{i+1})$.
Con $k$ uniforme se recupera $\alpha(T_{i+1} - 2T_i + T_{i-1})/\Delta x^2$. En el tiempo, Euler
explícito. *Mostrar:* `Malla1D`, `conductancia_interior`, bucle de `resolver()`.

**C2. ¿Por qué celdas y no nodos?** Porque en 3D cada vóxel será una celda, el balance conserva
energía exactamente y los materiales distintos/moldajes se tratan como resistencias en serie.

**C3. ¿Qué es el número de Fourier?** $Fo = \alpha\Delta t/\Delta x^2$: relación entre el paso
de tiempo y el tiempo que tarda el calor en cruzar una celda.

**C4. Estabilidad del esquema explícito**
La actualización es $T_i^{n+1} = (1 - \Delta t\,S_i/C_i)\,T_i^n + \dots$ Si el coeficiente es
negativo, la celda responde al revés y aparece un serrucho que crece. Condición:
$\Delta t \le C_i/S_i$ en todas las celdas → $Fo \le 1/2$ en el interior (1D), $1/4$ en 2D,
$1/6$ en 3D. Von Neumann: el modo serrucho se multiplica por $g = 1 - 4Fo$; si $Fo > 1/2$,
$|g| > 1$. *Mostrar:* `paso_tiempo_maximo()`, `fig05_estabilidad.png`, `test_esquema_inestable_explota`.

**C5. ¿Por qué en el borde Dirichlet el límite es Fo ≤ 1/3?** Porque la superficie está a
$\Delta x/2$ del centro de la celda: su conductancia es $2k/\Delta x$ y la suma de
conductancias de esa celda es $3k/\Delta x$. *Mostrar:* `test_paso_maximo_borde_dirichlet_es_fourier_un_tercio`.

**C6. ¿Qué pasa si piden un dt inestable?** `resolver()` lanza `ValueError` explicando el máximo.
Por defecto usa el 90 % del paso máximo.

**C7. ¿Estable significa correcto?** No. Estable = los errores no crecen. Que sea correcto se
demuestra con la verificación.

---

### Bloque D — Verificación

**D1. ¿Cómo saben que el solver está bien?** 15 tests contra respuestas conocidas: solución
exacta del seno, conservación de energía, balance con flujo y fuente, aumento adiabático lineal,
estados estacionarios con Robin y con dos materiales, límites de estabilidad.
*Mostrar:* `pytest -v`.

**D2. ¿Qué significa orden 2 en el espacio?** Al dividir $\Delta x$ por 2 el error se divide por
$2^2 = 4$. Lo medimos: orden 2,00 (fig06). En el tiempo: orden 1 (Euler explícito).

**D3. ¿Verificar vs validar?** Verificar: el programa resuelve bien la ecuación. Validar: la
ecuación representa la realidad (comparando con mediciones de Ju y Lei, más adelante).

---

### Recorrido sugerido para la demo (≈ 5 min)

1. `README.md`: estructura del repo y cómo instalar.
2. `config/` → `caso_base.yaml` y `plan_v02.yaml` (geometría docente, propiedades provisionales y cinco vertidas).
3. `python scripts/01_crear_modelo_ifc.py` y `python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc`.
4. `mcplan/solver1d.py`: docstring del módulo → `Robin` → `paso_tiempo_maximo` → bucle de `resolver`.
5. `pytest -v`.
6. Figuras fig01 (exacta), fig05 (estabilidad), fig06 (convergencia), fig04 (física del problema).


---

## E04 — 2D/3D, condiciones de borde y fuentes (miércoles 28 de octubre)

> Objetivo del enunciado: *demostrar dominio del solver térmico general antes de introducir
> hidratación y clima real.* Teoría: [Parte D](1_teoria.md#parte-d-solver-2d3d-mismo-método-en-tres-direcciones).
> Código: `mcplan/solver3d.py` (los números de línea de abajo son de ese archivo).

### Checklist individual (antes del miércoles)

- [ ] `pytest -v` → **43 passed**.
- [ ] `python scripts/06_verificacion_3d.py` y mirar fig07 a fig10.
- [ ] Poder escribir **en la pizarra** el balance de un voxel y el criterio $Fo \le 1/6$.
- [ ] Abrir `fundaciones/fundacion.ifc` en Bonsai y mostrar las 5 vertidas con sus propiedades.
- [ ] Responder las 9 preguntas de abajo sin leer.

### Lo que el repo debe mostrar (según el enunciado)

| Pide | Dónde mostrarlo |
|---|---|
| 2D | `Malla3D.uniforme(..., (n, n, 1), ...)`, fig07, `test_solucion_exacta_seno[2]` |
| 3D | `mcplan/solver3d.py`, fig08–fig10, `test_solucion_exacta_seno[3]` |
| Dirichlet / Neumann / Robin | las mismas clases del 1D, una por cara: `bordes={"x-": ..., "z+": ...}` |
| Fuente volumétrica | argumento `fuente=` (línea 340), `test_fuente_uniforme_adiabatica_aumento_lineal` |
| Modelo BIM del caso base | `fundaciones/fundacion.ifc` en Bonsai + `python scripts/02_inspeccionar_ifc.py` |
| Método de Schmidt | [teoría D.11](1_teoria.md), `test_en_el_limite_es_el_metodo_de_schmidt` |

### Las 9 preguntas posibles

**E1. Escribir el balance de un voxel.**
Acumulación = lo que entra − lo que sale + lo que se genera, en cada una de las 3 direcciones:

$$
\rho c_p V\,\frac{dT}{dt} = (F_{x,atrás} - F_{x,adel}) + (F_{y,atrás} - F_{y,adel}) + (F_{z,abajo} - F_{z,arriba}) + Q\,V
$$

con $F = G\,(T_{celda} - T_{vecina})$ y $G = kA/d$ (mismo material). Con $k$ uniforme y celdas
cúbicas, dividiendo por $\rho c_p V$ queda la clásica de 7 puntos:
$\dfrac{dT}{dt} = \dfrac{\alpha}{\Delta x^2}\left(\sum_{6\ vecinas} T_v - 6T\right) + \dfrac{Q}{\rho c_p}$.
*Mostrar:* líneas 321–341 (`balance += ...` y `T = T + dt * (balance + Q * V) / C`).

**E2. Derivar el criterio de estabilidad.**
Euler explícito con la forma anterior:
$T^{n+1} = (1 - 6Fo)\,T^n + Fo\sum_{6} T_v^n$. La nueva temperatura es una mezcla de la propia y
de las vecinas. Si el peso de la propia, $(1 - 6Fo)$, fuera **negativo**, la celda respondería al
revés (más caliente ahora → más fría después, y se pasa): aparece un serrucho que crece. Exigimos
$1 - 6Fo \ge 0$ → $Fo \le 1/6$. En 2D son 4 vecinas → $1/4$; en 1D → $1/2$. En general:
$\Delta t \le C/S$ (capacidad / suma de conductancias), que también cubre bordes y materiales distintos.
*Mostrar:* `paso_tiempo_maximo()` (línea 198), `test_paso_maximo_es_fourier_un_cuarto_en_2d_y_un_sexto_en_3d`.

**E3. Explicar flujo y signos.**
$F$ es positivo hacia $+x$, $+y$, $+z$, y $F = G\,(T_{atrás} - T_{adelante})$: el calor va de
caliente a frío (el signo menos de Fourier). En los bordes, $q_{entrante} > 0$ significa que entra
calor. Por la cara de inicio eso es $F = +q$; por la cara final es $F = -q$.
*Mostrar:* líneas 326, 329–335, 337. Teoría D.7.

**E4. Identificar una superficie adiabática.**
Es una superficie por donde no cruza calor: `Neumann(0)`, $\partial T/\partial n = 0$. Puede ser
una aislación ideal o un **plano de simetría** (permite simular 1/4 de la fundación). En el código,
toda cara no nombrada en `bordes` es adiabática. En un gráfico, las isotermas llegan perpendiculares
a ella (la base en fig09). *Mostrar:* `_completar_bordes()` (línea 174). Teoría D.9.

**E5. Implementar conceptualmente la convección.**
Ley de Newton: $q = h\,(T_{aire} - T_{sup})$. En una malla de celdas, la superficie no es una
incógnita, así que se arma un **circuito de resistencias en serie** desde el centro de la celda
hasta el aire: media celda $\Delta x/2k$ + moldaje/aislación $R_{extra}$ + película de aire $1/h$.
$G_b = A / (\Delta x/2k + R_{extra} + 1/h)$ y $q = G_b\,(T_{aire} - T_{celda})$. Error típico:
imponer $T_{aire}$ como Dirichlet en la superficie (equivale a $h = \infty$).
*Mostrar:* clase `Robin` y `_conductancia_borde()` en `solver1d.py`; `test_estacionario_robin_con_moldaje`.

**E6. Distinguir fuente volumétrica de flujo superficial.**
$Q$ [W/m³] se genera **dentro** de todo el volumen (hidratación): suma $Q\,V$ en cada celda. $q$
[W/m²] entra **por una cara** (sol): suma $q\,A$ solo en las celdas del borde. Tienen unidades
distintas y van en lugares distintos del balance. La convección no es fuente: es intercambio por
el borde. *Mostrar:* línea 340 vs. `Neumann(q)`; los dos tests (fuente pareja vs. cara que se
calienta primero). Teoría D.8.

**E7. Explicar cómo tratar dos materiales.**
Cada celda tiene su $k$, $\rho$, $c_p$. Entre dos celdas distintas el calor cruza dos medias celdas
**en serie**: $G = A/(d/2k_1 + d/2k_2)$, que equivale a la media armónica
$k = 2k_1k_2/(k_1+k_2)$. No se usa el promedio porque una capa aislante controla todo el paso. Así
se conserva el flujo en la interfaz. *Mostrar:* `Malla3D.por_materiales()` (línea 112),
`conductancia_interior()`, `test_interfaz_hormigon_sobre_suelo` (quiebre del perfil en la interfaz).

**E8. Explicar por qué las condiciones de borde pueden dominar el resultado.**
*Mostrar la fig10:* el mismo hormigón y la misma fuente; solo con cambiar el borde, el diferencial
núcleo–superficie pasa de 41 °C (desnudo) a 11 °C (5 cm de aislación) y a 0 °C (todo aislado),
mientras el núcleo casi no cambia (73–78 °C). La razón: el calor penetra solo ~1 m en
días ($\sqrt{\alpha t}$), así que el núcleo se comporta casi adiabático y la superficie la decide
el borde. Además $Bi = hL/k \approx 5 \gg 1$: el aire saca calor más rápido de lo que el interior
lo repone. Por eso moldaje, aislación y retiro son decisiones del plan. Teoría D.10.

**E9. Mostrar el avance BIM de la fundación de referencia.**
Geometría del caso base: planta **octagonal** de 18 m entre caras, borde de 1,10 m, 3,30 m al
centro y pedestal Ø7 m × 1 m (altura total 4,30 m, 636,8 m³), dividida en **5 vertidas**
(`IfcFooting`) con su pset `UANDES_HormigonMasivo` (ID, hormigón, inicio, temperatura de colocación,
moldaje, aislación, retiro). Se genera por script desde `config/` y se edita en Bonsai; Python
detecta cambios de propiedades y de geometría. **Próximo paso (semana 3):** voxelizar este IFC
para que cada celda del solver sepa si es hormigón (y de qué vertida), suelo o aire.
*Mostrar:* Bonsai (File → Open IFC Project) y `python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc`.

### Si preguntan por el método de Schmidt

Es nuestro mismo esquema explícito, usado justo en el límite: con $Fo = 1/2$ (1D), $1/4$ (2D) o
$1/6$ (3D), el peso de la propia celda se anula y **la nueva temperatura es el promedio de las
vecinas**. Se calculaba a mano. Nosotros usamos el 90 % del límite (con margen) y Robin y
materiales generales. *Mostrar:* `test_en_el_limite_es_el_metodo_de_schmidt`. Teoría D.11.

### Recorrido sugerido para la demo (≈ 5 min)

1. Bonsai: la fundación octagonal con sus 5 vertidas y propiedades (E9).
2. `mcplan/solver3d.py`: docstring (balance y caras) → bucle de `resolver()` (líneas 321–341).
3. `pytest tests/test_solver3d.py -v`: la batería de verificación.
4. fig07 (exacta 2D) → fig08 (orden 2) → fig10 (los bordes dominan).


---

## E05 — Lab: batería de verificación térmica (jueves 29 de octubre)

> Enunciado: *ejecutar en vivo una batería de tests 1D/2D/3D y condiciones de borde. El
> ayudante seleccionará tests al azar y puede cambiar parámetros. También se revisará
> conservación de energía.*
> Archivo del lab: [`tests/test_e05_bateria.py`](../tests/test_e05_bateria.py).

### Antes de entrar (2 minutos)

```bash
cd "C:\Users\Administrador\OneDrive\Escritorio\UANDES\202602\Metodos computacionales\P2\P2-Grupo7"
git pull
C:\Users\Administrador\.venvs\p2-grupo7\Scripts\Activate.ps1
pytest tests/test_e05_bateria.py -v -s
```

Deben salir **7 passed**, cada uno con sus números. Deja abierto `tests/test_e05_bateria.py` en el
editor, en el bloque **PARÁMETROS DEL LAB**.

### Los 7 casos: qué correr y qué decir

`-s` hace que se vean los números; sin `-s` solo sale PASSED/FAILED.

| # | Caso del enunciado | Comando | Qué comprueba | Qué mirar en los números |
|---|---|---|---|---|
| 1 | Solución 1D exacta | `pytest -k exacta_1d -s` | $T = T_0 + A\sin(\pi x/L_x)\,e^{-\alpha\lambda t}$, $\lambda = (\pi/L_x)^2$ | error < 1 % de A; λ medido ≈ λ teórico |
| 2 | Seno 2D | `pytest -k seno_2d -s` | fórmula 2D del enunciado, $\lambda = (\pi/L_x)^2 + (\pi/L_y)^2$ | igual que el 1 |
| 3 | Seno 3D | `pytest -k seno_3d -s` | fórmula 3D, $\lambda$ con tres términos | igual que el 1 |
| 4 | Energía adiabática | `pytest -k energia -s` | todas las caras aisladas: $E$ constante | variación ~$10^{-10}$ J (redondeo) |
| 5 | Fuente uniforme | `pytest -k fuente -s` | $dT/dt = q_0/(\rho c_p)$ en todas las celdas; $E = E_0 + q_0 V t$ | T final = exacta; energía ganada = $q_0Vt$ |
| 6 | Neumann | `pytest -k neumann -s` | entra $q$ por la cara de arriba: $E = E_0 + qAt$ | energía ganada = $qAt$; arriba más caliente que abajo |
| 7 | Robin | `pytest -k robin -s` | (a) estacionario = resistencias en serie; (b) el calor va hacia el aire más frío | error del perfil ~$10^{-13}$; la energía baja si el aire está más frío |

**Cómo se "mide" λ (casos 1 a 3):** el test ve cuánto quedó del seno en el resultado del solver
(su amplitud) y despeja $\lambda = -\ln(\text{amplitud}/A)/(\alpha t)$. Si da lo mismo que la
fórmula $\sum(\pi/L_i)^2$, el seno se apagó a la velocidad que dice la teoría.

**Por qué el tiempo por defecto:** si `T_FINAL = None`, se simula hasta que queda ~60 % del seno
($t = 0{,}5/(\alpha\lambda)$). Así siempre hay algo que comparar, sea cual sea el material o el dominio.

**Ojo con la notación:** en el enunciado $q_0$ es la fuente **volumétrica** [W/m³] (caso 5). El
flujo por una cara es otra cosa, [W/m²], y en el archivo se llama `q_SUP` (caso 6).

### Cuando el ayudante cambie parámetros

1. Cambiar el valor en el bloque **PARÁMETROS DEL LAB** y guardar (Ctrl+S).
2. Volver a correr el mismo comando.
3. **Antes de correr, decir qué esperas.** Eso es lo que demuestra que entiendes.

Lo que pasa en cada cambio (números reales del caso 2D):

| Cambio | Qué predecir | Resultado |
|---|---|---|
| base (N = 20, L = 1 m, A = 40) | — | error 0,0196 °C (0,049 % de A), λ = 19,74 |
| `N = 10` (celdas el doble de grandes) | orden 2: error **×4** | error 0,0775 °C (×3,95) ✅ |
| `Lx = 2.0` | su término en λ se divide por 4: $\pi^2/4 + \pi^2 = 12{,}34$; se enfría más lento | λ teórico 12,34, medido 12,36 ✅ |
| `A = 80` | la ecuación es lineal: el error se duplica, pero el % no cambia | error 0,039 °C (0,049 % de A) ✅ |
| `T0 = 0` | solo desplaza todo; nada cambia | mismos errores ✅ |
| hormigón (`k, rho, cp = 2.0, 2400.0, 900.0`) | el material solo cambia la **escala de tiempo** (α ~100 veces menor → ~100 veces más lento) | t = 27 360 s en vez de 253 s; mismo % de error ✅ |
| `DT = 10` (paso muy grande) | **el solver se niega**: es inestable | `ValueError: dt = 9.742 s supera el paso máximo estable 4.167 s` ✅ |

Al terminar, para volver a los valores originales (**no subir los cambios del lab**):

```bash
git restore tests/test_e05_bateria.py
```

### Si algo "falla" en vivo: casi siempre es correcto

| Lo que ves | Por qué pasa | Qué decir |
|---|---|---|
| `ValueError: ... supera el paso máximo estable` | `DT` mayor que $C/S$ (Fo > 1/2, 1/4 o 1/6) | "Es inestable; el solver lo detecta antes de calcular. Con `DT = None` usa el 90 % del máximo." |
| Falla el caso 3D con un `DT` que en 2D funcionaba | en 3D el límite es más estricto (6 vecinas; esquinas con 3 caras Dirichlet: $Fo \le 1/9$) | "Más vecinas → paso más chico." |
| "el seno ya se apagó: no se mide λ" | `T_FINAL` muy largo: queda menos de 0,1 % del seno | "Todo llegó a $T_0$; no queda decaimiento que medir. El error sigue comprobándose." |
| Fallan los casos 1–3 con `N` muy chico (4 o menos) | con tan pocas celdas el error supera el 1 % de A | "Es correcto que falle: la malla no resuelve el seno. Con N = 6 ya pasa." |
| Neumann o Robin: la energía cambia | por esas caras entra o sale calor | "No tiene que conservarse: tiene que cambiar exactamente lo que entra ($qAt$)." |

### Preguntas sobre conservación de energía

**¿Qué es E y cómo se calcula?** $E = \sum_i \rho c_p V\,T_i$ [J], la energía interna medida
desde 0 °C. *Mostrar:* `Resultado3D.energia()` en `mcplan/solver3d.py`.

**¿Por qué se conserva exactamente?** Por cómo está armado el método de volúmenes finitos: el
flujo por una cara interior se **resta** a una celda y se **suma** exactamente igual a su vecina.
Al sumar todas las celdas, los flujos interiores se cancelan de a pares y solo quedan los de los
bordes. Si los bordes son adiabáticos, no entra ni sale nada: $E$ no cambia.

**¿Por qué no da exactamente 0 sino ~$10^{-10}$ J?** Por el **redondeo** del computador (unos
16 dígitos). Sobre $E \approx 4\cdot10^5$ J, $10^{-10}$ J es un error relativo de
$10^{-16}$: el mínimo posible. Equivale a $10^{-14}$ °C de temperatura promedio.

**¿Se conserva con Dirichlet?** No, y no debe: por las caras a $T_0$ sale calor. Lo que se
comprueba ahí es la solución exacta. Con Neumann y con fuente se comprueba el balance completo:
energía ganada = $qAt$ o $q_0Vt$, exacto.

**¿Qué pasaría si un signo estuviera mal?** La energía aumentaría con todo aislado, o un bloque
caliente al aire frío se calentaría. Los casos 4 y 7b lo detectan.
