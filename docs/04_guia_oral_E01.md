# 04 — Guía de estudio para la revisión oral E01 (miércoles 21 de octubre)

> En la revisión **se puede preguntar a cualquier integrante** y se pide **navegar por código
> real**. Cada respuesta indica qué mostrar. Las respuestas son guías: hay que poder decirlas
> con palabras propias.
>
> Lecturas previas exigidas: guías del curso G00, G01, G02 y G05.

## Checklist individual (cada integrante, antes del miércoles)

- [ ] Clonar el repo, crear el entorno e instalar (`README.md` → *Instalación*).
- [ ] Correr `pytest` y ver los 23 tests pasar.
- [ ] Correr los scripts 01 a 04 y mirar cada figura en `resultados/figuras/`.
- [ ] Abrir el IFC en Bonsai, cambiar una propiedad, guardar como otra versión y compararla
      con `scripts/02_inspeccionar_ifc.py --comparar`.
- [ ] Cambiar `dt` en un ejemplo para que `Fo > 1/2` y ver el error / la explosión.
- [ ] Poder responder **todas** las preguntas de abajo sin leer.

---

## Bloque A — BIM / IFC

**A1. ¿Qué es BIM?**
Una forma de trabajar con un modelo digital de **información** del proyecto: cada objeto es un
elemento constructivo con identidad, tipo, relaciones y propiedades (no solo un dibujo). En
nuestro caso describe qué se construye y cuándo: la fundación, cada vertida, su material, inicio
planificado, temperatura de colocación, moldaje y aislación.
*Mostrar:* `modelos/fundacion_G7_plan_v01.ifc` en Bonsai o la salida de `scripts/02_inspeccionar_ifc.py`.

**A2. ¿Qué es IFC?**
El formato **abierto** e internacional (ISO 16739, buildingSMART) para intercambiar modelos BIM.
Es un archivo de texto con entidades (`IfcFooting`, `IfcSite`…), relaciones y property sets.
Al ser abierto, cualquier programa lo lee: Bonsai, Revit o nuestro Python.
*Mostrar:* abrir el `.ifc` en un editor de texto y ubicar la línea `IFCFOOTING(...'V01 - losa'...)`.

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
*Mostrar:* `src/mcplan/__init__.py` (regla: ningún módulo importa Blender) y `bim.leer_vertidas()`.

**A6. ¿Dónde están las propiedades de cada vertida?** En el property set `MC_Vertida` de cada
`IfcFooting` (`VertidaID`, `TipoHormigon`, `InicioPlanificado`, `TemperaturaColocacion_C`,
`Moldaje`, `Aislacion`, `RetiroMoldaje_h`). Las del material en los psets estándar
`Pset_MaterialCommon` y `Pset_MaterialThermal`.
*Mostrar:* `bim.crear_modelo_fundacion()`, paso 5.

**A7. ¿Por qué `MC_Vertida` y no `Pset_Vertida`?** El prefijo `Pset_` está reservado para los
property sets estándar de buildingSMART; los propios deben usar otro prefijo.

**A8. ¿Qué es el GlobalId y por qué lo fijamos?** Es la identidad única y permanente de cada
objeto IFC. IfcOpenShell lo genera al azar; nosotros lo calculamos con `uuid5` a partir de una
clave (`"vertida/V01"`) para que regenerar el modelo no cambie las identidades.
*Mostrar:* `bim._guid_estable()` y `test_globalid_estable_al_regenerar`.

**A9. ¿Qué dimensiones son datos y cuáles supuestos?** Datos del enunciado: 18 × 18 m, 3,30 m,
Ø 7 m. Supuestos provisionales (documentados, no silenciosos): 3,30 m como altura total y
pedestal de 1,00 m. *Mostrar:* `config/fundacion.yaml` y [05](05_decisiones_y_supuestos.md).

---

## Bloque B — Física

**B1. Significado físico de k, ρ, c_p, α**
- $k$ [W/(m·K)]: qué tan fácil **conduce** calor (flujo por unidad de gradiente).
- $\rho$ [kg/m³]: masa por volumen.
- $c_p$ [J/(kg·K)]: energía para subir 1 K a 1 kg; mide cuánto calor **almacena**.
- $\alpha = k/(\rho c_p)$ [m²/s]: **qué tan rápido se propaga** un cambio de temperatura.
  Para el hormigón $\alpha \approx 9{,}3\cdot10^{-7}$ m²/s: en 14 días el calor penetra solo
  $\sqrt{\alpha t} \approx 1$ m, por eso el núcleo de 3,30 m se calienta casi como si estuviera aislado.
*Mostrar:* clase `Material` en `solver1d.py`; [02 §2](02_fisica_conduccion.md).

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

## Bloque C — Discretización y estabilidad

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

## Bloque D — Verificación

**D1. ¿Cómo saben que el solver está bien?** 15 tests contra respuestas conocidas: solución
exacta del seno, conservación de energía, balance con flujo y fuente, aumento adiabático lineal,
estados estacionarios con Robin y con dos materiales, límites de estabilidad.
*Mostrar:* `pytest -v`.

**D2. ¿Qué significa orden 2 en el espacio?** Al dividir $\Delta x$ por 2 el error se divide por
$2^2 = 4$. Lo medimos: orden 2,00 (fig06). En el tiempo: orden 1 (Euler explícito).

**D3. ¿Verificar vs validar?** Verificar: el programa resuelve bien la ecuación. Validar: la
ecuación representa la realidad (comparando con mediciones de Ju y Lei, más adelante).

---

## Recorrido sugerido para la demo (≈ 5 min)

1. `README.md`: estructura del repo y cómo instalar.
2. `config/` → los 3 YAML (dato vs supuesto `[PROVISIONAL]`).
3. `python scripts/01_crear_modelo_ifc.py` y `python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc`.
4. `src/mcplan/solver1d.py`: docstring del módulo → `Robin` → `paso_tiempo_maximo` → bucle de `resolver`.
5. `pytest -v` (23 tests).
6. Figuras fig01 (exacta), fig05 (estabilidad), fig06 (convergencia), fig04 (física del problema).
