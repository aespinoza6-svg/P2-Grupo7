# Teoría del proyecto

Todo lo que hay que saber explicar, en el orden en que conviene estudiarlo. Cada parte dice
**qué hicimos, por qué y dónde está en el código**.

| Parte | Tema | Responde a… |
|---|---|---|
| [A](#parte-a-bim-ifc-blender-bonsai-e-ifcopenshell) | BIM, IFC, Blender, Bonsai e IfcOpenShell | ¿Qué es BIM/IFC? ¿Qué rol tienen Blender y Bonsai? ¿Por qué el solver no depende de Blender? |
| [B](#parte-b-física-de-la-conducción-de-calor) | Física de la conducción de calor | Significado de k, ρ, cp, α; balance de energía; condición inicial; Dirichlet, Neumann, Robin |
| [C](#parte-c-solver-1d-discretización-estabilidad-y-verificación) | Solver 1D (discretización, estabilidad y verificación) | ¿Cómo discretizamos? ¿Por qué es estable? ¿Cómo sabemos que está bien? |

---

## Parte A: BIM, IFC, Blender, Bonsai e IfcOpenShell

> Código: [`mcplan/bim.py`](../mcplan/bim.py) ·
> Modelo: [`modelos/fundacion_G7_plan_v01.ifc`](../modelos/fundacion_G7_plan_v01.ifc) ·
> Tests: [`tests/test_bim.py`](../tests/test_bim.py)

### A.1 Conceptos

| Concepto | Qué es | Rol en nuestro proyecto |
|---|---|---|
| **BIM** (Building Information Modeling) | Una forma de trabajar con un **modelo digital de información** del proyecto: cada objeto es un elemento constructivo con **identidad, tipo, relaciones y propiedades**, no solo líneas o geometría. | Describe **qué se construye y cuándo**: la fundación, sus vertidas, el material, inicio planificado, temperatura de colocación, moldaje, aislación, retiro. **No** describe el campo de temperaturas (eso es resultado del solver). |
| **IFC** (Industry Foundation Classes) | **Formato abierto** e internacional (ISO 16739, buildingSMART) para intercambiar modelos BIM. Es un archivo de texto `.ifc`. | Es el **contrato** entre el modelo y el cálculo: Bonsai lo escribe, nuestro Python lo lee. |
| **Blender** | Programa de modelado 3D de código abierto (mallas, render). Por sí solo no sabe qué es una zapata ni qué es IFC. | Es la plataforma gráfica. |
| **Bonsai** (ex BlenderBIM) | Complemento de Blender que lo convierte en una herramienta **BIM nativa IFC**: lo que se edita *es* el IFC. Agrega clases IFC, árbol espacial, property sets, materiales, planificación. | Herramienta para **ver y editar** el modelo (geometría y propiedades de cada vertida). |
| **IfcOpenShell** | Librería de código abierto (C++ con interfaz Python) para leer, escribir y procesar geometría IFC. Bonsai está construido sobre ella. | La usamos desde Python para **crear** el IFC desde `config/` y **leerlo** (propiedades, volúmenes, cotas). |

### A.2 Arquitectura y por qué el solver no depende de Blender

```
Bonsai / Blender  ──►  archivo IFC  ──►  IfcOpenShell (Python)  ──►  geometría + propiedades
   (editar)          (fuente de verdad)       mcplan.bim                  ──► voxelización ──► solver
```

1. **Separación de responsabilidades**: Bonsai edita, el IFC guarda, Python calcula.
2. **Reproducibilidad**: un tercero (o el cluster HPC, que no tiene pantalla) puede correr todo
   con `pip install` + Python, sin instalar Blender.
3. **El estado interno de Blender no es un formato estable** (escena, selección, deshacer,
   objetos sin exportar). El IFC sí: es un estándar abierto, de texto, versionable en git.
4. **Testeable**: `pytest` lee el IFC sin abrir ninguna interfaz gráfica.
5. **Trazabilidad**: cada vertida tiene un `GlobalId` permanente; se puede comparar qué cambió
   entre dos versiones del modelo (`scripts/02_inspeccionar_ifc.py --comparar`).

Regla del código: **ningún módulo de `mcplan` importa `bpy` (Blender)**.

### A.3 Cómo es un archivo IFC por dentro

Es texto en formato STEP: cada línea es una entidad con un número `#id` y sus atributos.
Extracto de nuestro modelo:

```
#32=IFCRECTANGLEPROFILEDEF(.AREA.,$,$,18.,18.);                      perfil 18 x 18 m
#33=IFCFOOTING('1WBX6NOJnMUuIOfTD8k_Tl',$,'V01 - losa',...,.PAD_FOOTING.);   la vertida V01
#40=IFCEXTRUDEDAREASOLID(#32,#38,#39,1.15);                          perfil extruido 1,15 m
#50=IFCRELDEFINESBYPROPERTIES('2OBW...',$,$,$,(#33),#49);             V01 ↔ su property set
#51=IFCPROPERTYSINGLEVALUE('VertidaID',$,IFCLABEL('V01'),$);
#54=IFCPROPERTYSINGLEVALUE('TemperaturaColocacion_C',$,IFCREAL(20.),$);
```

Conceptos que aparecen:

- **Entidades (clases)**: `IfcProject`, `IfcSite`, `IfcBuilding`, `IfcBuildingStorey`, `IfcFooting`, `IfcMaterial`…
- **GlobalId**: identificador único de 22 caracteres. Es la **identidad** del objeto: no cambia
  aunque se edite su nombre, forma o propiedades.
- **Relaciones** (también son entidades): `IfcRelAggregates` (proyecto → sitio → edificio → nivel),
  `IfcRelContainedInSpatialStructure` (el nivel contiene las zapatas),
  `IfcRelDefinesByProperties` (objeto ↔ property set), `IfcRelAssociatesMaterial` (objeto ↔ material).
- **Property sets (psets)**: grupos de propiedades nombre–valor. Los que empiezan con `Pset_`
  son **estándar** de buildingSMART; los propios deben usar otro prefijo (nosotros: `MC_`).
- **Geometría**: la forma de cada vertida es un **perfil 2D extruido** (`IfcExtrudedAreaSolid`):
  rectángulo para la losa, círculo para el pedestal.
- **Unidades**: declaradas en el archivo (`IfcUnitAssignment`). Usamos **metros**
  (IfcOpenShell por defecto usaría milímetros: una losa de "18" sería de 18 mm).

### A.4 Nuestro modelo

```
IfcProject  "P2 Grupo 7 - Fundación aerogenerador"
└── IfcSite  "Renaico - Parque Eólico Vergara"   lat -37.6921, lon -72.624   pset MC_Sitio
    └── IfcBuilding  "Aerogenerador"
        └── IfcBuildingStorey  "Nivel fundación"
            ├── IfcFooting  "V01 - losa"       z = 0,00 → 1,15 m   372,6 m³
            ├── IfcFooting  "V02 - losa"       z = 1,15 → 2,30 m   372,6 m³
            └── IfcFooting  "V03 - pedestal"   z = 2,30 → 3,30 m    38,5 m³
```

**Cada vertida es una `IfcFooting`** (clase IFC para fundaciones) con:

| Propiedad (pset `MC_Vertida`) | Tipo | Significado |
|---|---|---|
| `VertidaID` | texto | Identificador del plan (V01, V02…) |
| `TipoHormigon` | texto | Tipo de hormigón (G30 provisional) |
| `InicioPlanificado` | texto ISO 8601 | Fecha y hora de inicio = `fecha_inicio` del plan + `inicio_h` |
| `TemperaturaColocacion_C` | real [°C] | Temperatura del hormigón al colocarlo (condición inicial) |
| `Moldaje` | texto | Tipo de moldaje |
| `Aislacion` | texto | Aislación de la superficie |
| `RetiroMoldaje_h` | real [h] | Horas desde el inicio de esa vertida hasta retirar el moldaje |

El **material** `Hormigon G30` lleva sus propiedades térmicas en los psets **estándar**
`Pset_MaterialCommon.MassDensity` y `Pset_MaterialThermal.ThermalConductivity /
SpecificHeatCapacity`, para que cualquier software BIM las entienda.

#### A.4.1 De dónde sale el modelo

El IFC **no se dibuja a mano**: lo genera `scripts/01_crear_modelo_ifc.py` desde dos archivos:

| Archivo | Contiene |
|---|---|
| `config/caso_base.yaml` | sitio (coordenadas oficiales del curso), geometría y hormigón (supuestos marcados `[PROVISIONAL]`) |
| `config/plan_v01.yaml` | las vertidas: cotas, inicio, temperatura, moldaje… (**ejemplo, no recomendación**) |

Ventajas: cambiar una decisión (p. ej. la fecha de inicio o el número de capas) es cambiar un
YAML y regenerar, sin tocar código. `validar_plan()` rechaza planes con huecos o traslapes.

#### A.4.2 Identidades estables

IfcOpenShell inventa un `GlobalId` aleatorio para cada objeto nuevo. Si lo dejáramos así, cada
vez que se regenera el IFC la "misma" vertida tendría otra identidad. Por eso `_guid_estable()`
calcula el `GlobalId` a partir de una clave fija (p. ej. `"vertida/V01"`) con `uuid5`:
**regenerar produce exactamente el mismo archivo** (`test_globalid_estable_al_regenerar`).

### A.5 Abrir y editar el modelo en Bonsai

> Los nombres de menús pueden variar levemente entre versiones de Bonsai; ver también los
> tutoriales del curso (Prof. Rino).

1. Instalar Blender (4.2 o superior) y luego Bonsai desde *Edit → Preferences → Get Extensions*
   (buscar "Bonsai").
2. Abrir `modelos/fundacion_G7_plan_v01.ifc` con *File → Open IFC Project*.
3. En el árbol espacial (panel de Bonsai) aparecen Sitio → Edificio → Nivel → 3 `IfcFooting`.
4. Seleccionar una vertida → pestaña de propiedades del objeto → panel de *Property Sets*:
   ahí está `MC_Vertida`. Editar, por ejemplo, `TemperaturaColocacion_C`.
5. Guardar **como una versión nueva** (*Save IFC Project As…*), p. ej.
   `modelos/fundacion_G7_plan_v01_bonsai.ifc`. No sobrescribir el `v01` generado.
6. Volver a Python y comparar:

```bash
python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc --comparar modelos/fundacion_G7_plan_v01_bonsai.ifc
```

Salida esperada (si se cambió la temperatura de V02):

```
DIFERENCIAS  fundacion_G7_plan_v01.ifc  ->  fundacion_G7_plan_v01_bonsai.ifc
  MODIFICADA V02.TemperaturaColocacion_C: 20.0 -> 15.0
```

Este es el **round-trip IFC** (Bonsai → IFC → Python detecta el cambio) que se demuestra en el
laboratorio del jueves (E02). El test `test_cambio_de_propiedad_se_detecta` hace lo mismo
automáticamente (editando con IfcOpenShell, que es lo que Bonsai usa por dentro).

### A.6 Pendiente

- Reemplazar los supuestos `[PROVISIONAL]` por la geometría oficial del profesor
  ([3_supuestos.md](3_supuestos.md)).
- Practicar en Bonsai (cada integrante): abrir, editar una propiedad, guardar, comparar.
- Semana 3: voxelizar la geometría del IFC (`geometria_vertidas()` ya entrega volúmenes y cotas).

---

## Parte B: Física de la conducción de calor

> Objetivo: entender **qué ecuación resolvemos y por qué**, antes de ver cómo se discretiza
> (Parte C).

### B.1 El problema en una frase

El hormigón **genera calor** al hidratarse el cemento. Ese calor viaja por **conducción**
hacia las superficies y desde ahí sale al **aire** (convección) y al **suelo**. Si el calor se
genera más rápido de lo que sale, el núcleo se calienta (riesgo **T_max > 70 °C**) mientras la
superficie se enfría con el ambiente (riesgo **ΔT núcleo–superficie > 20 °C**, que produce
tensiones y fisuras).

### B.2 Propiedades térmicas

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

#### B.2.1 Longitud de difusión: por qué el hormigón masivo se calienta tanto

En un tiempo $t$, un cambio de temperatura penetra aproximadamente $\ell \approx \sqrt{\alpha t}$:

| tiempo | 1 hora | 1 día | 7 días | 14 días |
|---|---|---|---|---|
| $\sqrt{\alpha t}$ (hormigón) | 6 cm | 28 cm | 75 cm | **1,06 m** |

En los 14 días del análisis el calor "viaja" del orden de **1 m**. El centro de una fundación de
3,30 m de espesor casi no se entera de que existen las superficies: se comporta casi como si
estuviera **aislado (adiabático)** y acumula casi todo el calor que genera.
Al revés, el tiempo para que el calor cruce media fundación (1,65 m) es
$t \approx L^2/\alpha \approx 34$ días.

#### B.2.2 Número de Biot: ¿basta una sola temperatura?

$$Bi = \frac{h\,L}{k} = \frac{\text{resistencia interna (conducción)}}{\text{resistencia externa (convección)}}$$

Con $h \approx 10$ W/(m²K), $L = 1{,}65$ m y $k = 2$: $Bi \approx 8$.
Si $Bi \ll 0{,}1$ el cuerpo tendría temperatura casi uniforme; con $Bi \gg 1$ hay **grandes
gradientes internos**. Por eso necesitamos un modelo con variación espacial (1D → 3D), y no
un modelo de "una sola temperatura".

### B.3 Balance de energía

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

En nuestro código el balance se aplica **a cada celda** (ver Parte C):
es exactamente esta ecuación, escrita para volúmenes pequeños.

### B.4 Condición inicial

La ecuación es de **primer orden en el tiempo**: necesita conocer $T(x, t=0)$ en todo el dominio.

- En el proyecto: cada vertida empieza a existir cuando se coloca, con su **temperatura de
  colocación** (propiedad `TemperaturaColocacion_C` del IFC). Antes de eso, esas celdas no son
  hormigón activo (activación de vertidas, semana 3).
- El suelo parte con su propia temperatura inicial (semana 4).
- En el código: argumento `T_inicial` de `resolver()` (un número o un valor por celda).

### B.5 Condiciones de borde

La ecuación es de **segundo orden en el espacio**: necesita una condición en **cada** borde.

| Tipo | Qué se impone | Ecuación (en el borde) | Ejemplo en el proyecto | Clase en el código |
|---|---|---|---|---|
| **Dirichlet** (1ª especie) | la **temperatura** | $T = T_b(t)$ | contacto con un medio de temperatura conocida; casos de verificación | `Dirichlet(T)` |
| **Neumann** (2ª especie) | el **flujo** de calor | $q_{\text{entrante}} = q(t)$ | $q = 0$: aislación perfecta o **plano de simetría**; radiación solar absorbida | `Neumann(q)` |
| **Robin** (3ª especie, mixta) | relación **flujo–temperatura** | $q_{\text{entrante}} = h\,(T_\infty - T_s)$ | superficie expuesta al aire, con o sin moldaje/aislación | `Robin(h, T_inf, R_extra)` |

Convención de signo del código: **$q$ positivo = calor que entra al hormigón.**

#### B.5.1 Robin con moldaje o aislación: resistencias en serie

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
  cambian $T_{max}$ (el núcleo es casi adiabático, B.2.1).

#### B.5.2 Onda térmica diaria

Si el aire oscila con período $P = 24$ h, la oscilación se amortigua con la profundidad como
$e^{-z/\delta}$, con $\delta = \sqrt{2\alpha/\omega} = \sqrt{\alpha P/\pi} \approx 16$ cm.
El ciclo día-noche afecta solo los primeros centímetros (justo donde se mide la superficie, a
25 mm) y no al núcleo ([fig03](../resultados/figuras/fig03_robin_ambiente_diario.png)).

### B.6 Resumen para el oral

1. Ecuación: $\rho c_p\,\partial_t T = \nabla\cdot(k\nabla T) + Q$, que es un **balance de energía**.
2. $k$ conduce, $\rho c_p$ almacena, $\alpha$ = rapidez de propagación ($\sqrt{\alpha t}$ ≈ 1 m en 14 días).
3. Condición inicial = temperatura de colocación; condiciones de borde = qué pasa en las superficies.
4. Dirichlet fija $T$, Neumann fija el flujo, Robin relaciona flujo con $T$ (aire, moldaje, aislación).

---

## Parte C: Solver 1D (discretización, estabilidad y verificación)

> Código: [`mcplan/solver1d.py`](../mcplan/solver1d.py) ·
> Tests: [`tests/test_solver1d.py`](../tests/test_solver1d.py) ·
> Física: Parte B

### C.1 Idea: un balance de energía por celda

El enunciado pide una "grilla cartesiana regular con balance térmico por celda". Eso es un
esquema de **volúmenes finitos centrados en celdas**: dividimos el dominio en celdas y a cada
una le aplicamos el balance *acumulación = entra − sale + genera*.

¿Por qué celdas (y no nodos en los bordes, como el código de ejemplo del curso)?

- En 3D cada **vóxel será una celda**: el solver 1D ya tiene la misma lógica que tendrá el 3D.
- **Conserva la energía exactamente**: lo que sale de una celda entra a su vecina.
- Materiales distintos (hormigón–suelo) y moldaje/aislación se tratan con **resistencias en serie**.

### C.2 Malla

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

### C.3 Balance de la celda $i$

$$C_i \frac{dT_i}{dt} = F_{i} - F_{i+1} + Q_i\,\Delta x$$

| Símbolo | Significado | Unidades | Código |
|---|---|---|---|
| $C_i = \rho_i c_{p,i} \Delta x$ | capacidad térmica de la celda | J/(m²K) | `malla.capacidad` |
| $F_j$ | flujo por la cara $j$, **positivo hacia $+x$** | W/m² | arreglo `F` (N+1 caras) |
| $F_i - F_{i+1}$ | lo que entra por la izquierda menos lo que sale por la derecha | W/m² | `F[:-1] - F[1:]` |
| $Q_i \Delta x$ | calor generado en la celda | W/m² | `Q * malla.dx` |

### C.4 Flujo entre celdas: conductancias

$$F_{i+1} = G_{i+1}\,(T_i - T_{i+1}), \qquad G = \frac{1}{R}, \qquad R = \frac{\Delta x}{2k_i} + \frac{\Delta x}{2k_{i+1}}$$

El calor va del centro de la celda $i$ a su cara (media celda) y de la cara al centro de la
celda $i+1$ (otra media celda): **dos resistencias en serie**. Código: `malla.conductancia_interior`.

- Si $k$ es igual en ambas: $G = k/\Delta x$ y se recupera la diferencia finita clásica:
  $\dfrac{dT_i}{dt} = \alpha\dfrac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2}$.
- Si los materiales son distintos, la resistencia en serie (media armónica de $k$) garantiza
  que el flujo sea continuo en la interfaz. Un promedio aritmético de $k$ sería incorrecto
  (por ejemplo, entre hormigón y aislación). Lo verifica `test_estacionario_dos_materiales_en_serie`.

### C.5 Bordes

La celda del borde está a $\Delta x/2$ de la superficie. Su conductancia hacia "afuera":

| Borde | Conductancia $G_b$ | Flujo que entra |
|---|---|---|
| Dirichlet $T_b$ | $\dfrac{1}{\Delta x/2k} = \dfrac{2k}{\Delta x}$ | $G_b\,(T_b - T_0)$ |
| Robin $(h, T_\infty, R_{extra})$ | $\dfrac{1}{\Delta x/2k + R_{extra} + 1/h}$ | $G_b\,(T_\infty - T_0)$ |
| Neumann $q$ | 0 (no depende de $T$) | $q$ |

Código: `_conductancia_borde()` y `_flujo_entrante()`. Como $F$ es positivo hacia $+x$:
por la izquierda `F[0] = +q_entrante` y por la derecha `F[-1] = -q_entrante`.

### C.6 Avance en el tiempo: Euler explícito

$$T_i^{n+1} = T_i^n + \frac{\Delta t}{C_i}\left(F_i^n - F_{i+1}^n + Q_i^n\Delta x\right)$$

**Explícito** = la temperatura nueva se calcula directamente con valores del paso anterior (no
hay que resolver un sistema de ecuaciones). Es simple y fácil de explicar, pero tiene un límite
de estabilidad para $\Delta t$. Precisión: orden 1 en el tiempo y orden 2 en el espacio (C.8).

En el código es una sola línea dentro del bucle de `resolver()`:
`T = T + dt * (F[:-1] - F[1:] + Q * malla.dx) / C`.

### C.7 Estabilidad

#### C.7.1 Criterio de coeficientes positivos (el que usa el código)

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

#### C.7.2 Análisis de von Neumann (la versión clásica)

Para una malla uniforme infinita, un error con forma $e^{\mathrm{i}\theta j}$ se multiplica en cada
paso por el factor de amplificación

$$g(\theta) = 1 - 4\,Fo\,\sin^2(\theta/2)$$

Se necesita $|g| \le 1$ para todo $\theta$. El peor caso es $\theta = \pi$ (el serrucho
$+,-,+,-$), con $g = 1 - 4Fo$: estable si $Fo \le 1/2$. Con $Fo = 0{,}55$, $g = -1{,}2$: el
serrucho **cambia de signo y crece un 20 % en cada paso**
([fig05](../resultados/figuras/fig05_estabilidad.png)).

#### C.7.3 Números para el hormigón ($\alpha \approx 9{,}3\cdot10^{-7}$ m²/s)

- 1D, $\Delta x = 5$ cm, borde Dirichlet ($Fo \le 1/3$): $\Delta t_{max} \approx 900$ s = 15 min.
- 3D, $\Delta x = 10$ cm ($Fo \le 1/6$): $\Delta t_{max} \approx 1800$ s = 30 min → 14 días ≈ 700 pasos.

**Estable no significa preciso**: la estabilidad solo garantiza que los errores no crezcan. La
precisión se mide con la verificación.

### C.8 Verificación

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

### C.9 Cómo usarlo

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

### C.10 Limitaciones actuales (próximas semanas)

- Es 1D: la extensión 2D/3D (semana 2) reutiliza la misma idea de celdas y conductancias.
- La fuente $Q$ no depende de la temperatura; la hidratación real sí (semana 3).
- $R_{extra}$ es constante; el retiro de moldaje requiere que cambie en el tiempo (semana 4).
- No hay activación de vertidas ni suelo explícito (semanas 3–4).
