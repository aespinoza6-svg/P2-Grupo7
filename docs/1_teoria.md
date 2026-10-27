# Teoría del proyecto

Todo lo que hay que saber explicar, en el orden en que conviene estudiarlo. Cada parte dice
**qué hicimos, por qué y dónde está en el código**.

| Parte | Tema | Responde a… |
|---|---|---|
| [A](#parte-a-bim-ifc-blender-bonsai-e-ifcopenshell) | BIM, IFC, Blender, Bonsai e IfcOpenShell | ¿Qué es BIM/IFC? ¿Qué rol tienen Blender y Bonsai? ¿Por qué el solver no depende de Blender? |
| [B](#parte-b-física-de-la-conducción-de-calor) | Física de la conducción de calor | Significado de k, ρ, cp, α; balance de energía; condición inicial; Dirichlet, Neumann, Robin |
| [C](#parte-c-solver-1d-discretización-estabilidad-y-verificación) | Solver 1D (discretización, estabilidad y verificación) | ¿Cómo discretizamos? ¿Por qué es estable? ¿Cómo sabemos que está bien? |
| [D](#parte-d-solver-2d3d-mismo-método-en-tres-direcciones) | Solver 2D/3D | ¿Qué cambia al pasar a 3D? ¿Por qué Fo ≤ 1/4 y 1/6? ¿Cómo se tratan las interfaces? ¿Signos, fuente vs. flujo, simetría? ¿Por qué dominan los bordes? ¿Qué es el método de Schmidt? |

---

## Parte A: BIM, IFC, Blender, Bonsai e IfcOpenShell

> Código: [`mcplan/bim.py`](../mcplan/bim.py) ·
> Modelo: [`fundaciones/fundacion.ifc`](../fundaciones/fundacion.ifc) ·
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
#32=IFCCLOSEDSHELL(...);                                              caras del prisma octagonal
#33=IFCFACETEDBREP(#32);                                              sólido IFC de la vertida
#34=IFCFOOTING('1WBX6NOJnMUuIOfTD8k_Tl',$,'Vertida_01',...);          la vertida 1
#50=IFCRELDEFINESBYPROPERTIES('2OBW...',$,$,$,(#34),#49);             V01 ↔ su property set
#51=IFCPROPERTYSINGLEVALUE('VertidaID',$,IFCINTEGER(1),$);
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
  son **estándar** de buildingSMART; el conjunto docente propio se llama `UANDES_HormigonMasivo`.
- **Geometría**: prismas octagonales, troncos octagonales y el cilindro superior se guardan
  como sólidos IFC faceteados cerrados (`IfcFacetedBrep`).
- **Unidades**: declaradas en el archivo (`IfcUnitAssignment`). Usamos **metros**
  (IfcOpenShell por defecto usaría milímetros: una losa de "18" sería de 18 mm).

### A.4 Nuestro modelo

```
IfcProject  "P2 Grupo 7 - Fundación aerogenerador"
└── IfcSite  "Renaico - Parque Eólico Vergara"   lat -37.6921, lon -72.624   pset MC_Sitio
    └── IfcBuilding  "Aerogenerador"
        └── IfcBuildingStorey  "Nivel fundación"
            ├── IfcFooting  "Vertida_01" z = 0,00 → 0,55 m   147,63 m³
            ├── IfcFooting  "Vertida_02" z = 0,55 → 1,10 m  147,63 m³
            ├── IfcFooting  "Vertida_03" z = 1,10 → 2,20 m  214,22 m³
            ├── IfcFooting  "Vertida_04" z = 2,20 → 3,30 m   88,92 m³
            └── IfcFooting  "Vertida_05" z = 3,30 → 4,30 m   38,37 m³
```

**Cada vertida es una `IfcFooting`** (clase IFC para fundaciones) con:

| Propiedad (pset `UANDES_HormigonMasivo`) | Tipo | Significado |
|---|---|---|
| `VertidaID` | entero | Número de vertida (1…5); el nombre IFC es `Vertida_01`… |
| `Elemento` | texto | Zona funcional: borde, transición o pedestal |
| `Forma` | texto | Tipo de sólido IFC generado |
| `TipoHormigon` | texto | Tipo de hormigón del ejemplo docente (M1) |
| `Inicio_h` | real [h] | Horas desde `fecha_inicio` |
| `InicioPlanificado` | texto ISO 8601 | Fecha y hora de inicio = `fecha_inicio` del plan + `inicio_h` |
| `TemperaturaColocacion_C` | real [°C] | Temperatura del hormigón al colocarlo (condición inicial) |
| `Moldaje` | texto | Tipo de moldaje |
| `Aislacion` | texto | Aislación de la superficie |
| `RetiroMoldaje_h` | real [h] | Horas desde el inicio de esa vertida hasta retirar el moldaje |

El **material** `Hormigon M1` lleva sus propiedades térmicas en los psets **estándar**
`Pset_MaterialCommon.MassDensity` y `Pset_MaterialThermal.ThermalConductivity /
SpecificHeatCapacity`, para que cualquier software BIM las entienda.

#### A.4.1 De dónde sale el modelo

El IFC **no se dibuja a mano**: lo genera `scripts/01_crear_modelo_ifc.py` desde dos archivos:

| Archivo | Contiene |
|---|---|
| `config/caso_base.yaml` | sitio, dimensiones de la simplificación geométrica docente y propiedades provisionales del hormigón |
| `config/plan_v02.yaml` | las cinco vertidas docentes: forma, cotas, inicio, temperatura y moldaje |

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
2. Abrir `fundaciones/fundacion.ifc` con *File → Open IFC Project*.
3. En el árbol espacial (panel de Bonsai) aparecen Sitio → Edificio → Nivel → 5 `IfcFooting`.
4. Seleccionar una vertida → pestaña de propiedades del objeto → panel de *Property Sets*:
   ahí está `UANDES_HormigonMasivo`. Editar, por ejemplo, `TemperaturaColocacion_C`.
5. Guardar **como una versión nueva** (*Save IFC Project As…*), p. ej.
   `fundaciones/fundacion_bonsai.ifc`. No sobrescribir el IFC generado.
6. Volver a Python y comparar:

```bash
python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc --comparar fundaciones/fundacion_bonsai.ifc
```

Salida esperada (si se cambió la temperatura de V02):

```
DIFERENCIAS  fundacion.ifc  ->  fundacion_bonsai.ifc
  MODIFICADA V02.TemperaturaColocacion_C: 18.0 -> 15.0
```

Este es el **round-trip IFC** (Bonsai → IFC → Python detecta el cambio) que se demuestra en el
laboratorio del jueves (E02). El test `test_cambio_de_propiedad_se_detecta` hace lo mismo
automáticamente (editando con IfcOpenShell, que es lo que Bonsai usa por dentro).

### A.6 Pendiente

- Confirmar si el profesor entrega un IFC geométrico oficial; si es distinto, reemplazará esta
  simplificación docente ([3_supuestos.md](3_supuestos.md)).
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

- Es 1D: la extensión 2D/3D está en la [Parte D](#parte-d-solver-2d3d-mismo-método-en-tres-direcciones).
- La fuente $Q$ no depende de la temperatura; la hidratación real sí (semana 3).
- $R_{extra}$ es constante; el retiro de moldaje requiere que cambie en el tiempo (semana 4).
- No hay activación de vertidas ni suelo explícito (semanas 3–4).

---

## Parte D: Solver 2D/3D (mismo método en tres direcciones)

Código: [`mcplan/solver3d.py`](../mcplan/solver3d.py). Tests: [`tests/test_solver3d.py`](../tests/test_solver3d.py).
Figuras: [`scripts/06_verificacion_3d.py`](../scripts/06_verificacion_3d.py).

### D.1 La idea: el 1D repetido en x, y, z

La ecuación ahora es

$$
\rho c_p \frac{\partial T}{\partial t} =
\frac{\partial}{\partial x}\!\left(k\frac{\partial T}{\partial x}\right) +
\frac{\partial}{\partial y}\!\left(k\frac{\partial T}{\partial y}\right) +
\frac{\partial}{\partial z}\!\left(k\frac{\partial T}{\partial z}\right) + Q
$$

y **no inventamos un método nuevo**: cada celda (voxel) hace el mismo balance que en 1D, pero
intercambia calor con **6 vecinas** (izquierda/derecha, adelante/atrás, abajo/arriba) en vez de 2.

```
            z+ (arriba)
             |
   x- ---- celda ---- x+        cada flecha = un flujo F = G·(T_celda - T_vecina)
             |
            z-          (y lo mismo con y- e y+, hacia adentro de la hoja)
```

### D.2 Balance de la celda

$$
C\,\frac{dT}{dt} = (F_{x,\text{atrás}} - F_{x,\text{adelante}}) + (F_{y,\text{atrás}} - F_{y,\text{adelante}}) + (F_{z,\text{abajo}} - F_{z,\text{arriba}}) + Q\,V
$$

| Símbolo | Significado | Unidad |
|---|---|---|
| $V = \Delta x\,\Delta y\,\Delta z$ | volumen de la celda | m³ |
| $C = \rho c_p V$ | capacidad térmica de la celda | J/K |
| $F$ | flujo por una cara, positivo hacia $+x$, $+y$ o $+z$ | **W** |

Diferencia con el 1D: allá trabajábamos "por m²" (W/m²) porque todas las caras tenían la misma
área. En 3D las caras $x$, $y$, $z$ pueden tener áreas distintas ($A_x = \Delta y\,\Delta z$,
$A_y = \Delta x\,\Delta z$, $A_z = \Delta x\,\Delta y$), así que usamos W y J/K.

### D.3 Conductancias e interfaces entre materiales

Entre dos celdas vecinas el calor atraviesa **dos medias celdas en serie**:

$$
G = \frac{A}{\dfrac{d}{2k_1} + \dfrac{d}{2k_2}} \qquad (d = \Delta x,\ \Delta y \text{ o } \Delta z)
$$

- Mismo material: $G = kA/d$ (la diferencia finita de siempre).
- Materiales distintos: equivale a usar la **media armónica** $k_{cara} = \dfrac{2k_1k_2}{k_1+k_2}$.
  ¿Por qué armónica y no el promedio? Porque las resistencias en serie se **suman**: si una de las
  dos capas aísla mucho, ella controla el paso de calor. El promedio aritmético lo sobrestimaría.

Esto es lo que permitirá, sin cambiar el solver, poner **hormigón sobre suelo**, **vertida nueva
sobre vertida vieja** o **aislación**: basta con decir qué material ocupa cada celda
(`Malla3D.por_materiales`). Es justo lo que entregará la voxelización del IFC en la semana 3.

### D.4 Bordes: 6 caras con nombre

| Cara | Dónde | Ejemplo en la fundación |
|---|---|---|
| `"x-"`, `"x+"` | $x = 0$, $x = L_x$ | costados (moldaje → Robin con $R_{extra}$) |
| `"y-"`, `"y+"` | $y = 0$, $y = L_y$ | costados |
| `"z-"` | $z = 0$, abajo | contacto con el suelo |
| `"z+"` | $z = L_z$, arriba | superficie expuesta al aire |

Se usan las **mismas clases del 1D** (`Dirichlet`, `Neumann`, `Robin`): la física del borde no
cambia por estar en 3D. **Una cara que no se indica queda aislada** (`Neumann(0)`), y si se escribe
mal el nombre de una cara el programa da error (no la deja aislada en silencio).

**2D** = una sola capa de celdas en $z$ ($n_z = 1$) con las caras $z$ aisladas: no fluye calor en $z$.

### D.5 Estabilidad: por qué 1/4 en 2D y 1/6 en 3D

El criterio es **el mismo del 1D** (coeficientes positivos): $\Delta t \le C/S$ en toda celda,
con $S$ = suma de las conductancias de la celda. En malla uniforme ($\Delta x = \Delta y = \Delta z$),
celda interior:

| | vecinas | $S$ | límite |
|---|---|---|---|
| 1D | 2 | $2kA/\Delta x$ | $Fo \le 1/2$ |
| 2D | 4 | $4kA/\Delta x$ | $Fo \le 1/4$ |
| 3D | 6 | $6kA/\Delta x$ | $Fo \le 1/6$ |

Intuición: la celda entrega calor a más vecinas en cada paso; si el paso es muy grande, entrega
más de lo que tiene y "se pasa" (oscila). Con más vecinas, el paso permitido es más chico.

Ojo con las **esquinas con bordes Dirichlet**: la superficie está a $\Delta x/2$ (conductancia
$2kA/\Delta x$). Una esquina 3D con tres caras Dirichlet tiene $S = 3\cdot 1 + 3\cdot 2 = 9$ en
unidades de $kA/\Delta x$, entonces $Fo \le 1/9$. El código lo calcula solo (`paso_tiempo_maximo`).

Para el hormigón ($\alpha \approx 9{,}3\cdot10^{-7}$ m²/s) con celdas de 10 cm en 3D:
$\Delta t \le (0{,}1)^2/(6\cdot 9{,}3\cdot10^{-7}) \approx 1800$ s ≈ 30 min. Para 14 días son
~700 pasos: barato. Con celdas de 5 cm el paso cae a ~7,5 min (4 veces menos) y hay 8 veces más
celdas: **32 veces más trabajo**. Esto pesará en Monte Carlo (semana 5).

### D.6 Verificación (batería de la guía del curso, sección 33)

| Test (`tests/test_solver3d.py`) | Qué comprueba | Test de la guía |
|---|---|---|
| `test_solucion_exacta_seno[2]` / `[3]` | Seno con bordes fijos vs. solución exacta $e^{-\alpha\pi^2(1/L_x^2+1/L_y^2+1/L_z^2)t}$ | 2 y 3 |
| `test_convergencia_orden_2[2]` / `[3]` | Duplicar las celdas divide el error por ~4 (orden 2) | 2 y 3 |
| `test_conservacion_energia_bordes_adiabaticos` | Todo aislado: la energía no cambia (error ~1e-12) | 4 |
| `test_fuente_uniforme_adiabatica_aumento_lineal` | $dT/dt = Q/(\rho c_p)$ en todas las celdas | 5 |
| `test_neumann_flujo_entrante_suma_energia` | Entra $q$ por una cara: $E(t) = E_0 + qAt$ exacto | 6 |
| `test_robin_bloque_caliente_se_enfria_hacia_el_aire` | El calor sale hacia el aire más frío; ninguna celda baja de $T_{aire}$ | 7 |
| `test_estacionario_robin_con_moldaje` | Estacionario con moldaje vs. resistencias en serie | 7 |
| `test_3d_con_caras_laterales_aisladas_es_igual_al_1d` | Barra 3D con costados aislados = solver 1D (diferencia < 1e-10) | conecta con la semana 1 |
| `test_interfaz_hormigon_sobre_suelo` | Dos materiales: perfil lineal por tramos con quiebre en la interfaz | interfaces |
| `test_paso_maximo_es_fourier_*` | $Fo \le 1/4$ (2D) y $1/6$ (3D) | estabilidad |
| `test_dt_inestable_lanza_error` | El solver se niega a correr con $\Delta t$ inestable | estabilidad |
| `test_simetria_del_campo` | Cubo simétrico → campo simétrico (detecta errores de índices) | extra |
| `test_cara_mal_escrita_lanza_error` | Un nombre de cara mal escrito no pasa en silencio | extra |
| `test_en_el_limite_es_el_metodo_de_schmidt[2]` / `[3]` | En el límite de estabilidad, la nueva T es el promedio de las vecinas (ver D.11) | Schmidt (ICH) |

Figuras (`python scripts/06_verificacion_3d.py`, ~5 s):

- [fig07](../resultados/figuras/fig07_seno_2d.png): solver 2D, solución exacta y error (máx ~1e-3 °C sobre una amplitud de 40 °C).
- [fig08](../resultados/figuras/fig08_convergencia_2d_3d.png): orden observado **2,00 en 2D y 1,99 en 3D**
  (tabla en `resultados/convergencia_2d_3d.csv`).
- [fig09](../resultados/figuras/fig09_bloque_hormigon_3d.png): bloque de hormigón 3×3×2 m con una
  fuente **ilustrativa** (todavía no es la hidratación). Muestra lo que el 1D no puede: el núcleo
  llega a ~73 °C mientras una esquina apenas pasa los 23 °C, porque la esquina pierde calor por
  tres caras.
- [fig10](../resultados/figuras/fig10_bordes_dominan.png): el mismo bloque con tres bordes distintos
  (ver D.10).

### D.7 Flujos y signos

Una sola convención en todo el código: **un flujo $F$ es positivo cuando el calor va hacia $+x$,
$+y$ o $+z$**.

| Dónde | Fórmula | Cómo leerla | Código (`solver3d.py`) |
|---|---|---|---|
| Entre dos celdas | $F = G\,(T_{atrás} - T_{adelante})$ | Si la de atrás está más caliente, $F > 0$: el calor avanza hacia la fría. Es el signo menos de Fourier, $q = -k\,\partial T/\partial x$ | línea 326 |
| Balance de la celda | entra por su cara de atrás $-$ sale por su cara de adelante | $F_{atrás} - F_{adelante} > 0$ → la celda gana calor | línea 337 |
| Borde | $q_{entrante} = G_b\,(T_{exterior} - T_{celda})$ | Positivo = el calor **entra** al hormigón. En Robin, si el aire está más frío, sale negativo: se enfría | `_flujo_entrante()` en `solver1d.py` |
| Borde "de inicio" ($x-$, $y-$, $z-$) | $F = +q_{entrante}$ | lo que entra por ahí avanza hacia $+$ | línea 329 |
| Borde "final" ($x+$, $y+$, $z+$) | $F = -q_{entrante}$ | lo que entra por ahí avanza hacia $-$ | línea 329 |

Comprobación física (está en los tests): un bloque caliente al aire frío siempre **pierde**
energía y ninguna celda baja de la temperatura del aire (`test_robin_bloque_caliente_se_enfria_hacia_el_aire`).
Si un signo estuviera mal, el bloque se calentaría solo.

### D.8 Fuente volumétrica vs. flujo superficial

Las dos agregan calor, pero **no son lo mismo**:

| | Fuente volumétrica $Q$ | Flujo superficial $q$ |
|---|---|---|
| Unidad | **W/m³** | **W/m²** |
| Dónde actúa | **dentro** de cada celda, en todo el volumen | solo **por una cara** del borde |
| En el balance | $+\,Q \cdot V$ en todas las celdas | $+\,q \cdot A$ solo en las celdas que tocan esa cara |
| Ejemplo en el proyecto | **calor de hidratación** del cemento | **radiación solar** absorbida por la superficie |
| En el código | argumento `fuente=` de `resolver()` (línea 340) | condición `Neumann(q)` en una cara |
| Test | `test_fuente_uniforme_adiabatica_aumento_lineal`: todo el bloque sube parejo | `test_neumann_flujo_entrante_suma_energia`: se calienta primero la cara |

Error típico: poner el sol como fuente $Q$ (lo repartiría en todo el volumen, cuando en
realidad solo calienta la superficie) o poner la hidratación como flujo en el borde.
La **convección no es una fuente**: es un intercambio por la frontera (Robin).

### D.9 Superficies adiabáticas y simetría

**Adiabático** = no cruza calor: $q = 0$, `Neumann(0)`. En el código, **toda cara que no se
nombra en `bordes` queda adiabática**, y `bordes={}` es un bloque totalmente aislado (así lo usan
los tests de energía y de fuente).

Hay dos situaciones físicas que se modelan así:

1. **Aislación perfecta**: una idealización. La aislación real es Robin con $R_{extra}$, porque
   siempre deja pasar algo (fig10).
2. **Plano de simetría**: si el problema es simétrico respecto de un plano, la temperatura es
   igual a ambos lados y el gradiente perpendicular es cero, así que no cruza calor. La fundación
   octagonal es simétrica respecto de los planos $x$ e $y$ que pasan por su centro: se puede
   simular **solo 1/4** con esos dos cortes como caras adiabáticas, con 4 veces menos celdas.
   **Ojo:** solo vale si los bordes también son simétricos. El sol de la tarde, que calienta más
   una cara que la otra, rompe la simetría.

Cómo reconocer una superficie adiabática en un resultado: las líneas de igual temperatura llegan
**perpendiculares** a ella, porque no hay flujo que la cruce. En la fig09 eso se ve en la base.

### D.10 Por qué las condiciones de borde pueden dominar el resultado

La [fig10](../resultados/figuras/fig10_bordes_dominan.png) usa el **mismo bloque, el mismo hormigón
y la misma fuente** que la fig09. Solo cambia el borde (costados y cara superior; la base sigue aislada):

| Borde | Núcleo máx. | Diferencial núcleo − superficie máx. |
|---|---|---|
| Superficie desnuda al aire ($h = 10$ W/m²K) | 73,5 °C | **40,9 °C** (no cumple 20 °C) |
| Con 5 cm de aislación ($R \approx 1{,}4$ m²K/W) | 76,7 °C | **11,2 °C** (cumple) |
| Todo aislado (adiabático) | 78,5 °C | **0 °C** (todo sube igual) |

Por qué pasa:

- **El núcleo casi no se entera de los bordes.** En 7 días el calor viaja solo
  $\sqrt{\alpha t} \approx \sqrt{9{,}3\cdot10^{-7}\cdot 6\cdot10^5} \approx 0{,}75$ m. El centro de
  un elemento grueso se comporta casi como si estuviera aislado: su máximo cambia poco (73–78 °C).
- **La superficie depende totalmente del borde.** Con aire, la superficie queda cerca del aire.
  Con aislación, queda cerca del núcleo.
- **El diferencial es la resta de ambos**, así que lo controla el borde: pasa de 41 °C a 0 °C
  sin tocar el hormigón.
- **Número de Biot** $Bi = hL/k$: compara cuánto calor se lleva el aire con cuánto puede conducir
  el hormigón. Con $h = 10$, $L = 1$ m y $k = 2$: $Bi = 5 \gg 1$. El aire saca calor más rápido
  de lo que el interior lo repone, y se forma un gradiente fuerte. Con aislación,
  $h_{efectivo} = 1/(1/h + R) \approx 0{,}67$ y $Bi \approx 0{,}3$: gradiente suave.

**Conclusión para el proyecto:** dos criterios tiran en sentidos opuestos. Aislar baja el
diferencial (bueno), pero sube un poco $T_{max}$ (malo). Por eso moldaje, aislación y retiro son
**decisiones del plan** que hay que simular, no detalles.

### D.11 Método de Schmidt (Anexo 4 del ICH)

> Lectura del Anexo 4 pendiente. Esta sección usa la formulación estándar del método.

El método de Schmidt es un método de diferencias finitas explícito, pensado originalmente para
calcular **a mano o con una planilla**. Es **nuestro mismo esquema** usado justo en el límite de
estabilidad.

Nuestra actualización en 1D (malla uniforme, sin fuente) es:

$$
T_i^{n+1} = (1 - 2Fo)\,T_i^n + Fo\,(T_{i-1}^n + T_{i+1}^n)
$$

Schmidt elige el paso de tiempo para que $Fo = 1/2$. Así el término de la propia celda se anula y queda

$$
T_i^{n+1} = \frac{T_{i-1}^n + T_{i+1}^n}{2}
$$

es decir, **la nueva temperatura es el promedio de las dos vecinas**. Lo mismo vale en más dimensiones:

| | Fo de Schmidt | Regla |
|---|---|---|
| 1D | 1/2 | promedio de 2 vecinas |
| 2D | 1/4 | promedio de 4 vecinas |
| 3D | 1/6 | promedio de 6 vecinas |

Son **exactamente los límites de estabilidad** de la D.5. Lo comprobamos con
`test_en_el_limite_es_el_metodo_de_schmidt`: con $\Delta t$ en el límite, nuestro solver da el
promedio de las vecinas.

Diferencias con nuestro solver:

| | Schmidt clásico | Nuestro solver |
|---|---|---|
| Paso de tiempo | fijo por la malla ($\Delta t = \Delta x^2/2\alpha$ en 1D) | cualquiera hasta el límite (por defecto 90 %) |
| Materiales | uno, o tratamiento especial | cada celda el suyo (resistencias en serie) |
| Bordes | reglas para cada caso | Dirichlet/Neumann/Robin generales |
| Riesgo | justo en el límite: estable, pero sin margen | 10 % de margen |

Por qué el ICH lo menciona: muestra que un método simple de celdas basta para hormigón masivo y
que se extiende a 2D/3D, secuencia de vertidas, aislación y otros materiales. Ese es el enfoque
del curso.

### D.12 Cómo usarlo

```python
from mcplan.solver3d import Material, Malla3D, Dirichlet, Neumann, Robin, resolver

hormigon = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)
malla = Malla3D.uniforme((3.0, 3.0, 2.0), (30, 30, 20), hormigon)   # celdas de 10 cm

aire = Robin(h=10.0, T_inf=15.0)
res = resolver(
    malla,
    T_inicial=18.0,
    bordes={"x-": aire, "x+": aire, "y-": aire, "y+": aire, "z+": aire},  # z- (base) aislada
    t_final=7 * 86400.0,
    fuente=1000.0,                 # W/m³; también un arreglo por celda o una función Q(t, T)
)
res.T[-1][15, 15, 0]               # T del último instante en la celda (i, j, l)
res.energia(), res.fourier, res.dt
```

### D.13 Limitaciones actuales (próximas semanas)

- La malla es una **caja**: la fundación octagonal real necesita marcar qué celdas son hormigón,
  suelo o aire (**voxelización del IFC**, semana 3).
- Todas las vertidas existen desde $t = 0$; falta **activar cada vertida** a su hora (semana 3).
- La fuente $Q$ ya puede depender de $T$ (`Q(t, T)`), pero falta el **modelo de hidratación** (semana 3).
- El borde es igual en toda una cara; el retiro de moldaje y el clima real llegan en la semana 4.
