# 01 — BIM, IFC, Blender, Bonsai e IfcOpenShell

> Código: [`src/mcplan/bim.py`](../src/mcplan/bim.py) ·
> Modelo: [`modelos/fundacion_G7_plan_v01.ifc`](../modelos/fundacion_G7_plan_v01.ifc) ·
> Tests: [`tests/test_bim.py`](../tests/test_bim.py)

## 1. Conceptos

| Concepto | Qué es | Rol en nuestro proyecto |
|---|---|---|
| **BIM** (Building Information Modeling) | Una forma de trabajar con un **modelo digital de información** del proyecto: cada objeto es un elemento constructivo con **identidad, tipo, relaciones y propiedades**, no solo líneas o geometría. | Describe **qué se construye y cuándo**: la fundación, sus vertidas, el material, inicio planificado, temperatura de colocación, moldaje, aislación, retiro. **No** describe el campo de temperaturas (eso es resultado del solver). |
| **IFC** (Industry Foundation Classes) | **Formato abierto** e internacional (ISO 16739, buildingSMART) para intercambiar modelos BIM. Es un archivo de texto `.ifc`. | Es el **contrato** entre el modelo y el cálculo: Bonsai lo escribe, nuestro Python lo lee. |
| **Blender** | Programa de modelado 3D de código abierto (mallas, render). Por sí solo no sabe qué es una zapata ni qué es IFC. | Es la plataforma gráfica. |
| **Bonsai** (ex BlenderBIM) | Complemento de Blender que lo convierte en una herramienta **BIM nativa IFC**: lo que se edita *es* el IFC. Agrega clases IFC, árbol espacial, property sets, materiales, planificación. | Herramienta para **ver y editar** el modelo (geometría y propiedades de cada vertida). |
| **IfcOpenShell** | Librería de código abierto (C++ con interfaz Python) para leer, escribir y procesar geometría IFC. Bonsai está construido sobre ella. | La usamos desde Python para **crear** el IFC desde `config/` y **leerlo** (propiedades, volúmenes, cotas). |

## 2. Arquitectura y por qué el solver no depende de Blender

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

## 3. Cómo es un archivo IFC por dentro

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

## 4. Nuestro modelo

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

### 4.1 De dónde sale el modelo

El IFC **no se dibuja a mano**: lo genera `scripts/01_crear_modelo_ifc.py` desde tres archivos:

| Archivo | Contiene |
|---|---|
| `config/sitio.yaml` | coordenadas oficiales del curso para G7 |
| `config/fundacion.yaml` | geometría y propiedades del hormigón (con supuestos marcados `[PROVISIONAL]`) |
| `config/plan_v01.yaml` | las vertidas: cotas, inicio, temperatura, moldaje… (**ejemplo, no recomendación**) |

Ventajas: cambiar una decisión (p. ej. la fecha de inicio o el número de capas) es cambiar un
YAML y regenerar, sin tocar código. `validar_plan()` rechaza planes con huecos o traslapes.

### 4.2 Identidades estables

IfcOpenShell inventa un `GlobalId` aleatorio para cada objeto nuevo. Si lo dejáramos así, cada
vez que se regenera el IFC la "misma" vertida tendría otra identidad. Por eso `_guid_estable()`
calcula el `GlobalId` a partir de una clave fija (p. ej. `"vertida/V01"`) con `uuid5`:
**regenerar produce exactamente el mismo archivo** (`test_globalid_estable_al_regenerar`).

## 5. Abrir y editar el modelo en Bonsai

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

## 6. Pendiente

- Reemplazar los supuestos `[PROVISIONAL]` por la geometría oficial del profesor
  ([05_decisiones_y_supuestos.md](05_decisiones_y_supuestos.md)).
- Practicar en Bonsai (cada integrante): abrir, editar una propiedad, guardar, comparar.
- Semana 3: voxelizar la geometría del IFC (`geometria_vertidas()` ya entrega volúmenes y cotas).
