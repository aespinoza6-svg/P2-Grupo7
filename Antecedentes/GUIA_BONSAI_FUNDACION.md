# Guía paso a paso: modelar la fundación del capstone en Blender + Bonsai

> **Design by Simulation in Civil Engineering — Capstone de Hormigón Masivo**  
> Esta guía asume **cero experiencia previa en Blender**.

---

## 0. ¿Qué vas a lograr?

Al terminar esta guía deberías tener un archivo:

```text
foundation.ifc
```

que contiene cinco volúmenes BIM:

```text
Vertida_01
Vertida_02
Vertida_03
Vertida_04
Vertida_05
```

y que puede ser leído desde Python con IfcOpenShell.

El flujo que queremos es:

![Workflow](images/05_workflow.png)

**Importante:** Blender no será nuestro solver. Bonsai tampoco.

- **Blender**: interfaz 3D para crear y editar geometría.
- **Bonsai**: capa BIM/IFC dentro de Blender.
- **IFC**: archivo abierto que guarda geometría + identidad + información.
- **IfcOpenShell**: biblioteca Python que leerá el IFC.
- **nuestro código**: voxeliza y resuelve difusión térmica.

La fuente de verdad del BIM será el archivo **`.ifc`**, no el `.blend`.

---

# 1. Geometría que modelaremos

El caso del curso está inspirado en la fundación eólica de Ju & Lei (2019).

El artículo reporta:

- planta global de aproximadamente **18 m × 18 m**;
- espesor máximo de **3.30 m**;
- parte superior cilíndrica de aproximadamente **7 m de diámetro**.

La Figura 2 del paper muestra además que la fundación disminuye de espesor hacia el borde.

Sin embargo, el paper **no entrega todas las cotas necesarias para reconstruir de forma única el sólido 3D**. Por eso no intentaremos fingir una precisión que la fuente no tiene.

## 1.1 Geometría docente oficial de esta guía

Usaremos una simplificación reproducible:

- planta exterior octagonal;
- **18 m entre caras opuestas**;
- espesor en borde: **1.10 m**;
- espesor central: **3.30 m**;
- plataforma/pedestal central de **7.0 m de diámetro**;
- pedestal superior de **1.0 m de altura**, adoptado sólo para el ejercicio docente.

![Planta](images/02_planta_fundacion.png)

![Sección](images/03_seccion_fundacion.png)

> Si posteriormente el profesor entrega un IFC geométrico oficial distinto, **ese IFC reemplaza estas dimensiones**.

---

# 2. ¿Por qué cinco vertidas?

Para que el BIM represente el proceso constructivo, no queremos un único sólido.

La geometría docente se separará así:

![Vertidas](images/04_vertidas_explotadas.png)

| Vertida | Desde z [m] | Hasta z [m] | Forma |
|---|---:|---:|---|
| `Vertida_01` | 0.00 | 0.55 | prisma octagonal |
| `Vertida_02` | 0.55 | 1.10 | prisma octagonal |
| `Vertida_03` | 1.10 | 2.20 | tronco octagonal inferior |
| `Vertida_04` | 2.20 | 3.30 | tronco octagonal superior |
| `Vertida_05` | 3.30 | 4.30 | cilindro central |

Estas cinco piezas **no implican todavía que éste sea el plan constructivo óptimo**.

Son el modelo inicial sobre el cual el grupo podrá posteriormente discutir:

- unir vertidas;
- modificar intervalos;
- cambiar temperaturas de colocación;
- cambiar orden;
- cambiar aislación;
- estudiar otras estrategias.

---

# 3. Un poco de geometría para no equivocarnos

Blender pide el **radio al vértice** de un octágono, pero nosotros conocemos la distancia entre caras.

Para un octágono regular:

$$
R=\frac{a}{\cos(\pi/8)}
$$

donde $a$ es la apotema.

Para 18 m entre caras:

$$
a=9\text{ m}
$$

y:

$$
R_{\text{ext}}
=
\frac{9}{\cos(\pi/8)}
\approx
9.74153\text{ m}
$$

Para una plataforma central de 7 m entre caras del octágono auxiliar:

$$
a=3.5\text{ m}
$$

y:

$$
R_{\text{int}}
\approx
3.78837\text{ m}
$$

Para dividir la transición en dos niveles iguales usamos:

$$
R_{\text{medio}}
=
\frac{R_{\text{ext}}+R_{\text{int}}}{2}
\approx
6.76495\text{ m}
$$

Estos valores ya están en:

```text
course_geometry.json
```

---

# 4. Instalar Blender y Bonsai

## 4.1 Blender

Descarga Blender desde:

https://www.blender.org/

## 4.2 Bonsai

La instalación oficial actual es:

1. abrir Blender;
2. `Edit > Preferences`;
3. ir a `Get Extensions`;
4. habilitar acceso online si Blender lo pide;
5. buscar **Bonsai**;
6. presionar `Install`.

Referencia oficial:

https://docs.bonsaibim.org/quickstart/installation.html

Al instalarse correctamente debería aparecer en el menú:

```text
File > New IFC Project
File > Open IFC Project
```

---

# 5. Blender para alguien que nunca ha usado Blender

No necesitas aprender Blender completo.

Para este proyecto puedes sobrevivir con aproximadamente diez acciones.

## 5.1 Reconocer la pantalla

![Layout](images/01_blender_layout.png)

Las tres zonas relevantes son:

### 3D Viewport

La zona grande central.

Aquí:

- ves la fundación;
- seleccionas objetos;
- mueves la cámara;
- agregas geometría.

### Outliner

Normalmente arriba a la derecha.

Es el árbol/lista de objetos.

### Properties

Normalmente abajo a la derecha.

Bonsai añade aquí paneles como:

- Project Overview;
- Object Information;
- Geometry and Materials.

---

# 6. Navegar sin aprender atajos

La forma más segura para un principiante es usar el **Navigation Gizmo** de la esquina superior derecha del viewport.

La documentación oficial de Bonsai recomienda usarlo para:

- **Orbit**: girar la vista;
- **Zoom**;
- **Pan**;
- vistas superior/frontal/lateral;
- cambiar perspectiva/ortográfica.

Referencia:

https://docs.bonsaibim.org/quickstart/explore_model.html

## 6.1 Con mouse

Si tienes mouse de tres botones:

| Acción | Control |
|---|---|
| orbitar | botón central |
| pan | `Shift` + botón central |
| zoom | rueda |
| seleccionar | click izquierdo |

## 6.2 Atajos que sí conviene conocer

| Acción | Tecla |
|---|---|
| mover objeto | `G` |
| rotar | `R` |
| escalar | `S` |
| panel numérico de transformaciones | `N` |
| vista superior | `Numpad 7` |
| vista frontal | `Numpad 1` |
| vista lateral | `Numpad 3` |
| cancelar operación | `Esc` |
| abrir opciones de la última operación | `F9` |

Para este proyecto **preferiremos escribir dimensiones numéricas**, no mover objetos “a ojo”.

---

# 7. Object Mode versus Edit Mode

Blender tiene dos modos que pueden confundirte.

## Object Mode

Trabajas con el objeto completo.

Ejemplo:

```text
mover toda Vertida_03
```

## Edit Mode

Modificas vértices, caras y aristas del objeto.

Ejemplo:

```text
mover un vértice de Vertida_03
```

Se alterna normalmente con:

```text
Tab
```

### Regla del curso

Durante esta guía:

> **permanece en Object Mode salvo que una instrucción diga explícitamente lo contrario.**

Nuestra geometría puede construirse enteramente con primitivas paramétricas.

---

# 8. Crear el proyecto IFC

Abre una sesión nueva de Blender.

Luego:

```text
File
> New IFC Project
> New Metric (m) Project
```

Bonsai crea un proyecto IFC4 y una jerarquía espacial básica:

```text
IfcProject
└── IfcSite
    └── IfcBuilding
        └── IfcBuildingStorey
```

La documentación oficial llama a esto el **Spatial Tree**.

Referencia:

https://docs.bonsaibim.org/quickstart/create_model.html

---

# 9. Guardar inmediatamente

Haz:

```text
File > Save IFC Project
```

y guarda como:

```text
foundation.ifc
```

## Muy importante

En un workflow Bonsai nativo:

- `.ifc` = modelo BIM;
- `.blend` = configuración auxiliar de Blender.

No dependas del `.blend` como fuente de verdad.

La documentación actual de Bonsai advierte explícitamente que el IFC es el archivo crítico del proyecto.

---

# 10. Antes de modelar: cursor y coordenadas

Queremos que el centro de la fundación sea:

$$
(x,y)=(0,0)
$$

y que la base sea:

$$
z=0.
$$

Cada vez que agregues una pieza:

1. colócala inicialmente cerca del origen;
2. presiona `N`;
3. usa el panel `Item / Transform`;
4. escribe exactamente:
   - `Location X`;
   - `Location Y`;
   - `Location Z`;
   - `Rotation Z`.

No confíes en arrastrar manualmente.

---

# 11. Vertida 01

## 11.1 Agregar un cilindro de 8 lados

En el viewport:

```text
Add
> Mesh
> Cylinder
```

Inmediatamente abre las opciones con:

```text
F9
```

Configura:

```text
Vertices = 8
Radius   = 9.74153
Depth    = 0.55
```

Un cilindro con 8 vértices es un prisma octagonal.

## 11.2 Ubicación

En el panel `N > Item`:

```text
Location X = 0
Location Y = 0
Location Z = 0.275
Rotation Z = 22.5°
```

¿Por qué z = 0.275?

Porque Blender posiciona la primitiva por su centro:

$$
z_c=\frac{0+0.55}{2}=0.275.
$$

¿Por qué 22.5°?

Para que las **caras**, no los vértices, queden alineadas con los ejes X/Y.

---

# 12. Convertir la geometría en un objeto IFC

Hasta ahora tienes geometría Blender.

Aún no necesariamente tienes un elemento BIM correctamente clasificado.

Selecciona el objeto.

Ve a:

```text
Properties
> Object Information
```

Busca la sección para asignar clase IFC.

Una opción razonable para este curso es:

```text
IfcFooting
```

y luego:

```text
Assign IFC Class
```

La documentación oficial explica que cualquier objeto que deba pertenecer al IFC debe recibir una **IFC Class**.

Referencia:

https://docs.bonsaibim.org/quickstart/create_model.html

## Nota semántica

Una modelación IFC más sofisticada podría crear:

```text
IfcFooting
└── IfcBuildingElementPart
    ├── Vertida_01
    ├── ...
```

Para este capstone **no lo exigiremos inicialmente**.

Nos importa que cada volumen:

- tenga identidad IFC;
- tenga geometría propia;
- tenga nombre;
- pueda leerse desde Python.

---

# 13. Nombrar correctamente la pieza

Con el objeto seleccionado:

```text
Properties
> Object Information
> Attributes
```

cambia:

```text
Name = Vertida_01
```

No te conformes sólo con renombrar el objeto en el Outliner.

Queremos que el **atributo IFC `Name`** sea `Vertida_01`.

---

# 14. Vertida 02

Repite el mismo procedimiento.

### Geometría

```text
Add > Mesh > Cylinder
Vertices = 8
Radius   = 9.74153
Depth    = 0.55
```

### Transformación

```text
X = 0
Y = 0
Z = 0.825
Rotation Z = 22.5°
```

porque:

$$
z_c=\frac{0.55+1.10}{2}=0.825.
$$

### BIM

```text
IFC Class = IfcFooting
Name      = Vertida_02
```

---

# 15. Vertida 03: primer tronco octagonal

Ahora usamos:

```text
Add
> Mesh
> Cone
```

No queremos una punta.

Queremos un **tronco de cono octagonal**.

En `F9`:

```text
Vertices = 8
Radius 1 = 9.74153
Radius 2 = 6.76495
Depth    = 1.10
```

Ubicación:

```text
X = 0
Y = 0
Z = 1.65
Rotation Z = 22.5°
```

porque:

$$
z_c=\frac{1.10+2.20}{2}=1.65.
$$

Luego:

```text
IFC Class = IfcFooting
Name      = Vertida_03
```

---

# 16. Vertida 04

Agrega otro:

```text
Add > Mesh > Cone
```

Parámetros:

```text
Vertices = 8
Radius 1 = 6.76495
Radius 2 = 3.78837
Depth    = 1.10
```

Ubicación:

```text
X = 0
Y = 0
Z = 2.75
Rotation Z = 22.5°
```

porque:

$$
z_c=\frac{2.20+3.30}{2}=2.75.
$$

Luego:

```text
IFC Class = IfcFooting
Name      = Vertida_04
```

---

# 17. Vertida 05: pedestal central

Agrega:

```text
Add > Mesh > Cylinder
```

Parámetros:

```text
Vertices = 48
Radius   = 3.50
Depth    = 1.00
```

No necesitamos 500 lados para que parezca circular.

48 es más que suficiente para el BIM de este curso.

Ubicación:

```text
X = 0
Y = 0
Z = 3.80
Rotation Z = 0°
```

porque:

$$
z_c=\frac{3.30+4.30}{2}=3.80.
$$

Luego:

```text
IFC Class = IfcFooting
Name      = Vertida_05
```

---

# 18. Tabla completa para revisar

Antes de continuar tu modelo debe coincidir con:

| Nombre | Primitive | Vertices | Radius 1 [m] | Radius 2 [m] | Depth [m] | Z centro [m] | Rot. Z |
|---|---|---:|---:|---:|---:|---:|---:|
| `Vertida_01` | Cylinder | 8 | 9.74153 | — | 0.55 | 0.275 | 22.5° |
| `Vertida_02` | Cylinder | 8 | 9.74153 | — | 0.55 | 0.825 | 22.5° |
| `Vertida_03` | Cone | 8 | 9.74153 | 6.76495 | 1.10 | 1.650 | 22.5° |
| `Vertida_04` | Cone | 8 | 6.76495 | 3.78837 | 1.10 | 2.750 | 22.5° |
| `Vertida_05` | Cylinder | 48 | 3.50000 | — | 1.00 | 3.800 | 0° |

---

# 19. Sanity check de volumen

Los volúmenes aproximados de la geometría docente son:

| Vertida | Volumen aproximado [m³] |
|---|---:|
| V1 | 147.63 |
| V2 | 147.63 |
| V3 | 214.22 |
| V4 | 88.92 |
| V5 | 38.48 |
| **Total** | **636.89** |

No necesitas obtener exactamente esos números desde Bonsai en esta etapa.

Sirven para detectar errores grandes como:

- usar milímetros en vez de metros;
- escribir radio = 18;
- duplicar una pieza;
- ubicar una vertida dos veces.

---

# 20. Inspeccionar el Spatial Tree

Ve a:

```text
Properties
> Project Overview
> Spatial Decomposition
```

Deberías ver una estructura equivalente a:

```text
IfcProject
└── IfcSite
    └── IfcBuilding
        └── IfcBuildingStorey
            ├── Vertida_01
            ├── Vertida_02
            ├── Vertida_03
            ├── Vertida_04
            └── Vertida_05
```

Bonsai asigna normalmente el elemento nuevo al contenedor espacial activo.

Si un objeto queda fuera del contenedor esperado, revisa:

```text
Object Information
> Spatial Container
```

---

# 21. Ocultar y aislar elementos

Para entender el modelo necesitas ser capaz de inspeccionarlo.

Algunos controles útiles:

```text
H
```

oculta selección.

```text
Alt + H
```

muestra objetos ocultos.

```text
Shift + H
```

oculta los no seleccionados.

La documentación oficial advierte que los hotkeys de Blender son **contextuales**: el mouse debe estar sobre el Viewport.

Si un atajo no funciona, usa:

```text
Object
> Show/Hide
```

---

# 22. ¿Qué información BIM agregaremos?

La geometría sola no basta.

Cada vertida deberá eventualmente guardar información como:

```text
VertidaID
Inicio_h
TemperaturaColocacion_C
TipoHormigon
Moldaje
Aislacion
RetiroMoldaje_h
```

Esto se agrupa en un **Property Set**.

## Ojo con el nombre

Los Property Sets estándar de buildingSMART usan nombres como:

```text
Pset_WallCommon
```

El prefijo `Pset_` está reservado para sets estandarizados.

Para nuestro set custom usaremos:

```text
UANDES_HormigonMasivo
```

y **no**:

```text
Pset_HormigonMasivo
```

---

# 23. Metadata en el repositorio del Grupo 7

El IFC que genera este repositorio ya incluye la metadata de cada vertida. No es necesario
crear ni ejecutar un archivo aparte llamado `add_course_properties.py`.

Desde la raíz del repositorio, con el entorno virtual activado, ejecuta:

```powershell
python scripts/01_crear_modelo_ifc.py
python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc
```

El primer comando lee `config/caso_base.yaml` y `config/plan_v02.yaml`, crea cinco elementos
`IfcFooting` y guarda `fundaciones/fundacion.ifc`. Cada vertida tiene el Property Set
`UANDES_HormigonMasivo`, con `VertidaID`, `Elemento`, `Forma`, `TipoHormigon`, `Inicio_h`,
`InicioPlanificado`, `TemperaturaColocacion_C`, `Moldaje`, `Aislacion` y `RetiroMoldaje_h`.

El segundo comando muestra la jerarquía IFC, las propiedades y las cotas/volúmenes calculados
desde la geometría del archivo. El IFC también se puede abrir en Bonsai con **File > Open IFC
Project**.

## 24. Editar metadata en Bonsai

Si cambias propiedades en Bonsai, guarda una copia como:

```text
fundaciones/fundacion_bonsai.ifc
```

No sobrescribas el IFC generado. Para comparar la copia editada con el archivo base, desde la
raíz del repositorio ejecuta:

```powershell
python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc --comparar fundaciones/fundacion_bonsai.ifc
```

La salida identifica cambios de propiedades, geometría, nombres o identidad entre ambas versiones.

## 25. Qué revisar en la salida

- Esquema IFC4 y unidades métricas.
- Cinco elementos IFC `Vertida_01` a `Vertida_05`, IDs numéricos 1 a 5 y Property Set `UANDES_HormigonMasivo`.
- Caja global cercana a `18 × 18 × 4.30 m`.
- Cotas Z consecutivas: `0–0.55`, `0.55–1.10`, `1.10–2.20`, `2.20–3.30` y `3.30–4.30 m`.
- Volumen total aproximado de **636.8 m³**; puede variar unas décimas por la aproximación faceteada del cilindro.

La caja en X/Y mide 18 m porque el octágono está orientado con sus caras alineadas a los ejes.

## 26. Si editas la geometría

La geometría generada desde YAML se cambia en `config/caso_base.yaml` y `config/plan_v02.yaml`;
luego se regenera con el comando de la sección 23. Las modificaciones manuales de Bonsai deben
guardarse en `fundaciones/fundacion_bonsai.ifc` para conservar una copia reproducible del IFC de entrada.

---
# 27. Guardar correctamente

Durante el trabajo usa:

```text
File > Save IFC Project
```

No confundas esto con guardar simplemente un `.blend`.

La documentación de Bonsai es especialmente enfática en este punto.

El archivo que debe terminar en el repositorio es:

```text
geometry/foundation.ifc
```

o equivalente.

---

# 28. ¿Debemos georreferenciar el IFC?

Para la primera parte del capstone: **no**.

La geometría debe quedar alrededor del origen:

```text
(0,0,0)
```

El sitio climático del grupo se guarda por separado en:

```text
sitio / latitude / longitude
```

La voxelización térmica no necesita coordenadas UTM reales.

Esto evita trabajar en Blender con coordenadas globales muy grandes.

Bonsai sí soporta georreferenciación IFC4, pero no aporta nada esencial a esta primera etapa.

Referencia:

https://docs.bonsaibim.org/studio-unstable/guides/authoring/advanced_modeling/georeferencing.html

---

# 29. ¿Por qué no modelar moldajes todavía?

Porque estamos separando capas.

Primero:

```text
geometría del hormigón
```

Luego:

```text
identificación de superficies
```

Luego:

```text
condiciones térmicas
```

El moldaje puede representarse posteriormente de dos maneras:

### Opción A — BIM explícito

Un objeto geométrico de moldaje.

### Opción B — propiedad/condición de borde

```text
cara lateral:
    formwork = plywood_18mm
```

Para nuestro solver térmico, la segunda alternativa probablemente será suficiente.

---

# 30. ¿Por qué no hacer booleanas complicadas?

Porque el propósito de Bonsai aquí no es evaluar destreza artística con Blender.

Necesitamos una geometría:

- cerrada;
- reproducible;
- separada por vertidas;
- fácil de leer con IfcOpenShell;
- fácil de voxelizar.

Primitivas simples son una ventaja.

---

# 31. Problemas típicos de principiante

## 31.1 Mi objeto desapareció

Probablemente:

- lo ocultaste;
- estás muy lejos;
- estás mirando desde un ángulo extraño.

Prueba:

1. selecciónalo en Outliner;
2. usa `View > Frame Selected`;
3. revisa X/Y/Z en `N > Item`.

---

## 31.2 Todo mide 18 mil metros

Confundiste:

```text
18 m
```

con:

```text
18000 m
```

El proyecto se crea con:

```text
New Metric (m) Project
```

Todos los valores de esta guía están en **metros**.

---

## 31.3 El octágono se ve “en punta”

Olvidaste:

```text
Rotation Z = 22.5°
```

---

## 31.4 La pieza empieza bajo z = 0

Recuerda que Blender posiciona la primitiva en su centro.

Para un sólido de altura $h$ que parte en $z_0$:

$$
z_c=z_0+\frac{h}{2}
$$

---

## 31.5 Cambié el nombre en el Outliner pero Python no lo encuentra

Asegúrate de editar el atributo IFC:

```text
Object Information
> Attributes
> Name
```

---

## 31.6 Cerré Blender y mi modelo BIM no aparece al abrir el `.blend`

El `.blend` no es la fuente principal del modelo IFC.

Abre:

```text
File > Open IFC Project
```

y selecciona tu `.ifc`.

---

## 31.7 No encuentro exactamente el mismo botón que muestra esta guía

Bonsai está en desarrollo activo y algunos nombres/paneles pueden cambiar.

Busca el concepto:

- `Project Overview`;
- `Object Information`;
- `IFC Class`;
- `Spatial Container`;
- `Property Sets`.

Si cambió la UI, consulta la documentación oficial correspondiente a tu versión.

---

# 32. Qué entregar para demostrar que esta etapa está lista

## BIM

- `foundation.ifc`;
- cinco vertidas;
- nombres IFC correctos;
- clases IFC asignadas.

## Evidencia

Una captura 3D donde se distingan las piezas.

Una vista superior.

Una vista lateral.

## Python

Ejecutar:

```bash
python code/check_ifc.py foundation.ifc
```

y guardar la salida o incorporarla al README.

## README

Debe indicar:

```text
Bonsai version:
Blender version:
IFC schema:
unidades:
simplificaciones geométricas:
```

---

# 33. Checklist final

- [ ] Instalé Blender.
- [ ] Instalé Bonsai.
- [ ] Creé `New Metric (m) Project`.
- [ ] Guardé `foundation.ifc`.
- [ ] La fundación está centrada en X=Y=0.
- [ ] La base comienza en z=0.
- [ ] Existen cinco objetos.
- [ ] Sus nombres IFC son `Vertida_01` ... `Vertida_05`.
- [ ] Tienen IFC Class.
- [ ] El octágono exterior mide aproximadamente 18 m entre caras.
- [ ] El espesor del borde es 1.10 m.
- [ ] El espesor central llega a 3.30 m.
- [ ] El pedestal tiene diámetro 7 m.
- [ ] Puedo ocultar y mostrar vertidas.
- [ ] Puedo cerrar y volver a abrir el IFC.
- [ ] Python puede encontrar las cinco vertidas.
- [ ] Conozco la diferencia entre geometría Blender y objeto IFC.
- [ ] Sé que el `.ifc` es la fuente de verdad.

---

# 34. Siguiente paso del capstone

Después de esta guía:

```text
foundation.ifc
      ↓
IfcOpenShell
      ↓
triangulación
      ↓
inside / outside
      ↓
voxel grid
      ↓
solver térmico
```

El siguiente problema ya no es BIM.

Es:

> **¿Cómo transformamos estos cinco sólidos IFC en una grilla cartesiana conservando la identidad de cada vertida?**

---

# 35. Referencias oficiales

## Bonsai

Documentación general:

https://docs.bonsaibim.org/

Instalación:

https://docs.bonsaibim.org/quickstart/installation.html

Crear un modelo IFC:

https://docs.bonsaibim.org/quickstart/create_model.html

Explorar IFC y navegar la interfaz:

https://docs.bonsaibim.org/quickstart/explore_model.html

Geometría y representaciones:

https://docs.bonsaibim.org/studio/guides/authoring/understanding_ifc/geometry_and_representations.html

## IfcOpenShell

Documentación:

https://docs.ifcopenshell.org/

Procesamiento geométrico:

https://docs.ifcopenshell.org/ifcopenshell-python/geometry_processing.html

API de Property Sets:

https://docs.ifcopenshell.org/autoapi/ifcopenshell/api/pset/index.html

## IFC

buildingSMART:

https://technical.buildingsmart.org/standards/ifc/

---

# 36. Videos para partir

## Bonsai BIM in Blender in 20 Minutes

Tutorial general para principiante:

https://www.youtube.com/watch?v=tyYxw5m-ovs

Incluye:

- instalación;
- abrir IFC;
- navegación;
- leer información BIM;
- modificar objetos.

## Bonsai Tutorial — Beginner Project — Part 1

Ifc Architect:

https://www.youtube.com/watch?v=Jc896Sob2bU

Es más largo de lo que necesitamos, pero es útil para ver cómo se crean y clasifican objetos en Bonsai.

---

# 37. Nota sobre la fuente geométrica

La fundación del capstone se inspira en:

Ju, Y. & Lei, H. (2019),  
**Actual Temperature Evolution of Thick Raft Concrete Foundations and Cracking Risk Analysis**,  
*Advances in Materials Science and Engineering*, Article ID 7029671.

El paper constituye la referencia física y geométrica principal, pero la forma usada en esta guía es una **simplificación docente explícita** porque la publicación no entrega una definición CAD completa e inequívoca del sólido.

