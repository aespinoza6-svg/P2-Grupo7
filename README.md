# Proyecto 2 — Hormigón masivo

Métodos Computacionales 2026 · Universidad de los Andes

## Grupo 7

### Integrantes

- Antonia Espinoza Rubio
- Natalia Paz Godoy Jaramillo
- Pedro Pablo Castillo Muratto

### Sitio asignado

- **Lugar:** Renaico, Región de La Araucanía, Chile
- **Proyecto eólico de referencia:** Parque Eólico Vergara
- **Coordenadas meteorológicas:** latitud `-37.692100`, longitud `-72.624000` (ver [`config/caso_base.yaml`](config/caso_base.yaml))

## Descripción

El proyecto consiste en desarrollar un flujo reproducible para analizar y recomendar un plan de construcción por vertidas para una fundación de hormigón masivo. El análisis considera la generación de calor por hidratación, la transferencia térmica, las condiciones ambientales y la incertidumbre.

El objetivo es comparar decisiones constructivas y evaluar el riesgo térmico durante los 14 días posteriores al inicio de la primera vertida.

```
Bonsai / Blender  ->  IFC  ->  IfcOpenShell  ->  geometría + propiedades  ->  voxelización  ->  solver propio
                                                                                        ->  Monte Carlo  ->  decisión
```

## Criterios térmicos del curso

- Temperatura máxima del hormigón: **70 °C**.
- Diferencial máximo núcleo–superficie: **20 °C**.

Estos límites son criterios de control del riesgo definidos para el proyecto docente; cumplirlos no constituye una demostración de ausencia de agrietamiento.

## ¿Por dónde empiezo?

1. Este README: qué hay en cada carpeta y cómo correr todo.
2. [`docs/1_teoria.md`](docs/1_teoria.md): BIM/IFC, física y solver, explicados.
3. El código que se explica en el oral: [`mcplan/solver1d.py`](mcplan/solver1d.py) y [`mcplan/bim.py`](mcplan/bim.py).
4. [`docs/2_guia_oral.md`](docs/2_guia_oral.md): preguntas probables con sus respuestas.

## Estructura del repositorio

```
P2-Grupo7/
├── README.md          <- este archivo
├── requirements.txt   <- dependencias (versiones exactas probadas)
├── pyproject.toml     <- define el paquete "mcplan" (técnico, no se toca)
│
├── config/            ENTRADAS: datos y decisiones (se cambian aquí, no en el código)
│   ├── caso_base.yaml     sitio + geometría + hormigón
│   └── plan_v01.yaml      plan de vertidas (ejemplo)
│
├── mcplan/            EL CÓDIGO (lo que se explica en el oral)
│   ├── solver1d.py        solver térmico 1D (balance de energía por celdas)
│   ├── analiticas.py      soluciones exactas para verificar el solver
│   └── bim.py             crear y leer el modelo IFC
│
├── scripts/           LO QUE SE EJECUTA, en orden
│   ├── 01_crear_modelo_ifc.py
│   ├── 02_inspeccionar_ifc.py
│   ├── 03_figuras_solver_1d.py
│   └── 04_convergencia_1d.py
│
├── tests/             verificación automática (pytest)
├── modelos/           (generado) modelo IFC, se abre en Bonsai
├── resultados/        (generado) figuras y tablas
│
└── docs/
    ├── 1_teoria.md        BIM/IFC + física + solver
    ├── 2_guia_oral.md     preguntas y respuestas para los orales
    ├── 3_supuestos.md     decisiones y supuestos provisionales
    └── informe_E03.md     informe de verificación 1D
```

Lo marcado como *(generado)* no se edita a mano: se regenera con los scripts.

## Instalación

Requisitos: **Python 3.12** y Git. (No se necesita Blender para correr el código; Bonsai
solo se usa para ver y editar el modelo IFC.)

**Windows (PowerShell):**

```bash
git clone https://github.com/aespinoza6-svg/P2-Grupo7.git
cd P2-Grupo7
py -3.12 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

**macOS / Linux:**

```bash
git clone https://github.com/aespinoza6-svg/P2-Grupo7.git
cd P2-Grupo7
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

`requirements.txt` instala numpy, matplotlib, pyyaml, ifcopenshell, pytest y nuestro paquete
`mcplan` (en modo editable).

> **Nota Windows:** si `import ifcopenshell` falla con *"El nombre del archivo o la extensión es
> demasiado largo"*, la ruta del entorno es demasiado larga. Crear el entorno en una ruta corta
> (p. ej. `py -3.12 -m venv C:\venvs\p2g7`). Además conviene no crear `.venv` dentro de una
> carpeta sincronizada por OneDrive (son cientos de MB).

## Uso

Siempre **desde la carpeta del repositorio** y con el entorno activado:

```bash
pytest -v                                                                 # 1. verificar (23 tests)
python scripts/01_crear_modelo_ifc.py                                     # 2. generar el IFC desde config/
python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc   # 3. leer el IFC sin Blender
python scripts/03_figuras_solver_1d.py                                    # 4. figuras del solver 1D
python scripts/04_convergencia_1d.py                                      # 5. estudio de convergencia
```

Para comparar un IFC editado en Bonsai con el original:

```bash
python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc --comparar modelos/fundacion_G7_plan_v01_bonsai.ifc
```

## Entregas

| ID | Fecha | Actividad | Estado | Commit evaluado |
|---|---|---|---|---|
| E01 | mié 21 oct | Oral: BIM + difusión + solver 1D | ✅ listo | *pendiente* |
| E02 | jue 22 oct | Lab: round-trip IFC + demo 1D | 🟡 falta practicar la edición en Bonsai | *pendiente* |
| E03 | vie 23 oct | Verificación 1D → [informe](docs/informe_E03.md) | ✅ listo | *pendiente* |
| E04–E15 | 28 oct – 20 nov | 2D/3D, vertidas, hidratación, clima, Monte Carlo | ⬜ | |

Texto para Canvas (cada entrega):

```
Repositorio: https://github.com/aespinoza6-svg/P2-Grupo7
Commit evaluado: https://github.com/aespinoza6-svg/P2-Grupo7/commit/<SHA>
```

<details>
<summary>Nota opcional E01 (93 palabras)</summary>

La geometría es provisional: interpretamos 3,30 m como altura total y supusimos un pedestal de
1,00 m hasta recibir la geometría oficial (supuestos documentados en docs/3_supuestos.md). El IFC
se genera por script con IfcOpenShell y se abre/edita en Bonsai; el round-trip editando en Bonsai
se mostrará en E02. El plan de vertidas del IFC es un ejemplo, no una recomendación. El solver 1D
(Dirichlet, Neumann, Robin, fuente) está verificado con solución exacta, conservación de energía
y convergencia (orden 2 en espacio, 1 en tiempo). Faltan: calibrar las propiedades del hormigón e
incorporar hidratación.

</details>

## Forma de trabajo con Git

1. **`main` siempre funciona.** Cada tarea se hace en una rama (`e04-solver-2d`) y entra por
   **Pull Request** revisado por otro integrante: así todos conocen todo el código.
2. **Antes de cada commit:** `pytest`.
3. **Mensajes de commit:** `tipo: descripción`, con tipo `feat` (nuevo), `fix` (corrección),
   `test`, `docs`, `refactor` (reorganizar sin cambiar resultados) o `chore` (configuración).
4. **Commit evaluado:** unir a `main`, obtener el SHA con `git rev-parse HEAD`, marcarlo con
   `git tag -a E01 -m "Entrega E01"` + `git push origin E01`, y completar la tabla de entregas.

## Repositorio

URL: https://github.com/aespinoza6-svg/P2-Grupo7
