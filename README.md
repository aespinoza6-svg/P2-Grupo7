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
- **Coordenadas meteorológicas:** latitud `-37.692100`, longitud `-72.624000` (ver [`config/sitio.yaml`](config/sitio.yaml))

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

## Estado de avance

| Semana | Entrega | Estado |
|---|---|---|
| W1 | E01 · BIM + difusión + solver 1D | ✅ IFC inicial, solver 1D, 23 tests, guía del oral |
| W1 | E02 · round-trip IFC + demo 1D | 🟡 comparación de versiones lista; falta practicar la edición en Bonsai |
| W1 | E03 · verificación 1D | 🟡 convergencia lista (fig06); falta el informe (≤ 800 palabras) |
| W2 | E04–E06 · 2D/3D + condiciones de borde | ⬜ |
| W3 | E07–E09 · voxelización + vertidas + hidratación | ⬜ |
| W4 | E10–E12 · entorno, moldaje, clima | ⬜ |
| W5 | E13–E15 · incertidumbre, Monte Carlo, HPC | ⬜ |

## Estructura del repositorio

```
P2-Grupo7/
├── README.md               <- este archivo: qué es, cómo instalar, cómo correr
├── pyproject.toml          <- define el paquete Python "mcplan" y sus dependencias
├── requirements.txt        <- versiones exactas probadas (entorno reproducible)
│
├── config/                 <- ENTRADAS: todo lo que define el caso (sin tocar código)
│   ├── sitio.yaml              coordenadas oficiales del grupo
│   ├── fundacion.yaml          geometría y material (supuestos marcados [PROVISIONAL])
│   └── plan_v01.yaml           plan de vertidas de ejemplo
│
├── src/mcplan/             <- CÓDIGO del proyecto (paquete importable)
│   ├── solver1d.py             solver térmico 1D explícito (balance por celdas)
│   ├── soluciones_analiticas.py  soluciones exactas para verificar
│   └── bim.py                  crear / leer el modelo IFC con IfcOpenShell
│
├── scripts/                <- PROGRAMAS a ejecutar, numerados en orden
│   ├── 01_crear_modelo_ifc.py
│   ├── 02_inspeccionar_ifc.py
│   ├── 03_demo_solver_1d.py
│   └── 04_convergencia_1d.py
│
├── tests/                  <- verificación automática (pytest)
├── modelos/                <- modelos IFC (se abren en Bonsai)
├── resultados/             <- SALIDAS: figuras y tablas generadas por los scripts
├── docs/                   <- explicación de teoría, código y decisiones (leer en orden)
└── entregas/               <- texto de cada entrega (repositorio + commit evaluado)
```

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

Con el entorno activado, desde la raíz del repositorio:

```bash
pytest -v                                                        # 1. verificar (23 tests)
python scripts/01_crear_modelo_ifc.py                            # 2. generar el IFC desde config/
python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc   # 3. leer el IFC sin Blender
python scripts/03_demo_solver_1d.py                              # 4. figuras del solver 1D
python scripts/04_convergencia_1d.py                             # 5. estudio de convergencia
```

Las figuras quedan en [`resultados/figuras/`](resultados/figuras/).
Para comparar un IFC editado en Bonsai con el original:

```bash
python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc --comparar modelos/fundacion_G7_plan_v01_bonsai.ifc
```

## Documentación

Ver [`docs/README.md`](docs/README.md). Para el oral E01, partir por
[`docs/04_guia_oral_E01.md`](docs/04_guia_oral_E01.md).

## Forma de trabajo

Ramas + Pull Requests revisados por otro integrante; cada entrega se marca con un tag sobre el
commit evaluado. Detalle en [`docs/06_flujo_de_trabajo_git.md`](docs/06_flujo_de_trabajo_git.md).

## Repositorio

URL: https://github.com/aespinoza6-svg/P2-Grupo7
