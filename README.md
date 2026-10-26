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
│   ├── plan_v01.yaml      plan provisional antiguo (referencia histórica)
│   └── plan_v02.yaml      geometría docente de cinco vertidas
│
├── mcplan/            EL CÓDIGO (lo que se explica en el oral)
│   ├── solver1d.py        solver térmico 1D (balance de energía por celdas)
│   ├── solver3d.py        el mismo método en 2D/3D (6 vecinas, 6 caras de borde)
│   ├── analiticas.py      soluciones exactas para verificar el solver
│   └── bim.py             crear y leer el modelo IFC
│
├── scripts/           LO QUE SE EJECUTA, en orden
│   ├── 01_crear_modelo_ifc.py
│   ├── 02_inspeccionar_ifc.py
│   ├── 03_figuras_solver_1d.py
│   ├── 04_convergencia_1d.py
│   ├── 05_jugar_barra.py       laboratorio: perillas de la ecuación 1D
│   └── 06_verificacion_3d.py   figuras de verificación 2D/3D
│
├── tests/             verificación automática (pytest)
├── fundaciones/       BIM y presentación visual
│   ├── fundacion.ifc                    IFC fuente, regenerable desde config/
│   ├── fundacion_araucania.blend        escena visual de Blender/Bonsai
│   └── vista_fundacion_araucania.png    render de esa escena
├── resultados/        (generado) figuras y tablas
├── Antecedentes/      material complementario del curso y bibliografía
│
└── docs/
    ├── 1_teoria.md        BIM/IFC + física + solver 1D + solver 2D/3D
    ├── 2_guia_oral.md     preguntas y respuestas para los orales
    ├── 3_supuestos.md     decisiones y supuestos provisionales
    └── informe_E03.md     informe de verificación 1D
```

Lo marcado como *(generado)* no se edita a mano: se regenera con los scripts.

## Antecedentes del curso

La carpeta [`Antecedentes/`](Antecedentes/) reúne material de referencia para el equipo:

- presentación inicial del Proyecto 2;
- guías G00–G08 y material complementario sobre difusión térmica;
- ejemplos didácticos de difusión 1D;
- artículo de Ju y Lei (2019), referencia experimental del caso base.

Estos archivos apoyan el estudio y la revisión oral. El código ejecutable y los resultados del
proyecto se mantienen en `mcplan/`, `scripts/`, `tests/`, `fundaciones/` y `resultados/`.

## Instalación

Requisitos: **Python 3.12** y Git. No se necesita Blender para ejecutar el código o generar el
IFC. Bonsai se usa para ver y editar el IFC; la escena `fundacion_araucania.blend` añade
materiales y ambientación visual, pero no modifica ni reemplaza `fundacion.ifc`.

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
pytest -v                                                                 # 1. verificar el proyecto
python scripts/01_crear_modelo_ifc.py                                     # 2. generar el IFC de cinco vertidas
python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc            # 3. leer el IFC sin Blender
python scripts/03_figuras_solver_1d.py                                    # 4. figuras del solver 1D
python scripts/04_convergencia_1d.py                                      # 5. estudio de convergencia
python scripts/06_verificacion_3d.py                                      # 6. verificación 2D/3D
```

Para comparar un IFC editado en Bonsai con el original, primero guárdalo como
`fundaciones/fundacion_bonsai.ifc` usando **Save IFC Project As…**. Luego ejecuta:

```bash
python scripts/02_inspeccionar_ifc.py fundaciones/fundacion.ifc --comparar fundaciones/fundacion_bonsai.ifc
```

## Entregas

| ID | Fecha | Actividad | Estado | Commit evaluado |
|---|---|---|---|---|
| E01 | mié 21 oct | Oral: BIM + difusión + solver 1D | ✅ listo | *pendiente* |
| E02 | jue 22 oct | Lab: round-trip IFC + demo 1D | 🟡 falta practicar la edición en Bonsai | *pendiente* |
| E03 | vie 23 oct | Verificación 1D → [informe](docs/informe_E03.md) | ✅ listo | *pendiente* |
| E04 | mié 28 oct | Oral: 2D/3D + condiciones de borde → [teoría, Parte D](docs/1_teoria.md#parte-d-solver-2d3d-mismo-método-en-tres-direcciones) | 🟡 en curso | |
| E05 | jue 29 oct | Lab: batería de verificación térmica (`pytest -v`, `scripts/06`) | 🟡 en curso | |
| E06 | vie 30 oct | Dossier de verificación + BIM del caso | ⬜ | |
| E07–E15 | 4 nov – 20 nov | vertidas, hidratación, clima, Monte Carlo | ⬜ | |

Texto para Canvas (cada entrega):

```
Repositorio: https://github.com/aespinoza6-svg/P2-Grupo7
Commit evaluado: https://github.com/aespinoza6-svg/P2-Grupo7/commit/<SHA>
```

<details>
<summary>Nota opcional E01 (93 palabras)</summary>

La geometría sigue la simplificación docente documentada en `Antecedentes/GUIA_BONSAI_FUNDACION.md`:
planta octagonal de 18 m entre caras, transición de 1,10 a 3,30 m y pedestal superior de 1,00 m.
El IFC de cinco vertidas se genera por script con IfcOpenShell en `fundaciones/fundacion.ifc` y
se abre/edita en Bonsai; el round-trip editando en Bonsai se muestra en E02. El plan es una geometría docente, no una
recomendación constructiva optimizada. El solver 1D
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
