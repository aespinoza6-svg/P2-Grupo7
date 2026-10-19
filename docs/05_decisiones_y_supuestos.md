# 05 — Registro de decisiones y supuestos

> Regla del enunciado: *"los grupos no deben inventar silenciosamente dimensiones faltantes"*.
> Todo lo que no es un dato del curso queda registrado aquí, con su estado.
> Al cambiar algo: actualizar la fila, no borrarla (trazabilidad).

## Supuestos provisionales (pendientes de validar)

| ID | Supuesto | Dónde está | Por qué | Estado |
|---|---|---|---|---|
| S01 | El "espesor 3,30 m" es la **altura total** (losa + pedestal) | `config/fundacion.yaml` | El artículo no lo define inequívocamente | Pendiente: geometría oficial del profesor |
| S02 | Altura del pedestal = **1,00 m** | `config/fundacion.yaml` | Idem | Pendiente |
| S03 | Losa = prisma recto 18 × 18 m sin taludes; pedestal = cilindro recto centrado; sin pernos, armaduras ni ductos | `config/fundacion.yaml`, `bim.py` | Simplificación para la semana 1 | Pendiente |
| S04 | Hormigón: $k = 2{,}0$ W/(m·K), $\rho = 2400$ kg/m³, $c_p = 900$ J/(kg·K) | `config/fundacion.yaml` | Valores típicos (los del código de ejemplo del curso) | Calibrar en semana 3 (Ju y Lei, ICH) |
| S05 | Plan v01: losa en 2 capas + pedestal, inicios a 0/72/144 h, 20 °C, moldaje de madera, retiro a 72 h | `config/plan_v01.yaml` | Solo para tener un IFC con vertidas; **no es la recomendación** | Se reemplaza por planes simulados (semanas 3–5) |
| S06 | Fecha de inicio = 2026-12-01 08:00 (−03:00) | `config/plan_v01.yaml` | La real es la fecha de la defensa final (TBD) | Cambiar cuando se conozca |
| S07 | Figuras 03 y 04: base adiabática, clima ficticio (onda de 24 h) y fuente ilustrativa $Q = Q_0 e^{-t/\tau}$ | `scripts/03_demo_solver_1d.py` | Ilustrar la física antes de tener clima real e hidratación | Reemplazar en semanas 3–4 |

## Decisiones de modelación y de software

| ID | Decisión | Alternativas consideradas | Justificación |
|---|---|---|---|
| D01 | Paquete Python `mcplan` en `src/`, instalable con `pip install -e .` | Scripts sueltos | Scripts, tests y (luego) el cluster importan el mismo código, sin trucos de rutas |
| D02 | Volúmenes finitos **centrados en celdas** | Diferencias finitas con nodos en los bordes (código de ejemplo del curso) | En 3D cada vóxel es una celda; conserva energía exactamente; interfaces y moldajes como resistencias en serie |
| D03 | Euler **explícito** | Implícito (Crank–Nicolson, Euler implícito) | Simple y explicable; con $\Delta x$ de 5–10 cm el paso estable (15–30 min) es razonable. Reevaluar si el 3D resulta lento |
| D04 | Conductancia entre celdas con **resistencias en serie** (media armónica de $k$) | Promedio aritmético de $k$ | Garantiza continuidad del flujo entre materiales distintos |
| D05 | Moldaje y aislación como $R_{extra}$ en la condición Robin | Celdas explícitas de moldaje | Capas delgadas con poca capacidad térmica: una resistencia basta. Revisar si se necesita su inercia |
| D06 | `dt` por defecto = 90 % del máximo estable; el solver **rechaza** un `dt` inestable | Dejar que el usuario elija sin control | Evita resultados basura silenciosos |
| D07 | Límite de estabilidad por **coeficientes positivos** ($\Delta t \le C_i/S_i$) | Fórmula fija $Fo \le 1/2$ | Funciona con bordes, materiales distintos y en 3D sin cambiar el código |
| D08 | IFC generado por script desde YAML (IFC4, metros), editable en Bonsai | Modelar solo a mano en Bonsai | Reproducible y parametrizable; Bonsai sigue siendo el editor |
| D09 | Cada vertida = `IfcFooting` (PAD_FOOTING) con pset propio `MC_Vertida` | Una `IfcFooting` dividida en partes (`IfcBuildingElementPart`) | Más simple de leer, editar en Bonsai y explicar. Reevaluar si el profesor pide otra estructura |
| D10 | Propiedades del material en psets estándar (`Pset_MaterialThermal`, `Pset_MaterialCommon`) | Pset propio | Cualquier software BIM las entiende |
| D11 | `GlobalId` deterministas (uuid5 de una clave) | Aleatorios (por defecto) | Regenerar el IFC no cambia identidades; comparación entre versiones |

## Preguntas abiertas para el profesor / ayudante

1. ¿Cuándo estará la geometría oficial del curso? (resuelve S01–S03)
2. ¿Es aceptable representar cada vertida como una `IfcFooting` independiente (D09)?
3. ¿El modelo IFC puede generarse por script y editarse en Bonsai, o se espera que la geometría base se modele a mano en Bonsai?
