# Decisiones y supuestos

> Regla del enunciado: *"los grupos no deben inventar silenciosamente dimensiones faltantes"*.
> Todo lo que no es un dato del curso queda registrado aquí, con su estado.
> Al cambiar algo: actualizar la fila, no borrarla (trazabilidad).

## Supuestos provisionales (pendientes de validar)

| ID | Supuesto | Dónde está | Por qué | Estado |
|---|---|---|---|---|
| S01 | Planta exterior octagonal de 18 m entre caras opuestas | `config/caso_base.yaml`, `mcplan/bim.py` | Simplificación reproducible definida en la guía docente | Aplicada |
| S02 | Espesor del borde 1,10 m; transición lineal hasta 3,30 m en el centro | `config/caso_base.yaml`, `mcplan/bim.py` | La publicación no entrega un CAD inequívoco; se usa la sección de la guía docente | Aplicada; transición octagonal faceteada |
| S03 | Plataforma octagonal central de 7 m entre caras y pedestal cilíndrico Ø7 m × 1 m; sin pernos, armaduras ni ductos | `config/caso_base.yaml`, `mcplan/bim.py` | Simplificación docente; la altura total resultante es 4,30 m | Aplicada; geometría del caso base (S01–S03) |
| S04 | Hormigón: $k = 2{,}0$ W/(m·K), $\rho = 2400$ kg/m³, $c_p = 900$ J/(kg·K) | `config/caso_base.yaml` | Valores típicos (los del código de ejemplo del curso) | Calibrar en semana 3 (Ju y Lei, ICH) |
| S05 | Plan docente de cinco vertidas, inicios a 0/6/12/18/24 h, 18 °C, hormigón M1, moldaje/aislación TBD y retiro a 72 h | `config/plan_v02.yaml` | Valores de ejemplo de la guía; **no son recomendación constructiva** | Aplicado; reemplazar al estudiar planes térmicos |
| S06 | Fecha de inicio = 2026-12-01 08:00 (−03:00) | `config/plan_v02.yaml` | Fecha demostrativa; la fecha real se confirmará más adelante | Provisional |
| S07 | Figuras 03, 04 y 09–10: base adiabática, clima ficticio o aire constante y fuente ilustrativa $Q = Q_0 e^{-t/\tau}$ | `scripts/03_figuras_solver_1d.py`, `scripts/06_verificacion_3d.py` | Ilustrar la física antes de tener clima real e hidratación | Reemplazar en semanas 3–4 |

## Decisiones de modelación y de software

| ID | Decisión | Alternativas consideradas | Justificación |
|---|---|---|---|
| D01 | Paquete Python en la carpeta `mcplan/`, instalable con `pip install -e .` | Scripts sueltos | Scripts, tests y (luego) el cluster importan el mismo código, sin trucos de rutas |
| D02 | Volúmenes finitos **centrados en celdas** | Diferencias finitas con nodos en los bordes (código de ejemplo del curso) | En 3D cada vóxel es una celda; conserva energía exactamente; interfaces y moldajes como resistencias en serie |
| D03 | Euler **explícito** | Implícito (Crank–Nicolson, Euler implícito) | Simple y explicable; con $\Delta x$ de 5–10 cm el paso estable (15–30 min) es razonable. Reevaluar si el 3D resulta lento |
| D04 | Conductancia entre celdas con **resistencias en serie** (media armónica de $k$) | Promedio aritmético de $k$ | Garantiza continuidad del flujo entre materiales distintos |
| D05 | Moldaje y aislación como $R_{extra}$ en la condición Robin | Celdas explícitas de moldaje | Capas delgadas con poca capacidad térmica: una resistencia basta. Revisar si se necesita su inercia |
| D06 | `dt` por defecto = 90 % del máximo estable; el solver **rechaza** un `dt` inestable | Dejar que el usuario elija sin control | Evita resultados basura silenciosos |
| D07 | Límite de estabilidad por **coeficientes positivos** ($\Delta t \le C_i/S_i$) | Fórmula fija $Fo \le 1/2$ | Funciona con bordes, materiales distintos y en 3D sin cambiar el código |
| D08 | IFC generado por script desde YAML (IFC4, metros), editable en Bonsai | Modelar solo a mano en Bonsai | Reproducible y parametrizable; Bonsai sigue siendo el editor |
| D09 | Cada vertida = `IfcFooting` (PAD_FOOTING) con pset propio `UANDES_HormigonMasivo` | Una `IfcFooting` dividida en partes (`IfcBuildingElementPart`) | Sigue el nombre de la guía Bonsai y es sencillo de leer y editar |
| D10 | Propiedades del material en psets estándar (`Pset_MaterialThermal`, `Pset_MaterialCommon`) | Pset propio | Cualquier software BIM las entiende |
| D11 | `GlobalId` deterministas (uuid5 de una clave) | Aleatorios (por defecto) | Regenerar el IFC no cambia identidades; comparación entre versiones |
| D12 | Octágonos y transiciones como sólidos IFC faceteados | Perfiles IFC extruidos | Representa las caras inclinadas y permite leer/triangular con IfcOpenShell |

## Preguntas abiertas para el profesor / ayudante

1. ~~¿Cuándo estará la geometría oficial del curso?~~ Resuelta: el grupo adopta la geometría octagonal de S01–S03 (la del modelo en Bonsai) como geometría del caso base.
2. ¿Es aceptable representar cada vertida como una `IfcFooting` independiente (D09)?
3. ¿El modelo IFC puede generarse por script y editarse en Bonsai, o se espera que la geometría base se modele a mano en Bonsai?
