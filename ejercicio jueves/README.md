# Ejercicio del jueves: difusión térmica 3D

Bloque de hormigón de **3 × 3 × 0,60 m**, inicialmente a **20 °C**, con una fuente volumétrica uniforme de calor:

`Q(t) = Qmax/2 · [1 − cos(2πt/P)]`, para `0 ≤ t < P`; cero después.

Se usa `Qmax = 1300 W/m³`, `P = 2 días`, `k = 2 W/(m·K)`, `ρ = 2400 kg/m³` y `cp = 900 J/(kg·K)`.

## Visor HTML

Descargar `visor_beta_3D.html` y abrirlo en un navegador. Funciona localmente, sin Python ni conexión. Permite cambiar β, reproducir la evolución, comparar los cuatro casos y girar el relieve térmico 3D. No guarda resultados automáticamente. Los botones de descarga son opcionales.

La altura del relieve representa la **temperatura en el plano horizontal z = 0,30 m**; no es una deformación mecánica del hormigón. Los cortes verticales corresponden al plano y = 1,50 m.

## Condiciones de borde

Cada caso aplica su condición a las seis caras:

1. Dirichlet a 20 °C.
2. Dirichlet a 30 °C.
3. Adiabático: `∂T/∂n = 0`.
4. Robin: `∂T/∂n = β(T − 20 °C)`, con normal exterior. β tiene unidades de 1/m. β = 0 coincide con adiabático; β negativo representa convección hacia aire a 20 °C, con `h = −kβ`.

## Resultados de referencia, β = 0

| Caso | Máxima central [°C] | Centro al día 2 [°C] |
|---|---:|---:|
| Dirichlet a 20 °C | 42,22 | 30,12 |
| Dirichlet a 30 °C | 51,53 | 39,95 |
| Adiabático | 72,00 | 72,00 |
| Robin, β = 0 | 72,00 | 72,00 |

Centro geométrico: (1,50; 1,50; 0,30) m. Malla de 30 × 30 × 6 celdas de 10 cm. Paso temporal de 120 s. La fuente se integra exactamente en cada paso; la conducción se avanza mediante Euler explícito. Las temperaturas se representan en centros de celda y se interpolan al centro geométrico.

Se verificó el caso adiabático frente a su solución exacta y los casos Dirichlet con una malla de 5 cm. El motor numérico del HTML se contrastó con el solver Python del repositorio.

## Ejecutar el código Python

Desde la raíz del repositorio, después de instalar `requirements.txt`:

```powershell
.\.venv\Scripts\python.exe -X utf8 "ejercicio jueves/ejercicio_jueves_3D.py"
```

Para una ejecución rápida sin regenerar películas:

```powershell
.\.venv\Scripts\python.exe -X utf8 "ejercicio jueves/ejercicio_jueves_3D.py" --sin-movies
```

El código Python crea una subcarpeta de resultados por ejecución. El visor HTML calcula en memoria y no crea carpetas.

## Evidencia incluida

- `movie_cortes_comparacion.gif`: cuatro cortes verticales sincronizados.
- `movie_*.gif`: películas individuales.
- `comparacion_condiciones.png`, `fuente_temperaturas_energia.png`, `cortes_verticales_dia2.png`: gráficos.
- `*_historia.csv`: temperaturas, fuente y balance de energía.
- `*_campos_3D.npz`: campos 3D, tiempos y coordenadas.
- `resumen_casos.csv` y `resumen.json`: resultados resumidos.

Este ejercicio es independiente del modelo IFC de la fundación.
