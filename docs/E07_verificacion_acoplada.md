# E07: resultados de verificación del prototipo acoplado

Corridas del 3 de noviembre de 2026. Son resultados numéricos con parámetros pedagógicos; no una predicción validada para Renaico.

## Configuración común

Grilla de 20 × 20 × 5 cm, 317 972 celdas de hormigón, cinco vertidas colocadas a 0, 6, 12, 18 y 24 horas a 18 °C. Propiedades: k=2 W/(m K), rho=2400 kg/m³ y cp=900 J/(kg K). Hidratación exponencial por edad: elevación final de 45 °C y tau=30 h, procedentes del ejemplo docente. Aire constante a 15 °C, h=10 W/(m² K), base adiabática y contacto perfecto.

## Horizonte de 14 días

Con paso máximo de 600 s y 2016 pasos, el máximo global fue **62,1640 °C a las 177,1667 h**, aproximadamente el día **7,38** desde la primera vertida. El residuo energético relativo máximo fue 1,50 × 10⁻¹⁵. El máximo de 72 h (57,3016 °C) todavía se encontraba en la etapa de calentamiento.

El residuo pequeño comprueba el balance de la implementación. No comprueba que los parámetros físicos representen una mezcla real.

## Sensibilidad temporal durante 72 horas

| Medida | Resultado |
|---|---:|
| Máxima con paso máximo de 600 s | 57,301632 °C |
| Máxima con paso máximo de 300 s | 57,297846 °C |
| Diferencia entre máximas | 0,003786 °C |
| Mayor diferencia absoluta entre celdas al finalizar | 0,010790 °C |
| Diferencia media absoluta del campo final | 0,004661 °C |
| Mayor diferencia entre máximas en las salidas horarias | 0,005803 °C |

Se mantienen geometría, parámetros y condiciones de borde. Esta comparación mide sensibilidad temporal a grilla fija durante 72 h; no estima el error espacial, no establece por sí sola un orden de convergencia y no verifica el paso de 300 s durante los 14 días completos.

## Reproducir

Desde la raíz del repositorio, con la voxelización ya generada:

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py --horizonte-h 336 --max-dt-s 600 --salida resultados/e07_simulacion_14dias
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py --horizonte-h 72 --max-dt-s 300 --salida resultados/e07_sensibilidad_dt300
.\.venv\Scripts\python.exe -X utf8 scripts/10_comparar_pasos_e07.py
```

Las carpetas indicadas se reutilizan; no se crea una carpeta por fecha de ejecución. Los resúmenes guardan los parámetros y la huella de la grilla. Si se modifica el IFC, hay que regenerar la voxelización y repetir las corridas que se deseen comparar.

## Qué falta para una conclusión térmica del proyecto

Calibrar propiedades e hidratación, incorporar clima del sitio, suelo y protecciones, estudiar la precisión espacial y evaluar el diferencial núcleo–superficie con el muestreo requerido. Que este ejemplo quede bajo 70 °C no demuestra cumplimiento completo del proyecto.

Para preparar la explicación individual: [guía oral E07](E07_guia_oral.md).
