# E07: recorrido de demostración y preguntas

## Qué estamos demostrando

El mínimo E07 es conectar el modelo BIM con una grilla identificable y formular las vertidas y la hidratación. Además se implementó un primer prototipo de difusión en el dominio que crece.

La grilla, el diseño de hidratación, la activación y el solver se pueden explicar por separado. Los resultados térmicos son demostrativos: aún faltan calibración, meteorología, suelo y protecciones.

## Recorrido de unos cinco minutos

1. Abrir `fundaciones/fundacion.ifc` en Bonsai. Mostrar las cinco vertidas, sus nombres, propiedades y horarios.
2. Abrir `mcplan/voxelizacion.py`: `extraer_mallas()` lee vértices y triángulos IFC en coordenadas globales. `voxelizar_mallas()` asigna el ID a cada centro de celda.
3. Mostrar `resultados/e07_voxelizacion/voxelizacion_cortes.png` y `volumenes_por_vertida.csv`. Explicar el exterior y los escalones.
4. Abrir `mcplan/construccion.py`: identidad geométrica, máscara temporal y temperatura son datos distintos. Mostrar la edad `t_global - inicio_s`.
5. Abrir `mcplan/hidratacion.py` y `config/hidratacion_e07.yaml`: la fuente se obtiene derivando una curva adiabática; los parámetros son pedagógicos y tienen una fuente documentada.
6. Mostrar `resultados/e07_construccion/activacion_vertidas_cortes.png` y explicar una colocación.
7. Si solicitan el acoplamiento: `mcplan/solver_vertidas.py`, `SistemaActivo.actualizar()`, `potencias_W()` y el avance en `resolver_vertidas()`. Mostrar un corte de temperatura y el balance energético.

No es necesario memorizar cada línea; sí localizarla y explicar sus entradas, unidades y efecto.

## 1. ¿Qué es voxelizar?

Representar la geometría mediante una grilla de pequeños volúmenes. Para cada celda decidimos si su centro está dentro del hormigón y a qué vertida pertenece.

El criterio de centro completo es una aproximación: algunas celdas cortadas por la superficie se consideran enteras y otras quedan fuera. Por eso hay escalones y error de volumen.

## 2. ¿Cómo saben si un punto está dentro?

Leemos la superficie triangulada del IFC y comprobamos que sea cerrada. Sumamos los ángulos sólidos orientados de sus triángulos vistos desde el punto. En el interior el total es aproximadamente ±4π; fuera, cero. La caja envolvente solo permite evitar comprobaciones innecesarias.

Esto no supone que el sólido sea una caja o convexo. Hay una prueba con una cavidad que debe quedar vacía.

## 3. ¿Qué significa vertida_id?

Es la etiqueta geométrica de la celda. Los valores 1 a 5 identifican las vertidas del IFC y 0 indica que está fuera del hormigón. El 0 todavía no distingue aire de suelo.

Una celda de la vertida 3 mantiene ID=3 antes y después de colocarla. Lo que cambia con el tiempo es su máscara `activa`.

## 4. ¿Por qué celdas de 20 × 20 × 5 cm?

Buscamos una primera grilla manejable. Los 5 cm en altura permiten que las juntas a 0,55; 1,10; 2,20; 3,30 y 4,30 m coincidan con caras de celda. No es obligatorio que las celdas sean cúbicas: el balance y la estabilidad usan Δx, Δy y Δz.

## 5. ¿Cómo estiman el volumen?

Contamos las celdas de cada vertida y multiplicamos por Δx Δy Δz. Comparamos ese valor con el volumen de los triángulos IFC.

En la malla inicial: IFC ≈636,775 m³ y vóxeles ≈635,944 m³, diferencia −0,131 %. Eso verifica una aproximación geométrica; no demuestra precisión térmica ni ausencia de grietas.

## 6. ¿Qué pasa si una celda cae justo en una junta?

Se aplica una regla determinista: el interior estricto tiene prioridad sobre una superficie; si el centro solo está sobre superficies compartidas, gana el ID menor. No se cuenta dos veces. Un centro estrictamente interior a dos vertidas causa un error de traslape.

## 7. ¿Qué hace la activación temporal?

El hormigón aparece en su horario y empieza con su temperatura de colocación. Las temperaturas de las vertidas anteriores se conservan. Antes de existir no aporta capacidad térmica ni absorbe calor como un bloque frío.

En la hora 12, las edades de las tres primeras vertidas son 12, 6 y 0 horas. La cuarta y la quinta no generan calor todavía.

## 8. ¿Qué hacen si una colocación ocurre dentro de un paso de tiempo?

Dividimos el paso exactamente en ese horario. Si el paso propuesto va de la hora 5 a la 7, hacemos 5→6, colocamos la vertida y luego hacemos 6→7.

## 9. ¿Cómo diseñaron la hidratación?

Partimos del ejemplo de elevación adiabática de la guía:

`theta_ad(a) = theta_inf (1-exp(-a/tau))`

Derivamos para obtener la potencia por volumen:

`q_h(a) = rho cp theta_inf/tau exp(-a/tau)` [W/m³].

La edad a empieza al colocar cada vertida. Elevación θ en °C y fuente q en W/m³ son magnitudes distintas. Multiplicar q por el volumen da potencia en W.

## 10. ¿De dónde salen 45 °C y 30 h?

Del ejemplo pedagógico de la sección 24 de la guía de difusión del curso. Están declarados como no calibrados. No proceden de un ensayo de nuestra mezcla ni del pulso de coseno del ejercicio independiente.

El modelo simple libera más potencia al comienzo y luego decae. Su forma tampoco es una representación calibrada de inducción y pico real. La dependencia con la temperatura se incorporaría mediante una edad equivalente o una cinética justificada.

## 11. ¿Qué cambia térmicamente cuando colocan una vertida encima de otra?

La superficie superior antigua deja de estar expuesta al ambiente en la zona cubierta y pasa a ser una interfaz de conducción con el nuevo hormigón. Las zonas no cubiertas siguen expuestas.

El prototipo usa contacto perfecto y conductancias de resistencias en serie. Aún falta representar moldajes, retiros y detalles de juntas si resultan relevantes.

## 12. ¿Cómo verifican la energía si aumenta el dominio?

La energía almacenada se compara con tres aportes: la que entra con el hormigón nuevo, la generada por hidratación y la que entra o sale por los bordes. La conducción interna se cancela entre celdas vecinas.

No se exige que la energía almacenada permanezca constante cuando se incorpora material o existe una fuente.

## 13. ¿Cómo eligen el paso temporal?

En cada configuración activa sumamos las conductancias S_i y calculamos C_i=ρ_i cp_i V_i. Usamos el criterio suficiente `dt <= min(C_i/S_i)`. Se actualiza al colocar hormigón porque cambian las interfaces y superficies expuestas.

La fuente depende de la edad, no de la temperatura, y su energía se integra analíticamente en cada paso. Si se añade una cinética térmicamente acoplada habrá que revisar también esa integración y su estabilidad.

## 14. ¿Un máximo menor a 70 °C significa que el plan ya está aprobado?

No. El ejemplo tiene hidratación no calibrada, aire constante y base aislada. Todavía no evalúa el diferencial núcleo–superficie a 25 mm, clima real, incertidumbre ni confiabilidad. Es una demostración del funcionamiento numérico.

Además, el máximo de las primeras 72 h no tiene por qué ser el máximo del horizonte de 14 días.

## Comandos para la demostración

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/07_voxelizar_ifc.py --comparar
.\.venv\Scripts\python.exe -X utf8 scripts/08_vertidas_hidratacion.py
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py
.\.venv\Scripts\python.exe -X utf8 -m pytest tests/test_voxelizacion.py tests/test_construccion_hidratacion.py tests/test_solver_vertidas.py -v
```

Para el oral pueden mostrar resultados ya guardados y ejecutar una prueba pequeña. La capacidad de explicarlos individualmente sigue siendo central.

## Lab: cambios en vivo

Antes del lab, activen su entorno y corran una vez el caso base, así queda listo y el ayudante no espera:

```powershell
python -X utf8 scripts/11_lab_vertidas.py
```

| Si el ayudante pide | Comando | Qué mostrar |
|---|---|---|
| Cambiar Δx y reportar el error de volumen | `python -X utf8 scripts/11_lab_vertidas.py --solo-voxelizar --dx .1 --dy .1` | La línea `Error de volumen` |
| Ocultar una vertida | `python -X utf8 scripts/11_lab_vertidas.py --solo-voxelizar --ocultar V05` | `vista_3d_sin_V05.png` |
| Mostrarla de nuevo | `python -X utf8 scripts/11_lab_vertidas.py --solo-voxelizar` | `vista_3d_todas.png` |
| Cambiar hora de inicio | `python -X utf8 scripts/11_lab_vertidas.py --inicio V05=36` | La tabla del programa y `comparacion_72h.png` |
| Cambiar temperatura de colocación | `python -X utf8 scripts/11_lab_vertidas.py --temp V05=10` | Igual |
| Probar un IFC editado en Bonsai | `python -X utf8 scripts/11_lab_vertidas.py --ifc-cambiado fundaciones/fundacion_lab.ifc` | Igual |

Los resultados quedan en `resultados/lab/`, en una carpeta por malla y otra por cambio. Git ignora esa carpeta.

**¿Cómo demuestran que la secuencia está en los inputs?** El script copia el IFC, cambia la propiedad en la copia y vuelve a correr los mismos scripts 07 y 09. El solver no se edita. El horario entra por `mcplan/construccion.py`, que lee `Inicio_h` y `TemperaturaColocacion_C` de cada vertida.

**¿Ocultar una vertida cambia el cálculo?** No, solo cambia el dibujo. Para que una vertida no entre al cálculo, pongan su inicio después del horizonte, por ejemplo `--inicio V05=100`.

**¿Por qué Tmax casi no cambia con Δx?** A las 72 h el máximo está en el núcleo, lejos de la superficie escalonada. Refinar sí cambia el volumen y el área expuesta. Tampoco esperen que el error de volumen baje siempre al refinar: los errores de distintas celdas se compensan.

**¿Qué pasa si piden un orden imposible?** El script avisa con `AVISO`, pero igual simula. El solver no revisa el orden constructivo: esa revisión es del plan, y una vertida sin apoyo se debe descartar.
