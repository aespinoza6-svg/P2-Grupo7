# E07: activación temporal e hidratación

Este avance completa el **diseño inicial** de hidratación y la **estructura** de activación temporal solicitados por E07. Usa la grilla real obtenida del IFC. Las referencias adiabáticas de este documento no son simulaciones de difusión; para el prototipo acoplado posterior, ver [E07_simulacion_acoplada.md](E07_simulacion_acoplada.md).

## Tres conceptos separados

1. **Identidad geométrica:** `vertida_id` indica a qué vertida pertenece una celda. No cambia con el tiempo.
2. **Estado constructivo:** `activa` indica si ese hormigón ya fue colocado. Exterior y vertidas futuras no participan como hormigón activo.
3. **Edad del hormigón:** `a = t_global - t_colocación`. Cada vertida tiene su propio reloj de hidratación.

Por ejemplo, en la hora global 12, la vertida 1 tiene 12 h, la vertida 2 tiene 6 h y la vertida 3 acaba de colocarse. Las vertidas 4 y 5 aún no existen como hormigón colocado y no generan calor.

## Qué se lee del IFC

La grilla guarda `Inicio_h`, `TemperaturaColocacion_C` y `TipoHormigon`, además del identificador y GlobalId. `ProgramaVertidas` toma esos valores de la metadata extraída del IFC, no de una lista de horarios incrustada en el código.

El caso actual contiene inicios a **0, 6, 12, 18 y 24 horas**. Son un plan docente de ejemplo, no una recomendación constructiva. La fecha calendario del plan sigue siendo demostrativa; aquí el reloj se expresa en segundos desde el inicio del plan.

La versión inicial requiere `Inicio_h`. Si falta, muestra un error en lugar de inventar un horario. Para un IFC que solo contenga fechas absolutas habrá que añadir una conversión explícita y trazable de `InicioPlanificado` a horas relativas.

## Qué ocurre al colocar una vertida

`EstadoConstruccion.activar_hasta(t)` identifica celdas nuevas y les asigna su temperatura de colocación **una sola vez**. No reinicia las temperaturas de las vertidas anteriores. Las celdas exteriores y futuras mantienen `NaN` como temperatura no definida: no se representan como hormigón frío que absorba calor.

Esta operación no calcula el calentamiento entre horarios. En una integración térmica posterior, se actualizarán las temperaturas de las celdas ya activas mediante el balance energético.

`partir_intervalo(t0,t1)` devuelve los eventos de colocación dentro del paso. Si un paso propuesto va de la hora 5 a la 7, se divide en 5→6 y 6→7 para incorporar la vertida exactamente a la hora 6.

## Diseño de hidratación y origen de los parámetros

Se adopta el **ejemplo pedagógico de las secciones 23–25** de `Antecedentes/guias/DIFUSION_DE_CALOR_1D_2D_3D.md`:

\[
\theta_{ad}(a)=\theta_\infty(1-e^{-a/\tau}),\qquad a\geq0
\]

La elevación adiabática es un aumento acumulado de temperatura, no una potencia. La fuente equivalente es:

\[
q_h(a)=\rho c_p\frac{d\theta_{ad}}{da}
=\rho c_p\frac{\theta_\infty}{\tau}e^{-a/\tau}
\quad[\mathrm{W/m^3}]
\]

Antes de colocar: fuente cero. En el instante de colocación empieza la curva del nuevo hormigón.

Los valores **θ∞=45 °C** y **τ=30 h** proceden del ejemplo de la guía, no de una calibración del hormigón del proyecto. Se declaran en `config/hidratacion_e07.yaml` con estado `pedagogico_no_calibrado`. La densidad y calor específico se leen de `config/caso_base.yaml`; en E07 se supone el mismo material para las mezclas configuradas.

Con ρ=2400 kg/m³ y cp=900 J/(kg·K), la fuente inicial es **900 W/m³**. Este modelo tiene su máximo al comenzar y luego decae; no incluye un período de inducción ni un pico retardado. Es una simplificación docente que se debe revisar al calibrar.

No se reutiliza el pulso de coseno con Qmax=1300 del ejercicio independiente: aquel ejercicio no define los parámetros de hidratación de esta fundación.

## Balance de energía

La energía generada por unidad de volumen entre edades a0 y a1 se integra exactamente:

\[
\Delta e_h=\rho c_p[\theta_{ad}(a_1)-\theta_{ad}(a_0)]\quad[\mathrm{J/m^3}]
\]

Las edades previas a colocación se recortan a cero. Así, si una vertida empieza a la hora 6 y se evalúa el intervalo global 5→7, solo genera energía durante 6→7.

En una celda aislada: ΔT=Δe_h/(ρcp). En una celda con conducción habrá que sumar el intercambio por sus caras. La energía volumétrica se multiplica por el volumen de celda para obtener joules.

## Alcance térmico y pasos siguientes

El modelo inicial depende de la edad real desde colocación, **no de la temperatura**. Una extensión futura es una edad equivalente o cinética Arrhenius, con parámetros justificados y calibrados. Debe incorporarse después de verificar el modelo simple y su balance energético. Cambiar la curva por otra no debería exigir reescribir la lógica de colocación.

El prototipo posterior `mcplan/solver_vertidas.py` conecta la difusión con la máscara activa mediante estas reglas:

- Cara entre dos celdas activas: conducción entre materiales, con sus resistencias.
- Cara entre una celda activa y una vertida futura: superficie expuesta, no contacto con hormigón aún inexistente.
- Al colocar la nueva vertida: esa superficie pasa a ser una interfaz entre hormigones.
- El exterior recibe condiciones ilustrativas por orientación. La clasificación detallada de suelo, ambiente y protecciones sigue pendiente.

No debe enviarse el arreglo con `NaN` directamente al solver actual, ni rellenarse todo el exterior con propiedades de hormigón.

## Ejecutar y mostrar

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/07_voxelizar_ifc.py
.\.venv\Scripts\python.exe -X utf8 scripts/08_vertidas_hidratacion.py
.\.venv\Scripts\python.exe -X utf8 -m pytest tests/test_construccion_hidratacion.py -v
```

El segundo script también se puede ejecutar con ▶ en VS Code y el intérprete `.venv`. Si el IFC cambió desde la generación de la grilla, solicita repetir la voxelización.

Se generan en `resultados/e07_construccion/`:

- `activacion_vertidas_cortes.png`: hormigón colocado por colores y vertidas futuras en gris.
- `hidratacion_y_activacion.png`: fuentes por edad, referencias adiabáticas por vertida y volumen colocado.
- `programa_vertidas.csv`: horarios y parámetros obtenidos de los datos.
- `historias_adiabaticas_por_vertida.csv`: fuentes y referencias energéticas por vertida.
- `estructura_estado_ultima_colocacion.npz`: geometría, máscara activa, temperaturas de colocación asignadas, edades y fuentes al último evento. No es un resultado de difusión.
- `resumen.json`: alcance, eventos y trazabilidad de los datos utilizados.

Las curvas de temperatura adiabática se calculan suponiendo **cada vertida aislada e independiente**. No son las temperaturas previstas del núcleo de la fundación.

## Qué explicar en el oral

“Conservamos la identidad de cada vóxel del IFC y agregamos una máscara que cambia con el tiempo. Cuando llega el horario de una vertida, asignamos su temperatura inicial sin reiniciar las anteriores. Su fuente se calcula con la edad desde colocación. La fuente proviene de derivar una curva adiabática; los parámetros actuales son los pedagógicos de la guía y faltan por calibrar.”
