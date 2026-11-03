# Primera simulación térmica del IFC construido por vertidas

## Estado y propósito

Este es un **prototipo numérico acoplado**, posterior a la estructura mínima E07. Conecta el IFC voxelizado, los horarios, las temperaturas de colocación, el modelo pedagógico de hidratación y la difusión 3D.

No es un Plan de Control de Temperatura ni una predicción validada: siguen pendientes la calibración del hormigón, el clima real, el suelo explícito, las protecciones y su retiro. No se evalúa todavía el diferencial núcleo–superficie a 25 mm ni confiabilidad.

## Archivos

- `mcplan/solver_vertidas.py`: balance en la geometría activa y eventos de colocación.
- `config/termica_e07.yaml`: ambiente y base de la demostración.
- `config/hidratacion_e07.yaml`: curva pedagógica, no calibrada.
- `scripts/09_simular_fundacion_vertidas.py`: ejecución y gráficos.
- `tests/test_solver_vertidas.py`: pruebas del acoplamiento.

El solver previo `solver3d.py` se conserva para los casos de caja y las verificaciones de semanas anteriores.

## Qué ocurre en cada paso

1. Incorporar las vertidas cuyo horario corresponde al tiempo actual.
2. Asignar su temperatura de colocación, manteniendo las temperaturas de las anteriores.
3. Determinar las caras entre celdas activas y las caras expuestas al exterior o a vertidas futuras.
4. Actualizar el límite temporal de estabilidad con las conductancias del dominio activo.
5. Calcular conducción y flujos de borde a partir de las temperaturas al inicio del paso.
6. Sumar la energía de hidratación del intervalo según la edad de cada vertida.
7. Actualizar exclusivamente las celdas activas.

El paso termina exactamente en cada nueva colocación. No se utiliza una vertida antes de su horario. Los tiempos solicitados de salida pueden acortar adicionalmente un paso.

## Cara expuesta e interfaz

Una cara entre hormigón colocado y una vertida futura recibe la condición exterior correspondiente a su orientación. Cuando se coloca la vecina, deja de ser una superficie exterior y se convierte en una interfaz de conducción. Se supone contacto perfecto, sin una resistencia adicional de junta.

El exterior y el hormigón aún no colocado tienen temperatura `NaN` y no aportan capacidad ni participan en conducción. El código emplea valores auxiliares finitos solo para operaciones vectorizadas; sus conductancias están anuladas donde no hay dos celdas activas. No se convierten en bloques fríos ficticios.

## Balance por celda y estabilidad

Para una celda activa:

\[
C_i\frac{dT_i}{dt}=\sum_j G_{ij}(T_j-T_i)+P_{borde,i}+q_{h,i}V_i,
\qquad C_i=\rho_i c_{p,i}V_i
\]

En una interfaz entre dos materiales:

\[
G_{ij}=\frac{A}{d/(2k_i)+d/(2k_j)}
\]

En Robin:

\[
G_b=\frac{A}{d/(2k_i)+R_{extra}+1/h},\qquad
P_b=G_b(T_{aire}-T_i)
\]

Neumann positivo representa potencia entrante `q A`. Las conductancias contribuyen al criterio suficiente de coeficientes positivos `dt <= min(C_i/S_i)` para las celdas activas. Neumann prescrito no suma una conductancia dependiente de temperatura.

La selección automática usa el menor entre el máximo de configuración y el 80 % del límite. Un `dt_s` explícito excesivo provoca un error. Se usa Euler explícito para conducción y bordes, con energía de hidratación integrada analíticamente dentro de cada paso.

## Balance global con material nuevo

En un dominio que crece, no basta con comparar energía almacenada y calor de hidratación. El hormigón nuevo introduce energía interna al colocarse:

\[
E(t)=E_{colocación}(t)+E_{hidratación}(t)+E_{bordes}(t)
\]

Todas las energías se refieren a 0 °C con propiedades constantes. `E_colocación` suma `rho cp V T_colocación` para cada celda incorporada; `E_hidratación` suma la energía generada; `E_bordes` es positiva si entra calor neto. La conducción entre celdas se cancela al sumar sus balances.

El residuo global se guarda para comprobar que no se crea ni pierde energía por errores de implementación. Eso es verificación numérica; no calibración o validación experimental.

## Supuestos de la primera demostración

- Aire constante a **15 °C** y `h=10 W/(m² K)` en caras laterales y superiores expuestas. Son valores ilustrativos, no datos meteorológicos del sitio.
- Las caras expuestas de orientación `z-` se consideran adiabáticas: es una simplificación de la base, no un modelo del suelo. Con la geometría actual no se representan voladizos.
- Contacto perfecto entre hormigones, sin moldajes, aislación o retiros.
- Hidratación exponencial de la guía, independiente de la temperatura, con elevación final 45 °C y tau 30 h.
- Malla cartesiana escalonada y colocación instantánea de cada vertida en su horario.
- Por defecto se demuestran **72 horas**, configurables. El horizonte final del proyecto será de 14 días.

El programa lee los horarios y temperaturas de las vertidas del IFC voxelizado. Si el IFC cambió desde la voxelización, solicita regenerar la grilla.

## Ejecutar

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py
```

Para extender a 14 días:

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/09_simular_fundacion_vertidas.py --horizonte-h 336
```

Pruebas:

```powershell
.\.venv\Scripts\python.exe -X utf8 -m pytest tests/test_solver_vertidas.py -v
```

## Resultados

En `resultados/e07_simulacion_acoplada/`:

- `temperaturas_por_vertida.png`: máximas y medias calculadas, y volumen activo.
- `cortes_termicos.png`: campos de temperatura con zonas no activas o exteriores en gris.
- `balance_energia.png` y `balance_energia.csv`: términos del balance y residuo.
- `historias_por_vertida.csv`: temperaturas mínima, media y máxima por vertida.
- `resultado_termico.npz`: temperatura final 3D, máscara activa, identidades, tiempos y cortes.
- `resumen.json`: máximo global, instante, parámetros y limitaciones.

Por defecto solo se guardan cortes durante la evolución y el campo 3D final para no almacenar centenares de campos completos. La API permite historiales 3D pequeños cuando se solicitan explícitamente.

## Verificaciones

Las pruebas comparan con el solver de caja cuando todo el dominio existe desde el inicio; comprueban que huecos y vertidas futuras no absorban calor; verifican las interfaces al activar, la fuente adiabática exacta, el aporte `qAt`, dos materiales, el rechazo de pasos inestables y el balance con incorporación de hormigón.

## Comparación de pasos y corrida de 14 días

El argumento `--max-dt-s` permite cambiar el máximo temporal sin editar el solver. El criterio de estabilidad puede reducirlo automáticamente. `--salida` permite conservar dos corridas comparables en carpetas fijas distintas.

Las ejecuciones de 14 días y de sensibilidad temporal de 72 h se documentan, con sus resultados y limitaciones, en [E07_verificacion_acoplada.md](E07_verificacion_acoplada.md). El script `scripts/10_comparar_pasos_e07.py` compara resultados guardados sin repetir las simulaciones.
