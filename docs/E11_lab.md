# E11: guion del lab de clima y visualización 3D

Jueves 12 de noviembre, 10:30–12:30.

Objetivo: conectar datos meteorológicos reales con las condiciones de borde y producir resultados 3D que se puedan inspeccionar.

## Antes de llegar

Desde `P2-Grupo7`, con el entorno activado:

```powershell
python -X utf8 scripts/14_simular_ambiente.py --horizonte-h 24 --vtk 12,24 --salida resultados/e11/clima_real_24h
```

Tarda unos 20 s y deja dos archivos:

- `resultados/e11/clima_real_24h/campo_12h.vtk`;
- `resultados/e11/clima_real_24h/campo_24h.vtk`.

Git no los sube porque pesan unos 15 MB cada uno: cada integrante los genera en su computador. Para abrirlos con PyVista, instálalo una vez:

```powershell
pip install -r requirements-visual.txt
```

Si quieres un pronóstico del día, corre `python -X utf8 scripts/13_descargar_pronostico.py` y después `python -X utf8 scripts/actualizar_visor_laboratorio.py`.

## Los 8 pedidos del ayudante

| # | Pedido | En el HTML (`visor_laboratorio.html`) | En Python |
|---|---|---|---|
| 1 | Mostrar el archivo meteorológico | **Ambiente → Clima y ambiente → «Ver archivo y metadata del snapshot»**: tabla de las primeras 24 h | Abrir `data/weather/2026-11-10T1532/forecast.json` |
| 2 | Mostrar la metadata | En el mismo panel, el bloque de metadata | `data/weather/2026-11-10T1532/metadata.json` |
| 3 | Graficar la temperatura ambiente | Pestaña **Gráficos → Clima usado**, con «Temperatura del aire a 2 m» | `resultados/e10/clima_real/clima.png` |
| 4 | Alimentar una BC con esa serie | Botón **«Usar clima real, moldaje y suelo»**. Las caras quedan de tipo **Ambiente** | `config/ambiente_e10.yaml` y `mcplan/bordes.py` |
| 5 | Ejecutar una simulación corta | Horizonte de 24 h en **Simulación** y **Calcular**, unos 4 s | El comando de «Antes de llegar» |
| 6 | Abrir resultados 3D | Pestaña **Vista 3D**, con animación en el tiempo | `python -X utf8 scripts/15_ver_vtk.py resultados/e11/clima_real_24h/campo_24h.vtk`, o ParaView |
| 7 | Hacer un slice o un clipping | **Clipping:** «Cortar en y = 0» en la Vista 3D. **Slice:** **Gráficos → Corte térmico**, con el plano XZ, YZ o XY y su posición | El script 15 muestra los dos: clip a la izquierda y slice con isotermas a la derecha. `--plano x --pos 2.5` cambia el plano |
| 8 | Identificar núcleo y superficie | **Vista 3D**: punto celeste = núcleo, verde = superficie. Debajo están sus coordenadas y temperaturas | El script 15 los marca y los imprime |

### Qué decir en cada paso

1. **Archivo:** es la respuesta de Open-Meteo guardada tal cual. Tiene horas en hora local y tres variables: temperatura a 2 m, viento a 10 m y radiación solar.
2. **Metadata:** dice de dónde salió el dato, con qué modelo, cuándo se descargó, las coordenadas y las unidades. Pedimos las coordenadas del curso (−37,6921; −72,624), y el modelo respondió con su celda más cercana (−37,645; −72,577). Es normal: el modelo tiene celdas de varios km.
3. **Gráfico:** la temperatura sube de día y baja de noche. Lo sombreado son horas que el pronóstico no cubre y que se completaron con histórico ERA5.
4. **BC:** cada cara al aire usa Robin, $q = (T_{aire}(t) - T)/(\Delta x/2k + R + 1/h)$. $T_{aire}(t)$ sale del snapshot, interpolado hora a hora.
5. **Simulación corta:** 24 h, con las vertidas entrando según el IFC.
6. **3D:** cada cubo es una celda; el color es su temperatura.
7. **Slice y clip:** el **clip** corta el sólido y deja ver el interior; el **slice** muestra solo el plano de corte.
8. **Núcleo y superficie:** el núcleo es la celda más caliente de la vertida con mayor diferencia. La superficie es su celda más fría con una cara al aire. Su diferencia es el diferencial que se compara con 20 °C.

## Si preguntan por viento o radiación: dónde entran

Las dos entran en la misma fórmula del borde, en `mcplan/bordes.py`:

$$q = \frac{T_{ref}(t) - T_{celda}}{\dfrac{\Delta x}{2k} + R(t) + \dfrac{1}{h(t)}}\,A$$

- **Viento V [m/s] → en h:** $h = h_{conv} + h_{rad}$, con $h_{conv} = 5{,}7 + 3{,}8\,V$ (función `h_conveccion`) y $h_{rad} = 4\varepsilon\sigma T^3$ (función `h_radiacion`). Más viento significa más h y más pérdida.
- **Radiación solar G [W/m²] → en $T_{ref}$**, solo en las caras que miran al cielo: $T_{ref} = T_{aire} + \alpha G/h$, la temperatura sol-aire. Está en `BordeAmbiental.coeficientes`. Así el sol se suma una sola vez y no se cuenta doble con la convección.

## Abrir el VTK en ParaView, si está disponible

1. **File → Open** y elegir `campo_24h.vtk`.
2. **Filters → Threshold**, con `activa` entre 1 y 1.
3. Color por `temperatura_C`.
4. **Filters → Clip** (plano Y) o **Filters → Slice**.
5. Para la superficie: **Threshold** con `superficie` = 1.

| Campo del VTK | Qué es |
|---|---|
| `temperatura_C` | Temperatura de la celda; vacía si no existe todavía |
| `vertida_id` | 1 a 5 para las vertidas y 90 para el suelo |
| `activa` | 1 si la celda ya está colocada en ese instante |
| `superficie` | 1 si es hormigón con una cara al aire |
