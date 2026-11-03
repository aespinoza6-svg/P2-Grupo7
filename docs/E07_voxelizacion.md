# E07: primer voxelizador de la fundación IFC

## Qué hace

El código lee la geometría del IFC, crea una grilla cartesiana y asigna una etiqueta a cada celda. No reconstruye la fundación a partir de las dimensiones del YAML: si se cambia la geometría y se guarda el IFC en Bonsai, el voxelizador lee esa nueva forma.

- `mcplan/voxelizacion.py`: extracción, clasificación geométrica, etiquetas y volumen.
- `scripts/07_voxelizar_ifc.py`: ejecución y generación de evidencia.
- `tests/test_voxelizacion.py`: pruebas geométricas.

Cada posición de `vertida_id[i,j,k]` corresponde a una celda en los ejes **(x,y,z)**, la misma convención que usa el solver 3D:

- **0:** fuera del hormigón. Todavía no significa automáticamente aire o suelo.
- **ID positivo:** identificador real de la vertida obtenido del IFC.

La etiqueta geométrica no cambia cuando una vertida aún no ha sido colocada. Más adelante se agregará una máscara temporal separada para representar su activación.

## Cómo decide qué celdas pertenecen al hormigón

1. IfcOpenShell triangula cada sólido y transforma los vértices a coordenadas globales en metros.
2. Se verifica que la malla sea cerrada y que los triángulos estén orientados consistentemente.
3. La grilla se ancla al mínimo de la caja global. Cada celda conserva el tamaño elegido; la última puede sobrepasar la caja del IFC.
4. Para cada centro candidato se suma el ángulo sólido orientado de los triángulos. Dentro de una superficie cerrada el total es aproximadamente ±4π; fuera es aproximadamente cero.
5. Se asigna el ID de la vertida que contiene el centro. Las cajas envolventes solo descartan candidatos; no se usan como sustituto de la forma real.

El método admite superficies inclinadas, cilindros faceteados y cavidades con sus superficies correctamente orientadas. Se utiliza la triangulación que entrega IfcOpenShell; por tanto hay también una aproximación geométrica de las superficies curvas.

## Juntas y traslapes

- Un centro estrictamente dentro de dos vertidas causa un error: no se oculta el traslape.
- Si un centro está en una superficie compartida, prevalece el interior estricto de una vertida sobre la superficie de otra.
- Si solo está sobre ambas superficies, se asigna al ID menor, de forma determinista. Nunca se cuenta dos veces.
- Se registra cuántos centros compartidos hubo. Esto permite identificar una malla cuya resolución resulta desfavorable para las juntas.

## Elección inicial de grilla

Por defecto: **Δx=Δy=0,20 m y Δz=0,05 m**. Es una grilla cartesiana regular aunque las celdas no sean cúbicas.

Las cotas de las cinco vertidas son múltiplos de 5 cm. Esta elección alinea las juntas horizontales con caras de celda. Una grilla cúbica de 20 cm es más pequeña, pero algunas juntas caen sobre centros y puede distorsionar los volúmenes por vertida. La elección se declara explícitamente; no se modifica la geometría para acomodarla a la grilla.

El volumen aproximado de cada vertida es:

`V_voxel = número_de_celdas_de_la_vertida × Δx × Δy × Δz`

Se compara con el volumen de la geometría triangulada del IFC, por vertida y total. Es una comprobación geométrica; no una validación experimental de la física térmica.

La superficie queda escalonada. Reducir el paso suele mejorar la representación, pero el error no tiene por qué bajar monótonamente para cada malla o vertida. La comparación mantiene Δz fijo y refina x/y, por lo que mide sensibilidad a esa resolución horizontal, no convergencia isotrópica completa.

## Cómo ejecutar

Desde la raíz del repositorio, después de instalar `requirements.txt`:

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/07_voxelizar_ifc.py
```

Para comparar las resoluciones horizontales de 40, 20 y 10 cm, con Δz=5 cm:

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/07_voxelizar_ifc.py --comparar
```

Para otro IFC o tamaño de celda:

```powershell
.\.venv\Scripts\python.exe -X utf8 scripts/07_voxelizar_ifc.py --ifc fundaciones/fundacion.ifc --dx .2 --dy .2 --dz .1
```

El script también se puede abrir y ejecutar en VS Code con el intérprete `.venv`; sus entradas predeterminadas se resuelven respecto del repositorio, no de la carpeta actual de la terminal.

## Archivos de salida

Se guardan en `resultados/e07_voxelizacion/`:

- `grilla_vertidas.npz`: origen, pasos, arreglo `vertida_id` y propiedades JSON de las vertidas, incluyendo `GlobalId`, `Inicio_h` y temperatura de colocación. No usa pickle.
- `grilla_vertidas.vtk`: etiquetas como datos por celda; abrir en ParaView o PyVista y filtrar `vertida_id >= 1` para ver solo hormigón. Las dimensiones VTK son de nodos y las etiquetas son `CELL_DATA`.
- `volumenes_por_vertida.csv`: conteos, volúmenes y errores relativos firmados.
- `resumen.json`: parámetros, huella SHA-256 del IFC leído y resultados.
- `voxelizacion_cortes.png`: cortes donde se ve el exterior y la identidad de las vertidas.
- `voxelizacion_3d.png`: vista de centros de celdas de la superficie, coloreados por ID. Es una vista muestreada para facilitar la lectura; los datos completos están en NPZ/VTK.
- `comparacion_mallas.csv`: solo al ejecutar `--comparar`.

## Pruebas y alcance

```powershell
.\.venv\Scripts\python.exe -X utf8 -m pytest tests/test_voxelizacion.py -v
```

Las pruebas cubren interior/exterior/superficie, una cavidad, una caja rotada y trasladada, juntas exactas, traslapes, mallas abiertas, pasos inválidos, límites de memoria, persistencia NPZ sin pickle y las cinco vertidas del IFC real.

Este avance **no** incluye todavía hidratación, suelo, activación temporal o condiciones térmicas en la superficie voxelizada. No se debe pasar toda la caja al solver como si fuera hormigón: primero se necesita tratar explícitamente la máscara `vertida_id > 0` y sus superficies expuestas.
