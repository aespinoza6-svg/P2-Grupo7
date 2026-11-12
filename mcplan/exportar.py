"""Exporta campos de la grilla a VTK para abrirlos en ParaView o PyVista.

Formato: VTK legacy BINARY, STRUCTURED_POINTS (una grilla cartesiana regular).
Cada celda del arreglo (x, y, z) es una celda VTK; x cambia mas rapido.
Campos por celda (CELL_DATA):
- temperatura_C: NaN fuera del hormigon/suelo o en vertidas aun no colocadas;
- vertida_id: 0 exterior, 1..5 vertidas, 90 suelo;
- activa: 1 si la celda existe en ese instante;
- superficie: 1 si es hormigon activo con una cara al aire (x±, y±, z+).
En ParaView: Threshold sobre activa = 1, luego Slice o Clip. En PyVista: scripts/15_ver_vtk.py.
"""
from pathlib import Path
import numpy as np

TIPOS = {np.dtype('float32'): 'float', np.dtype('int32'): 'int', np.dtype('uint8'): 'unsigned_char'}


def celdas_superficie(activa, vertida_id, excluir=(90,)):
    """Hormigon activo con al menos una cara hacia una celda inactiva o fuera de la grilla, salvo hacia abajo."""
    hormigon = activa & ~np.isin(vertida_id, excluir)
    vecino = np.pad(activa, 1, constant_values=False)
    sup = np.zeros(activa.shape, dtype=bool)
    nx, ny, nz = activa.shape
    for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1)):
        sup |= ~vecino[1 + dx:1 + dx + nx, 1 + dy:1 + dy + ny, 1 + dz:1 + dz + nz]
    return hormigon & sup


def escribir_vtk(ruta, origen, paso, campos: dict, titulo='mcplan'):
    """Escribe campos de forma (nx, ny, nz) como CELL_DATA de STRUCTURED_POINTS."""
    forma = next(iter(campos.values())).shape
    if any(np.shape(v) != forma for v in campos.values()) or len(forma) != 3:
        raise ValueError('Todos los campos deben tener la misma forma 3D.')
    nx, ny, nz = forma
    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open('wb') as f:
        cabeza = ['# vtk DataFile Version 3.0', titulo.replace('\n', ' ')[:250], 'BINARY', 'DATASET STRUCTURED_POINTS',
                  f'DIMENSIONS {nx + 1} {ny + 1} {nz + 1}', 'ORIGIN ' + ' '.join(f'{x:.9g}' for x in origen),
                  'SPACING ' + ' '.join(f'{x:.9g}' for x in paso), f'CELL_DATA {nx * ny * nz}']
        f.write(('\n'.join(cabeza) + '\n').encode('ascii'))
        for nombre, valores in campos.items():
            arr = np.asarray(valores)
            if arr.dtype == bool:
                arr = arr.astype(np.uint8)
            elif np.issubdtype(arr.dtype, np.floating):
                arr = arr.astype(np.float32)
            elif np.issubdtype(arr.dtype, np.integer) and arr.dtype != np.uint8:
                arr = arr.astype(np.int32)
            tipo = TIPOS[arr.dtype]
            f.write(f'SCALARS {nombre} {tipo} 1\nLOOKUP_TABLE default\n'.encode('ascii'))
            f.write(arr.ravel(order='F').astype(arr.dtype.newbyteorder('>')).tobytes())
            f.write(b'\n')
    return ruta


def campo_vtk(ruta, grilla, T, activa, t_h, extra=''):
    """Atajo para un instante de la simulacion por vertidas."""
    temp = np.where(activa, T, np.nan).astype(np.float32)
    campos = dict(temperatura_C=temp, vertida_id=grilla.vertida_id.astype(np.int32), activa=activa,
                  superficie=celdas_superficie(activa, grilla.vertida_id))
    return escribir_vtk(ruta, grilla.origen, grilla.paso, campos, f'mcplan t={t_h:g} h {extra}')


def leer_vtk(ruta):
    """Lector minimo de lo que escribe escribir_vtk (para tests, sin VTK instalado)."""
    datos = Path(ruta).read_bytes()
    i = 0
    def linea():
        nonlocal i
        j = datos.index(b'\n', i); s = datos[i:j].decode('ascii'); i = j + 1; return s
    for _ in range(4):
        linea()
    dims = [int(x) - 1 for x in linea().split()[1:]]
    origen = [float(x) for x in linea().split()[1:]]
    paso = [float(x) for x in linea().split()[1:]]
    n = int(linea().split()[1])
    inv = {v: k for k, v in TIPOS.items()}
    campos = {}
    while i < len(datos):
        partes = linea().split()
        if not partes:
            continue
        _, nombre, tipo, _ = partes
        linea()
        dt = inv[tipo].newbyteorder('>')
        arr = np.frombuffer(datos, dtype=dt, count=n, offset=i)
        i += n * dt.itemsize + 1
        campos[nombre] = arr.astype(dt.newbyteorder('=')).reshape(dims, order='F')
    return dict(forma=dims, origen=origen, paso=paso, campos=campos)
