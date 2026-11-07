"""Ejecuta el baseline E09 y extrae historias puntuales del corte guardado.

Desde la raiz: python -X utf8 scripts/12_baseline_e09.py
Para repetir solo las figuras: agregar --solo-figuras.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mcplan.voxelizacion import GrillaVoxeles

RAIZ = Path(__file__).resolve().parents[1]
SALIDA = RAIZ / 'resultados/e09'


def sha256(ruta):
    return hashlib.sha256(ruta.read_bytes()).hexdigest()


def leer_json(ruta):
    return json.loads(ruta.read_text(encoding='utf-8'))


def leer_csv(ruta):
    with ruta.open(encoding='utf-8-sig', newline='') as archivo:
        return list(csv.DictReader(archivo))


def escribir_csv(ruta, filas):
    with ruta.open('w', encoding='utf-8-sig', newline='') as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=list(filas[0]))
        escritor.writeheader()
        escritor.writerows(filas)


def ejecutar(script, *argumentos):
    subprocess.run([sys.executable, '-X', 'utf8', str(RAIZ / 'scripts' / script),
                    *map(str, argumentos)], cwd=RAIZ, check=True)


def figuras():
    voxeles = SALIDA / 'voxelizacion'
    simulacion = SALIDA / 'simulacion'
    grilla = GrillaVoxeles.cargar(voxeles / 'grilla_vertidas.npz')
    resumen = leer_json(simulacion / 'resumen.json')
    geometria = leer_json(voxeles / 'resumen.json')
    if sha256(voxeles / 'grilla_vertidas.npz') != resumen['grilla_sha256']:
        raise ValueError('La simulacion no corresponde a esta grilla.')
    if sha256(SALIDA / 'entradas/fundacion.ifc') != geometria['ifc_sha256']:
        raise ValueError('La copia del IFC cambio despues de voxelizar.')

    with np.load(simulacion / 'resultado_termico.npz', allow_pickle=False) as archivo:
        tiempos = archivo['tiempo_s'] / 3600
        cortes = archivo['cortes_xz_C'].copy()
        y_corte = float(archivo['y_corte_m'])
        final = archivo['T_final_C'].copy()
        ids = archivo['vertida_id'].copy()
    np.testing.assert_array_equal(ids, grilla.vertida_id)
    iy = int(np.argmin(abs(grilla.centros[1] - y_corte)))
    # Los puntos se definen en metros, se informa el centro real de la celda.
    pedidos = [('P1 base interior', 0., .275),
               ('P2 interior transicion', 0., 1.675),
               ('P3 interior pedestal', 0., 3.825),
               ('P4 cerca de coronacion', 0., 4.275),
               ('P5 transicion lateral', 6., 1.675)]
    sensores = []
    historias = []
    colores = ['#287bb5', '#d97922', '#258b62', '#985bad', '#cd4c44']
    figura, ejes = plt.subplots(2, 1, figsize=(10, 8), layout='constrained')
    balance = leer_csv(simulacion / 'balance_energia.csv')
    ejes[0].plot(tiempos / 24, [float(fila['T_max_C']) for fila in balance], color='#176e60', lw=2)
    ejes[0].scatter(resumen['tiempo_maximo_h'] / 24, resumen['maximo_global_C'], color='#cd4c44', zorder=3)
    ejes[0].axhline(70, color='#9e6662', ls='--', label='Limite docente 70 C')
    ejes[0].set(title='Maxima global del dominio activo', ylabel='Temperatura [C]')
    ejes[0].legend()
    for (nombre, x, z), color in zip(pedidos, colores):
        ix = int(np.argmin(abs(grilla.centros[0] - x)))
        iz = int(np.argmin(abs(grilla.centros[2] - z)))
        vid = int(ids[ix, iy, iz])
        if vid == 0:
            raise ValueError(f'{nombre} queda fuera del hormigon.')
        valores = cortes[:, iz, ix].astype(float)
        inicio = float(grilla.vertidas[vid]['Inicio_h'])
        if np.any(np.isfinite(valores[tiempos < inicio])):
            raise ValueError('Una historia contiene temperatura antes de la colocacion.')
        indice = int(np.nanargmax(valores))
        punto = dict(punto=nombre, x_m=float(grilla.centros[0][ix]), y_m=y_corte,
                     z_m=float(grilla.centros[2][iz]), vertida_id=vid, inicio_h=inicio,
                     T_max_C=float(valores[indice]), hora_max_h=float(tiempos[indice]),
                     T_final_C=float(valores[-1]))
        sensores.append(punto)
        for hora, valor in zip(tiempos, valores):
            historias.append(dict(tiempo_h=float(hora), punto=nombre,
                                  T_C=float(valor) if np.isfinite(valor) else None))
        ejes[1].plot(tiempos / 24, valores, label=nombre, color=color)
    for eje in ejes:
        eje.set_xlabel('Tiempo desde primera vertida [dias]')
        eje.grid(alpha=.2)
    ejes[1].set(title='Temperatura en centros de celdas fijos', ylabel='Temperatura [C]')
    ejes[1].legend(fontsize=8, ncol=2)
    figura.savefig(SALIDA / 'historias_puntos_y_maxima.png', dpi=160)
    plt.close(figura)
    escribir_csv(SALIDA / 'puntos_resumen.csv', sensores)
    escribir_csv(SALIDA / 'historias_puntos.csv', historias)
    # Comparacion simultanea de dos puntos definidos. No es maximo-minimo
    # global ni una busqueda completa de pares nucleo-superficie.
    p2 = sensores[1]
    p4 = sensores[3]
    ix2 = int(np.argmin(abs(grilla.centros[0] - p2['x_m'])))
    iz2 = int(np.argmin(abs(grilla.centros[2] - p2['z_m'])))
    ix4 = int(np.argmin(abs(grilla.centros[0] - p4['x_m'])))
    iz4 = int(np.argmin(abs(grilla.centros[2] - p4['z_m'])))
    diferencia = abs(cortes[:, iz2, ix2].astype(float) - cortes[:, iz4, ix4].astype(float))
    indice_diferencial = int(np.nanargmax(diferencia))
    diagnostico = dict(par='P2 interior transicion - P4 cerca de coronacion',
                       maximo_C=float(diferencia[indice_diferencial]),
                       tiempo_h=float(tiempos[indice_diferencial]),
                       alcance='Par fijo del corte; falta evaluar otros pares y sensibilidad espacial.')

    figura, ejes = plt.subplots(2, 2, figsize=(12, 7), layout='constrained')
    horas_mostrar = [24., 96., resumen['tiempo_maximo_h'], 336.]
    extension = (grilla.origen[0], grilla.origen[0] + ids.shape[0] * grilla.paso[0],
                 grilla.origen[2], grilla.origen[2] + ids.shape[2] * grilla.paso[2])
    mapa = plt.get_cmap('inferno').copy()
    mapa.set_bad('#e5e7e9')
    for eje, hora in zip(ejes.flat, horas_mostrar):
        indice = int(np.argmin(abs(tiempos - hora)))
        imagen = eje.imshow(np.ma.masked_invalid(cortes[indice]), origin='lower', extent=extension,
                            vmin=15, vmax=65, cmap=mapa, interpolation='nearest', aspect='equal')
        eje.set(title=f't = {tiempos[indice]:g} h ({tiempos[indice]/24:.2f} dias)', xlabel='x [m]', ylabel='z [m]')
        for sensor, color in zip(sensores, colores):
            eje.scatter(sensor['x_m'], sensor['z_m'], s=16, facecolor=color, edgecolor='white', linewidth=.5)
    figura.colorbar(imagen, ax=list(ejes.flat), label='Temperatura [C]', shrink=.75)
    figura.suptitle(f'Cortes del dominio 3D en y = {y_corte:.2f} m; gris = exterior o sin colocar')
    figura.savefig(SALIDA / 'cortes_con_puntos.png', dpi=160)
    plt.close(figura)

    # Tres planos ortogonales del campo final, en sus coordenadas fisicas.
    figura = plt.figure(figsize=(10, 7), layout='constrained')
    eje = figura.add_subplot(111, projection='3d')
    ejes_indice = [int(np.argmin(abs(grilla.centros[0]))), iy,
                  int(np.argmin(abs(grilla.centros[2] - 1.675)))]
    mapa = plt.get_cmap('inferno')
    norma = matplotlib.colors.Normalize(vmin=15, vmax=65)
    for normal, indice in enumerate(ejes_indice):
        variables = [direccion for direccion in range(3) if direccion != normal]
        a, b = np.meshgrid(grilla.centros[variables[0]], grilla.centros[variables[1]], indexing='ij')
        coordenadas = [None, None, None]
        coordenadas[normal] = np.full(a.shape, grilla.centros[normal][indice])
        coordenadas[variables[0]], coordenadas[variables[1]] = a, b
        campo = np.take(final, indice, axis=normal)
        colores_campo = mapa(norma(np.nan_to_num(campo, nan=15)))
        colores_campo[~np.isfinite(campo), 3] = 0
        eje.plot_surface(*coordenadas, facecolors=colores_campo, shade=False,
                         rstride=1, cstride=1, linewidth=0, antialiased=False)
    eje.set(xlabel='x [m]', ylabel='y [m]', zlabel='z [m]', title='Tres cortes del campo final a los 14 dias')
    eje.set_box_aspect((18, 18, 4.3))
    eje.view_init(24, -55)
    figura.colorbar(matplotlib.cm.ScalarMappable(norm=norma, cmap=mapa), ax=eje, shrink=.6, label='Temperatura [C]')
    figura.savefig(SALIDA / 'cortes_3d_final.png', dpi=160)
    plt.close(figura)
    manifiesto = dict(resumen=resumen, volumen=geometria, sensores=sensores,
                      diferencial_diagnostico=diagnostico,
                      archivos={str(ruta.relative_to(SALIDA)): sha256(ruta)
                                for ruta in sorted((SALIDA / 'entradas').iterdir())})
    (SALIDA / 'evidencia.json').write_text(json.dumps(manifiesto, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(sensores, indent=2, ensure_ascii=False), flush=True)
    print(json.dumps(diagnostico, indent=2, ensure_ascii=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solo-figuras', action='store_true')
    args = parser.parse_args()
    if not args.solo_figuras:
        entradas = SALIDA / 'entradas'
        entradas.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(RAIZ / 'fundaciones/fundacion.ifc', entradas / 'fundacion.ifc')
        for nombre in ['caso_base.yaml', 'hidratacion_e07.yaml', 'termica_e07.yaml']:
            shutil.copyfile(RAIZ / 'config' / nombre, entradas / nombre)
        ejecutar('07_voxelizar_ifc.py', '--ifc', entradas / 'fundacion.ifc',
                 '--dx', .2, '--dy', .2, '--dz', .05, '--salida', SALIDA / 'voxelizacion')
        ejecutar('09_simular_fundacion_vertidas.py', '--grilla', SALIDA / 'voxelizacion/grilla_vertidas.npz',
                 '--termica', entradas / 'termica_e07.yaml', '--hidratacion', entradas / 'hidratacion_e07.yaml',
                 '--horizonte-h', 336, '--max-dt-s', 600, '--salida', SALIDA / 'simulacion')
    figuras()


if __name__ == '__main__':
    main()
