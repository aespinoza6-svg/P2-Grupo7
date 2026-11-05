"""Lab: recorrido IFC -> voxeles -> simulacion por etapas, con cambios en vivo.

Desde la raiz y con el entorno instalado:
python -X utf8 scripts/11_lab_vertidas.py --solo-voxelizar --dx .1 --dy .1
python -X utf8 scripts/11_lab_vertidas.py --solo-voxelizar --ocultar V05
python -X utf8 scripts/11_lab_vertidas.py --inicio V05=36 --temp V05=10
python -X utf8 scripts/11_lab_vertidas.py --ifc-cambiado fundaciones/fundacion_lab.ifc

La idea que demuestra: la secuencia de vertidas esta en el IFC, no en el solver.
Este script nunca edita el IFC original ni el solver: copia el IFC, cambia las
propiedades en la copia y vuelve a correr los mismos scripts 07 y 09.
Todo queda en resultados/lab/, que git ignora.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mcplan import bim
from mcplan.voxelizacion import GrillaVoxeles

RAIZ = Path(__file__).resolve().parents[1]
COLORES = ['#287bb5', '#d97922', '#258b62', '#985bad', '#cd4c44']  # mismos de los scripts 07 y 09
CONFIGURACION = ['caso_base.yaml', 'termica_e07.yaml', 'hidratacion_e07.yaml']


# ---------- Leer lo que pide el usuario ----------

def leer_vertida(texto):
    """Convierte V02, v2 o 2 en el numero 2."""
    try:
        return bim._vertida_id_numero(texto)
    except ValueError:
        raise argparse.ArgumentTypeError(f'Vertida invalida: {texto!r}. Ejemplo: V02.') from None


def leer_cambio(texto):
    """Convierte 'V02=24' en (2, 24.0)."""
    if '=' not in texto:
        raise argparse.ArgumentTypeError(f'Falta el signo =: {texto!r}. Ejemplo: V02=24.')
    vertida, valor = texto.split('=', 1)
    try:
        numero = float(valor)
    except ValueError:
        raise argparse.ArgumentTypeError(f'El valor de {texto!r} no es un numero.') from None
    if not math.isfinite(numero):
        raise argparse.ArgumentTypeError(f'El valor de {texto!r} no es finito.')
    return leer_vertida(vertida), numero


def nombre_escenario(args):
    """Nombre de carpeta que dice que se cambio, por ejemplo V02_inicio24__V02_temp10."""
    if args.ifc_cambiado:
        return args.ifc_cambiado.stem
    partes = [f'V{v:02d}_inicio{x:g}' for v, x in args.inicio]
    partes += [f'V{v:02d}_temp{x:g}' for v, x in args.temp]
    return '__'.join(partes)


# ---------- Correr el pipeline, reutilizando los scripts que ya existen ----------

def correr(script, *argumentos):
    """Corre un script del repo y deja ver su salida. Si falla, se detiene aqui."""
    print(f'\n>>> {script} ' + ' '.join(map(str, argumentos)), flush=True)
    orden = [sys.executable, '-X', 'utf8', str(RAIZ / 'scripts' / script), *map(str, argumentos)]
    resultado = subprocess.run(orden, cwd=RAIZ)
    if resultado.returncode:
        sys.exit(f'{script} termino con error; se detiene el recorrido.')


def sha256(ruta):
    return hashlib.sha256(Path(ruta).read_bytes()).hexdigest()


def huella(ifc, paso, horizonte=None):
    """Resume todas las entradas de una corrida: si cambia una, se repite la corrida."""
    datos = dict(ifc=sha256(ifc), paso=list(paso))
    if horizonte is not None:
        datos['horizonte_h'] = horizonte
        datos['configuracion'] = {c: sha256(RAIZ / 'config' / c) for c in CONFIGURACION}
    return datos


def vigente(carpeta, datos, archivo):
    """True si la carpeta ya tiene el resultado de exactamente estas entradas."""
    ruta = carpeta / 'entradas.json'
    return (carpeta / archivo).is_file() and ruta.is_file() and json.loads(ruta.read_text(encoding='utf-8')) == datos


def anotar(carpeta, datos):
    (carpeta / 'entradas.json').write_text(json.dumps(datos, indent=2), encoding='utf-8')


def voxelizar(ifc, carpeta, paso, figuras=True):
    datos = huella(ifc, paso)
    if vigente(carpeta, datos, 'grilla_vertidas.npz'):
        print(f'Voxelizacion ya calculada, se reutiliza: {carpeta}', flush=True)
        return
    extras = [] if figuras else ['--sin-figuras']
    correr('07_voxelizar_ifc.py', '--ifc', ifc, '--dx', paso[0], '--dy', paso[1], '--dz', paso[2], '--salida', carpeta, *extras)
    anotar(carpeta, datos)


def simular(ifc, voxeles, carpeta, paso, horizonte):
    datos = huella(ifc, paso, horizonte)
    if vigente(carpeta, datos, 'resultado_termico.npz'):
        print(f'Simulacion ya calculada, se reutiliza: {carpeta}', flush=True)
        return
    correr('09_simular_fundacion_vertidas.py', '--grilla', voxeles / 'grilla_vertidas.npz',
           '--horizonte-h', horizonte, '--salida', carpeta)
    anotar(carpeta, datos)


def crear_ifc_modificado(original, destino, inicios, temperaturas):
    """Escribe una COPIA del IFC con las propiedades cambiadas. El original no se toca."""
    modelo = bim.abrir(original)
    for vertida, valor in inicios:
        bim.editar_propiedad(modelo, vertida, 'Inicio_h', valor)
    for vertida, valor in temperaturas:
        bim.editar_propiedad(modelo, vertida, 'TemperaturaColocacion_C', valor)
    destino.parent.mkdir(parents=True, exist_ok=True)
    modelo.write(str(destino))


# ---------- Revisiones y tablas ----------

def programa_de(grilla):
    """Horario y temperatura que realmente recibe el solver: vienen del IFC via la grilla."""
    return {v: (float(p['Inicio_h']), float(p['TemperaturaColocacion_C'])) for v, p in sorted(grilla.vertidas.items())}


def avisos_programa(grilla, horizonte=None):
    """Avisa de horarios raros. No cambia el calculo: el solver no revisa el orden constructivo."""
    programa = programa_de(grilla)
    arriba, abajo = grilla.vertida_id[:, :, 1:], grilla.vertida_id[:, :, :-1]
    contacto = (arriba > 0) & (abajo > 0) & (arriba != abajo)
    pares = {(int(a), int(b)) for a, b in zip(arriba[contacto], abajo[contacto])}
    avisos = [f'V{a:02d} entra a las {programa[a][0]:g} h, antes que V{b:02d} que esta debajo ({programa[b][0]:g} h): quedaria sin apoyo.'
              for a, b in sorted(pares) if programa[a][0] < programa[b][0]]
    if horizonte is not None:
        avisos += [f'V{v:02d} entra a las {i:g} h, despues del horizonte de {horizonte:g} h: no alcanza a simularse.'
                   for v, (i, _) in programa.items() if i >= horizonte]
    return avisos


def mostrar_avisos(avisos):
    for aviso in avisos:
        print(f'AVISO: {aviso}', flush=True)


def imprimir_volumen(carpeta):
    """Muestra el error total de volumen que guardo el script 07."""
    r = json.loads((carpeta / 'resumen.json').read_text(encoding='utf-8'))
    print(f"\nError de volumen, con paso {r['paso_m']} m: IFC {r['volumen_ifc_m3']:.3f} m3, "
          f"voxeles {r['volumen_voxel_m3']:.3f} m3, error {r['error_relativo_porcentaje']:+.3f} %.", flush=True)


def leer_historias(carpeta):
    """Lee historias_por_vertida.csv como {vertida: {'t','media','max'}}; NaN si aun no existe."""
    filas = {}
    with (carpeta / 'historias_por_vertida.csv').open(encoding='utf-8-sig', newline='') as f:
        for fila in csv.DictReader(f):
            d = filas.setdefault(int(fila['vertida_id']), dict(t=[], media=[], max=[]))
            d['t'].append(float(fila['tiempo_h']))
            d['media'].append(float(fila['T_media_C']) if fila['T_media_C'] else np.nan)
            d['max'].append(float(fila['T_max_C']) if fila['T_max_C'] else np.nan)
    return {v: {k: np.array(x) for k, x in d.items()} for v, d in filas.items()}


def maximo_por_hora(historias):
    """Temperatura mas alta de toda la fundacion en cada instante guardado."""
    t = next(iter(historias.values()))['t']
    pila = np.array([d['max'] for d in historias.values()])
    return t, np.array([np.nan if np.isnan(c).all() else np.nanmax(c) for c in pila.T])


# ---------- Figuras ----------

def vista_3d(grilla, ocultas, ruta):
    """Dibuja los voxeles visibles como un solido. Ocultar solo cambia este dibujo.

    Una cara se dibuja si su celda es visible, la celda vecina no lo es y la cara
    mira hacia el observador. Es la misma idea de las caras expuestas del solver.
    """
    visibles = (grilla.vertida_id > 0) & ~np.isin(grilla.vertida_id, list(ocultas))
    con_borde = np.pad(visibles, 1)
    origen, paso = np.asarray(grilla.origen), np.asarray(grilla.paso)
    # Con esta camara se ven las caras x+, y- y z+; el brillo distinto marca la orientacion.
    caras = [(0, 1, .82), (1, -1, .66), (2, 1, 1.)]
    poligonos, colores = [], []
    for eje, signo, brillo in caras:
        vecina = np.roll(con_borde, -signo, axis=eje)[1:-1, 1:-1, 1:-1]
        indices = np.column_stack(np.nonzero(visibles & ~vecina))
        u, v = [e for e in range(3) if e != eje]
        esquina = origen + indices * paso
        esquina[:, eje] += paso[eje] * (signo > 0)
        du, dv = np.eye(3)[u] * paso[u], np.eye(3)[v] * paso[v]
        poligonos.append(np.stack([esquina, esquina + du, esquina + du + dv, esquina + dv], axis=1))
        base = np.array([matplotlib.colors.to_rgb(COLORES[(i - 1) % len(COLORES)]) for i in sorted(grilla.vertidas)])
        posicion = np.searchsorted(sorted(grilla.vertidas), grilla.vertida_id[tuple(indices.T)])
        colores.append(base[posicion] * brillo)
    fig = plt.figure(figsize=(10, 7), layout='constrained')
    ax = fig.add_subplot(projection='3d')
    ax.add_collection3d(Poly3DCollection(np.concatenate(poligonos), facecolors=np.concatenate(colores), edgecolors='none'))
    limites = [(origen[e], origen[e] + paso[e] * grilla.vertida_id.shape[e]) for e in range(3)]
    nota = 'Todas las vertidas visibles' if not ocultas else 'Oculta en el dibujo, sigue en el cálculo: ' + ', '.join(f'V{v:02d}' for v in sorted(ocultas))
    ax.set(xlabel='x [m]', ylabel='y [m]', zlabel='z [m]', xlim=limites[0], ylim=limites[1], zlim=limites[2], title=nota)
    ax.set_box_aspect([b - a for a, b in limites])
    ax.view_init(elev=24, azim=-48)
    ax.legend(handles=[Patch(color=COLORES[(v - 1) % len(COLORES)], label=f'V{v:02d}') for v in sorted(grilla.vertidas) if v not in ocultas],
              loc='upper left', fontsize=9)
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def figura_comparacion(base, modificado, cambiadas, horizonte, ruta):
    """Original contra modificado: curvas arriba y cortes termicos abajo."""
    h_base, h_mod = leer_historias(base), leer_historias(modificado)
    datos_base = np.load(base / 'resultado_termico.npz')
    datos_mod = np.load(modificado / 'resultado_termico.npz')
    horas = [h for h in (12, 24, 48, horizonte) if h <= horizonte]
    horas = list(dict.fromkeys(horas))

    def corte(datos, hora):
        return datos['cortes_xz_C'][int(np.argmin(abs(datos['tiempo_s'] / 3600 - hora)))]

    def extension(datos, campo):
        o, d = datos['origen_m'], datos['paso_m']
        return (o[0], o[0] + d[0] * campo.shape[1], o[2], o[2] + d[2] * campo.shape[0])

    cortes = [corte(d, h) for d in (datos_base, datos_mod) for h in horas]
    vmin, vmax = min(np.nanmin(c) for c in cortes), max(np.nanmax(c) for c in cortes)
    fig = plt.figure(figsize=(13, 7.5), layout='constrained')
    malla = fig.add_gridspec(3, len(horas), height_ratios=[2.6, 1, 1])
    arriba = malla[0, :].subgridspec(1, 2)
    izq, der = fig.add_subplot(arriba[0, 0]), fig.add_subplot(arriba[0, 1])
    for vid in cambiadas:
        color = COLORES[(vid - 1) % len(COLORES)]
        izq.plot(h_base[vid]['t'], h_base[vid]['media'], color=color, ls='--', lw=1.5, alpha=.7, label=f'V{vid:02d} original')
        izq.plot(h_mod[vid]['t'], h_mod[vid]['media'], color=color, lw=2, label=f'V{vid:02d} modificado')
    izq.set(title='Temperatura media de las vertidas cambiadas', xlabel='Tiempo [h]', ylabel='T [°C]')
    izq.legend(fontsize=9)
    der.plot(*maximo_por_hora(h_base), color='#6b7280', ls='--', lw=1.5, label='Original')
    der.plot(*maximo_por_hora(h_mod), color=COLORES[0], lw=2, label='Modificado')
    der.set(title='Temperatura máxima de toda la fundación', xlabel='Tiempo [h]', ylabel='T [°C]')
    der.legend(fontsize=9)
    for ax in (izq, der):
        ax.grid(alpha=.2)
    mapa = plt.get_cmap('inferno').copy()
    mapa.set_bad('#e5e7e9')
    ejes = []
    for fila, (datos, nombre) in enumerate(((datos_base, 'Original'), (datos_mod, 'Modificado')), 1):
        for col, hora in enumerate(horas):
            ax = fig.add_subplot(malla[fila, col])
            campo = corte(datos, hora)
            imagen = ax.imshow(np.ma.masked_invalid(campo), origin='lower', extent=extension(datos, campo), cmap=mapa,
                               vmin=vmin, vmax=vmax, interpolation='nearest')
            ax.set(title=f'{nombre}, t = {hora:g} h', xlabel='x [m]', ylabel='z [m]')
            ejes.append(ax)
    fig.colorbar(imagen, ax=ejes, label='Temperatura [°C]', shrink=.8)
    fig.suptitle('Efecto de cambiar solo el IFC · supuestos pedagógicos · gris = aún sin colocar o exterior')
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def imprimir_efecto(base, modificado, cambiadas, horizonte):
    r_base = json.loads((base / 'resumen.json').read_text(encoding='utf-8'))
    r_mod = json.loads((modificado / 'resumen.json').read_text(encoding='utf-8'))
    print(f"\nTmax global: original {r_base['maximo_global_C']:.3f} °C a {r_base['tiempo_maximo_h']:.1f} h; "
          f"modificado {r_mod['maximo_global_C']:.3f} °C a {r_mod['tiempo_maximo_h']:.1f} h.", flush=True)
    h_base, h_mod = leer_historias(base), leer_historias(modificado)
    print('Vertida  hora [h]  T media original [°C]  T media modificado [°C]', flush=True)
    for vid in cambiadas:
        for hora in [h for h in (12, 24, 48, horizonte) if h <= horizonte]:
            a = np.interp(hora, h_base[vid]['t'], h_base[vid]['media'])
            b = np.interp(hora, h_mod[vid]['t'], h_mod[vid]['media'])
            texto_a = 'sin colocar' if np.isnan(a) else f'{a:.2f}'
            texto_b = 'sin colocar' if np.isnan(b) else f'{b:.2f}'
            print(f'V{vid:02d}     {hora:8g}  {texto_a:>21}  {texto_b:>23}', flush=True)


# ---------- Recorrido ----------

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ifc', type=Path, default=RAIZ / 'fundaciones' / 'fundacion.ifc', help='IFC original; no se modifica')
    p.add_argument('--dx', type=float, default=.2)
    p.add_argument('--dy', type=float, default=.2)
    p.add_argument('--dz', type=float, default=.05)
    p.add_argument('--inicio', type=leer_cambio, action='append', default=[], metavar='V02=24', help='Hora de inicio [h]; se puede repetir')
    p.add_argument('--temp', type=leer_cambio, action='append', default=[], metavar='V02=10', help='Temperatura de colocacion [°C]; se puede repetir')
    p.add_argument('--ifc-cambiado', type=Path, help='IFC ya editado, por ejemplo en Bonsai, para compararlo con --ifc')
    p.add_argument('--ocultar', type=leer_vertida, action='append', default=[], metavar='V05', help='Oculta una vertida solo en la vista 3D; se puede repetir')
    p.add_argument('--solo-voxelizar', action='store_true', help='Solo volumen y vista 3D, sin simular')
    p.add_argument('--horizonte-h', type=float, help='Por defecto, el de config/termica_e07.yaml')
    p.add_argument('--salida', type=Path, default=RAIZ / 'resultados' / 'lab')
    args = p.parse_args()

    paso = (args.dx, args.dy, args.dz)
    if min(paso) <= 0:
        p.error('--dx, --dy y --dz deben ser positivos.')
    if not args.ifc.is_file():
        p.error(f'No existe el IFC: {args.ifc}')
    if args.ifc_cambiado and not args.ifc_cambiado.is_file():
        p.error(f'No existe el IFC cambiado: {args.ifc_cambiado}')
    if args.ifc_cambiado and (args.inicio or args.temp):
        p.error('Use --ifc-cambiado o --inicio/--temp, no ambos.')
    if any(valor < 0 for _, valor in args.inicio):
        p.error('La hora de inicio no puede ser negativa.')
    hay_cambios = bool(args.inicio or args.temp or args.ifc_cambiado)
    if args.solo_voxelizar and hay_cambios:
        p.error('--solo-voxelizar no simula, asi que no puede mostrar un cambio de hora o temperatura.')
    horizonte = args.horizonte_h if args.horizonte_h is not None else float(bim.cargar_yaml(RAIZ / 'config' / 'termica_e07.yaml')['horizonte_h'])
    if not math.isfinite(horizonte) or horizonte <= 0:
        p.error('--horizonte-h debe ser positivo.')

    existentes = {bim._vertida_id_numero(v['VertidaID']) for v in bim.leer_vertidas(bim.abrir(args.ifc))}
    pedidas = {v for v, _ in args.inicio + args.temp} | set(args.ocultar)
    if pedidas - existentes:
        p.error('No existen en el IFC: ' + ', '.join(f'V{v:02d}' for v in sorted(pedidas - existentes)))
    if set(args.ocultar) >= existentes:
        p.error('No se pueden ocultar todas las vertidas.')

    huella_original = sha256(args.ifc)
    malla = args.salida / f'dx{args.dx:g}_dy{args.dy:g}_dz{args.dz:g}'
    vox_base = malla / 'base' / 'voxeles'
    voxelizar(args.ifc, vox_base, paso)
    imprimir_volumen(vox_base)
    grilla_base = GrillaVoxeles.cargar(vox_base / 'grilla_vertidas.npz')

    sufijo = 'todas' if not args.ocultar else 'sin_' + '_'.join(f'V{v:02d}' for v in sorted(set(args.ocultar)))
    ruta_3d = malla / f'vista_3d_{sufijo}.png'
    vista_3d(grilla_base, set(args.ocultar), ruta_3d)
    print(f'Vista 3D: {ruta_3d}', flush=True)
    if args.solo_voxelizar:
        mostrar_avisos(avisos_programa(grilla_base))
        return

    sim_base = malla / 'base' / f'simulacion_{horizonte:g}h'
    simular(args.ifc, vox_base, sim_base, paso, horizonte)
    if not hay_cambios:
        mostrar_avisos(avisos_programa(grilla_base, horizonte))
        print(f'\nResultados: {sim_base}', flush=True)
        return

    carpeta = malla / nombre_escenario(args)
    if args.ifc_cambiado:
        ifc_mod = args.ifc_cambiado
    else:
        ifc_mod = carpeta / 'fundacion_modificada.ifc'
        crear_ifc_modificado(args.ifc, ifc_mod, args.inicio, args.temp)
    vox_mod = carpeta / 'voxeles'
    voxelizar(ifc_mod, vox_mod, paso, figuras=False)
    grilla_mod = GrillaVoxeles.cargar(vox_mod / 'grilla_vertidas.npz')

    # Lo que sigue se lee de la grilla, que es lo que recibe el solver. Asi se ve que el dato viene del IFC.
    antes, despues = programa_de(grilla_base), programa_de(grilla_mod)
    print('\nPrograma leido del IFC (original -> modificado):', flush=True)
    for vid in sorted(despues):
        ia, ta = antes.get(vid, (float('nan'), float('nan')))
        ib, tb = despues[vid]
        marca = '  <- cambio' if (ia, ta) != (ib, tb) else ''
        print(f'V{vid:02d}: inicio {ia:g} -> {ib:g} h; colocacion {ta:g} -> {tb:g} °C{marca}', flush=True)
    if antes == despues:
        print('AVISO: el IFC cambiado tiene el mismo horario y temperaturas que el original. Revise que haya guardado la edicion.', flush=True)
    mostrar_avisos(avisos_programa(grilla_mod, horizonte))

    sim_mod = carpeta / f'simulacion_{horizonte:g}h'
    simular(ifc_mod, vox_mod, sim_mod, paso, horizonte)
    cambiadas = [v for v in sorted(despues) if antes.get(v) != despues[v] and v in antes] or sorted(despues)
    ruta_comparacion = carpeta / f'comparacion_{horizonte:g}h.png'
    figura_comparacion(sim_base, sim_mod, cambiadas, horizonte, ruta_comparacion)
    imprimir_efecto(sim_base, sim_mod, cambiadas, horizonte)
    if sha256(args.ifc) != huella_original:
        sys.exit('El IFC original cambio durante el recorrido; revise que nada mas lo este editando.')
    print(f'\nIFC original sin cambios. Comparacion: {ruta_comparacion}', flush=True)


if __name__ == '__main__':
    main()
