"""E11: abre un campo 3D exportado y muestra clip, slice, nucleo y superficie (PyVista).

python -X utf8 scripts/15_ver_vtk.py resultados/e11/clima_real_24h/campo_24h.vtk
python -X utf8 scripts/15_ver_vtk.py CAMPO.vtk --plano x --pos 2.5
python -X utf8 scripts/15_ver_vtk.py CAMPO.vtk --png figura.png     (sin ventana, guarda imagen)

Requiere PyVista: pip install -r requirements-visual.txt
Izquierda: CLIP, se corta el solido con un plano y se ve el interior.
Derecha: SLICE, solo el plano de corte, con isotermas cada 5 C.
Nucleo: la celda mas caliente de la peor vertida. Superficie: su celda mas fria con una
cara al aire (la misma definicion del visor HTML y de mcplan.suelo.diferencial_por_vertida).
En pantalla se marcan los de la cara cortada; la consola informa tambien los de todo el volumen.
El mismo VTK se abre en ParaView: Threshold (activa = 1) -> Clip o Slice, color temperatura_C.
"""
from pathlib import Path
import argparse
import numpy as np

EJES = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}


def nucleo_y_superficie(malla, centros, eje, pos, dentro=None):
    """Indices del nucleo y la superficie de la vertida con mayor diferencial.

    dentro: mascara opcional de celdas candidatas (p. ej. las que corta el plano).
    Entre celdas empatadas (la fundacion es simetrica) se elige la del lado que queda a la vista.
    """
    T = malla.cell_data['temperatura_C']
    vid = malla.cell_data['vertida_id']
    sup = malla.cell_data['superficie'] > 0
    lado = centros[:, eje] <= pos + 1e-9
    base = np.ones(len(T), bool) if dentro is None else dentro
    def elegir(m, mejor):
        cand = np.flatnonzero(m & (np.abs(T - mejor) < 1e-6))
        visibles = cand[lado[cand]]
        return int((visibles if len(visibles) else cand)[0])
    peor = None
    for v in np.unique(vid):
        if v == 0 or v >= 90:
            continue
        m = (vid == v) & base
        s = m & sup
        if not s.any():
            continue
        i_n = elegir(m, T[m].max())
        i_s = elegir(s, T[s].min())
        d = float(T[i_n] - T[i_s])
        if peor is None or d > peor[0]:
            peor = (d, int(v), i_n, i_s)
    return peor


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('vtk', type=Path)
    p.add_argument('--plano', choices=list(EJES), default='y', help='Normal del plano de corte')
    p.add_argument('--pos', type=float, default=0.0, help='Posicion del plano [m]')
    p.add_argument('--sin-suelo', action='store_true', help='Oculta el bloque de suelo')
    p.add_argument('--png', type=Path, help='Guarda una imagen en vez de abrir la ventana')
    args = p.parse_args()
    try:
        import pyvista as pv
    except ImportError:
        p.error('Falta PyVista: pip install -r requirements-visual.txt (o abra el VTK en ParaView).')
    if not args.vtk.is_file():
        p.error(f'No existe {args.vtk}. Genere uno con scripts/14_simular_ambiente.py --vtk 24.')

    grilla = pv.read(str(args.vtk))
    activa = grilla.threshold(value=(0.5, 1.5), scalars='activa', preference='cell')
    if args.sin_suelo:
        activa = activa.threshold(value=(0.5, 89.5), scalars='vertida_id', preference='cell')
    centros = activa.cell_centers().points
    eje = 'xyz'.index(args.plano)
    paso = grilla.spacing[eje]
    peor = nucleo_y_superficie(activa, centros, eje, args.pos)
    en_plano = np.abs(centros[:, eje] - args.pos) <= paso / 2 + 1e-4   # puntos en float32
    peor_plano = nucleo_y_superficie(activa, centros, eje, args.pos, dentro=en_plano)
    T = activa.cell_data['temperatura_C']
    print(f'Archivo: {args.vtk} | celdas activas: {activa.n_cells} | T de {np.nanmin(T):.2f} a {np.nanmax(T):.2f} C')
    if peor:
        d, v, i_n, i_s = peor
        print(f'Peor vertida: V{v:02d}. Nucleo {T[i_n]:.2f} C en {np.round(centros[i_n], 2)}; '
              f'superficie {T[i_s]:.2f} C en {np.round(centros[i_s], 2)}; diferencia {d:.2f} C.')
    if peor_plano:
        d, v, i_n, i_s = peor_plano
        print(f'En el plano {args.plano} = {args.pos:g} m: V{v:02d}, nucleo {T[i_n]:.2f} C y superficie {T[i_s]:.2f} C ({d:.2f} C).')

    normal = EJES[args.plano]
    origen = [args.pos if c else 0.0 for c in normal]
    recortado = activa.clip(normal=normal, origin=origen, invert=True)   # queda el lado de atras: se ve la cara cortada
    corte = activa.slice(normal=normal, origin=origen)
    rango = [float(np.nanmin(T)), float(np.nanmax(T))]
    pl = pv.Plotter(shape=(1, 2), off_screen=bool(args.png), window_size=(1500, 650))
    barra = dict(title='Temperatura [°C]', n_labels=5, fmt='%.0f')
    for col, (malla, nombre) in enumerate(((recortado, f'Clip en {args.plano} = {args.pos:g} m'),
                                           (corte, f'Slice en {args.plano} = {args.pos:g} m'))):
        pl.subplot(0, col)
        if malla.n_cells:
            pl.add_mesh(malla, scalars='temperatura_C', cmap='inferno', clim=rango, scalar_bar_args=barra, show_edges=False)
        if col == 1 and corte.n_cells:
            iso = corte.cell_data_to_point_data().contour(np.arange(np.ceil(rango[0] / 5) * 5, rango[1], 5), scalars='temperatura_C')
            if iso.n_points:
                pl.add_mesh(iso, color='white', line_width=1.5)
        marcas = peor_plano   # los puntos que se ven en la cara cortada y en el plano
        if marcas:
            for i, color, etiqueta in ((marcas[2], '#18b5d6', 'Núcleo'), (marcas[3], '#7cc242', 'Superficie')):
                pl.add_mesh(pv.Sphere(radius=0.3, center=centros[i]), color=color)
                pl.add_point_labels([centros[i]], [f'{etiqueta} {T[i]:.1f} °C'], font_size=12, point_size=1,
                                    shape_opacity=0.7, always_visible=True)
        pl.add_text(nombre, font_size=11)
        pl.show_axes()
        if col == 0:
            mira = np.array([1.0, 1.0, 0.8]); mira[eje] = 1.6   # mirar desde el lado del corte
            pl.camera.position = tuple(np.array(activa.center) + 42 * mira / np.linalg.norm(mira))
            pl.camera.focal_point = activa.center; pl.camera.up = (0, 0, 1)
        else:
            pl.camera_position = {'x': 'yz', 'y': 'xz', 'z': 'xy'}[args.plano]
    if args.png:
        pl.screenshot(str(args.png))
        print(f'Imagen: {args.png}')
    else:
        pl.show()


if __name__ == '__main__':
    main()
