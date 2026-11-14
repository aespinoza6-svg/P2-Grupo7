"""E10: simulacion por vertidas con clima real, suelo y moldaje/aislacion.

python -X utf8 scripts/14_simular_ambiente.py
python -X utf8 scripts/14_simular_ambiente.py --constante --salida resultados/e10/constante
python -X utf8 scripts/14_simular_ambiente.py --comparar-con resultados/e10/constante
python -X utf8 scripts/14_simular_ambiente.py --horizonte-h 24 --vtk 12,24 --salida resultados/e11/clima_real_24h
python -X utf8 scripts/14_simular_ambiente.py --vtk-maximos --salida resultados/e12/caso_14dias   (E12)

Entradas: grilla del script 07, config/ambiente_e10.yaml y un snapshot de data/weather/.
Todo lo que cambia el resultado queda en resumen.json (snapshot, fecha de descarga,
metodo de relleno, huellas SHA-256 y commit), para que cada figura sea trazable.
--constante reproduce el ambiente de E07 (aire 15 C, h = 10, sin sol, sin suelo ni moldaje).
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mcplan import bim
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import modelos_desde_configuracion
from mcplan.construccion import ProgramaVertidas
from mcplan.solver1d import Material
from mcplan.solver_vertidas import resolver_vertidas
from mcplan.clima import cargar_snapshot, ultimo_snapshot, serie_clima, serie_constante
from mcplan.bordes import caras_ambientales
from mcplan.suelo import SUELO_ID, agregar_suelo, modelo_suelo, nucleo_superficie_por_vertida
from mcplan.exportar import campo_vtk

RAIZ = Path(__file__).resolve().parents[1]
COLORES = ['#287bb5', '#d97922', '#258b62', '#985bad', '#cd4c44']


def celda(g, i):
    """Centro [m] e indices (i, j, k) de la celda de indice plano i. z = 0 es la base de la fundacion."""
    ijk = np.unravel_index(i, g.vertida_id.shape)
    xyz = np.asarray(g.origen, dtype=float) + (np.array(ijk) + 0.5) * np.asarray(g.paso, dtype=float)
    return dict(x_m=round(float(xyz[0]), 3), y_m=round(float(xyz[1]), 3), z_m=round(float(xyz[2]), 3),
                indice=[int(n) for n in ijk], vertida=int(g.vertida_id[ijk]))


def sha256(ruta):
    return hashlib.sha256(Path(ruta).read_bytes()).hexdigest()


def commit():
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    except OSError:
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=RAIZ / 'config' / 'ambiente_e10.yaml')
    p.add_argument('--grilla', type=Path, default=RAIZ / 'resultados' / 'e07_voxelizacion' / 'grilla_vertidas.npz')
    p.add_argument('--snapshot', help='Carpeta en data/weather/; por defecto la de la configuracion')
    p.add_argument('--inicio', help='Fecha ISO de la primera vertida; por defecto la de la configuracion')
    p.add_argument('--horizonte-h', type=float)
    p.add_argument('--constante', action='store_true', help='Ambiente de E07 para comparar')
    p.add_argument('--sin-suelo', action='store_true')
    p.add_argument('--comparar-con', type=Path, help='Otra salida de este script para superponer curvas')
    p.add_argument('--vtk', help='Horas a exportar como campo 3D, p. ej. 24,72 (ParaView o scripts/15_ver_vtk.py)')
    p.add_argument('--vtk-maximos', action='store_true', help='Exporta el campo 3D en el instante de T max y de ΔT max')
    p.add_argument('--salida', type=Path, default=RAIZ / 'resultados' / 'e10' / 'clima_real')
    args = p.parse_args()
    if not args.grilla.is_file():
        p.error('Falta la grilla: corra scripts/07_voxelizar_ifc.py.')
    amb = bim.cargar_yaml(args.config)
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    hc = bim.cargar_yaml(RAIZ / 'config' / 'hidratacion_e07.yaml')
    inicio = args.inicio or amb['inicio_primera_vertida']
    horizonte = float(args.horizonte_h or amb['horizonte_h'])

    # --- Geometria: hormigon del IFC + suelo
    g = GrillaVoxeles.cargar(args.grilla)
    s = amb['suelo']
    con_suelo = s['incluir'] and not args.sin_suelo and not args.constante
    if con_suelo:
        g = agregar_suelo(g, espesor_m=s['espesor_m'], margen_m=s['margen_m'], T_inicial_C=s['T_inicial_C'])
    modelos = modelos_desde_configuracion(caso, hc)
    h = caso['hormigon']
    materiales = {m: Material(m, h['conductividad_W_mK'], h['densidad_kg_m3'], h['calor_especifico_J_kgK']) for m in modelos}
    if con_suelo:
        modelos['SUELO'] = modelo_suelo(s['densidad_kg_m3'], s['calor_especifico_J_kgK'])
        materiales['SUELO'] = Material('SUELO', s['k_W_mK'], s['densidad_kg_m3'], s['calor_especifico_J_kgK'])
    programa = ProgramaVertidas(g, modelos)

    # --- Clima: snapshot (pronostico + cola historica) o constante
    if args.constante:
        clima = serie_constante(15.0, horizonte)
        snap = None
    else:
        carpeta_base = RAIZ / amb['clima']['carpeta']
        nombre = args.snapshot or amb['clima']['snapshot']
        carpeta = ultimo_snapshot(carpeta_base) if nombre == 'ultimo' else carpeta_base / nombre
        snap = cargar_snapshot(carpeta)
        clima = serie_clima(snap, inicio, horizonte, amb['clima']['metodo_cola'])

    # --- Bordes: moldaje en lados, aislacion arriba, suelo abajo
    sup, mol, ais = amb['superficie'], amb['moldaje'], amb['aislacion_superior']
    hormigon = [v for v in programa.vertidas if v != SUELO_ID]
    if args.constante:
        bordes = caras_ambientales(clima, h_fijo=10.0, absortividad=0.0)
        moldaje = aislacion = {}
    else:
        moldaje = {v: (float(mol['R_m2K_W']), float(g.vertidas[v]['RetiroMoldaje_h']) if mol['retiro_desde_ifc'] else None)
                   for v in hormigon} if mol['R_m2K_W'] else {}
        aislacion = {v: (float(ais['R_m2K_W']), ais['retiro_h']) for v in hormigon} if ais['R_m2K_W'] else {}
        bordes = caras_ambientales(clima, moldaje=moldaje, aislacion=aislacion,
                                   suelo=dict(id=SUELO_ID, T_fondo_C=s['T_fondo_C']) if con_suelo else None,
                                   absortividad=sup['absortividad_solar'], emisividad=sup['emisividad'],
                                   h_fijo=sup['h_fijo_W_m2K'])

    # --- Simulacion
    salidas = np.unique(np.r_[np.arange(0, horizonte, amb['intervalo_salida_h']), horizonte]) * 3600
    dif, tmax_h = [], []
    horas_vtk = sorted({float(x) for x in args.vtk.split(',')}) if args.vtk else []
    if any(not 0 <= x <= horizonte for x in horas_vtk):
        p.error('--vtk: las horas deben estar dentro del horizonte.')
    campos = {}
    ns_hist = []                         # nucleo y superficie de cada vertida en cada salida
    maximos = dict(T=(-np.inf, None, None), dT=(-np.inf, None, None))   # (valor, ubicacion, campo)
    def al_guardar(t, T, activa, sistema):
        th = t / 3600
        for hv in horas_vtk:
            if abs(th - hv) < 1e-6:
                campos[hv] = (T.copy(), activa.copy())
        ns = nucleo_superficie_por_vertida(T, activa, sistema)
        ns_hist.append((th, ns))
        d, vid = np.nan, None
        for v, (tn, ts, _, _) in ns.items():
            if not tn - ts <= d:
                d, vid = tn - ts, v
        hormigon_activo = activa & (g.vertida_id != SUELO_ID)
        dif.append((th, d, vid))
        tmax_h.append(float(T[hormigon_activo].max()) if hormigon_activo.any() else np.nan)
        # Donde ocurren los maximos: se guarda la celda (y el campo si se pide --vtk-maximos)
        if hormigon_activo.any() and tmax_h[-1] > maximos['T'][0]:
            i = int(np.flatnonzero(hormigon_activo)[np.argmax(T[hormigon_activo])])
            maximos['T'] = (tmax_h[-1], dict(t_h=th, celda=celda(g, i)), (T.copy(), activa.copy()) if args.vtk_maximos else None)
        if vid is not None and d > maximos['dT'][0]:
            tn, ts, i_n, i_s = ns[vid]
            maximos['dT'] = (d, dict(t_h=th, vertida=vid, T_nucleo_C=round(tn, 3), T_superficie_C=round(ts, 3),
                                     nucleo=celda(g, i_n), superficie=celda(g, i_s)),
                             (T.copy(), activa.copy()) if args.vtk_maximos else None)
    reloj = time.monotonic()
    print(f'Grilla {g.vertida_id.shape}, {np.count_nonzero(g.vertida_id)} celdas ({"con" if con_suelo else "sin"} suelo). '
          f'Clima: {clima.snapshot}, {clima.horas_pronostico} h de pronostico y {clima.horas_relleno} h de relleno ({clima.metodo_cola}).', flush=True)
    ultimo = [0.0]
    def progreso(msg):
        if time.monotonic() - ultimo[0] > 20 or msg.startswith('Colocacion'):
            print(msg, flush=True); ultimo[0] = time.monotonic()
    r = resolver_vertidas(programa, materiales, bordes, horizonte * 3600, max_dt_s=float(amb['max_dt_s']),
                          tiempos_guardado_s=salidas, progreso=progreso, al_guardar=al_guardar)
    escala = np.array([max(1., abs(f['energia_colocacion_J']) + abs(f['energia_hidratacion_J']) + abs(f['energia_bordes_J'])) for f in r.balance])
    relativo = float(np.max(np.abs([f['residuo_balance_J'] for f in r.balance]) / escala))
    assert relativo < 1e-8, 'No cierra el balance de energia.'

    # --- Salidas trazables
    out = args.salida; out.mkdir(parents=True, exist_ok=True)
    horas = r.t / 3600
    filas = [dict(tiempo_h=hh, T_max_hormigon_C=tm, dT_nucleo_superficie_C=d, vertida_peor=v,
                  T_aire_C=clima.en(hh * 3600)[0], viento_m_s=clima.en(hh * 3600)[1], ghi_W_m2=clima.en(hh * 3600)[2],
                  fuente_clima=clima.fuente[min(int(round(hh)), len(clima.fuente) - 1)])
             for (hh, d, v), tm in zip(dif, tmax_h)]
    with (out / 'historia.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)
    with (out / 'historias_por_vertida.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(r.por_vertida[0])); w.writeheader(); w.writerows(r.por_vertida)
    dT = np.array([f['dT_nucleo_superficie_C'] for f in filas], dtype=float)
    tm = np.array(tmax_h)
    i_dt, i_tm = int(np.nanargmax(dT)), int(np.nanargmax(tm))
    resumen = dict(
        alcance='Simulacion determinista con clima real, suelo y moldaje; hidratacion pedagogica (45 C, 30 h) sin calibrar.',
        inicio_primera_vertida=inicio, horizonte_h=horizonte, ambiente='constante_E07' if args.constante else 'clima_real',
        clima=dict(snapshot=clima.snapshot, retrieved_at=snap['metadata']['retrieved_at'] if snap else None,
                   source=snap['metadata']['source'] if snap else None, model=snap['metadata']['model'] if snap else None,
                   metodo_cola=clima.metodo_cola, horas_pronostico=clima.horas_pronostico, horas_relleno=clima.horas_relleno,
                   sha256_forecast=sha256(Path(snap['carpeta']) / 'forecast.json') if snap else None),
        suelo=s if con_suelo else None, moldaje={str(k): v for k, v in moldaje.items()},
        aislacion_superior={str(k): v for k, v in aislacion.items()}, superficie=sup,
        T_max_hormigon_C=float(tm[i_tm]), tiempo_T_max_h=float(horas[i_tm]),
        dT_nucleo_superficie_max_C=float(dT[i_dt]), tiempo_dT_max_h=float(horas[i_dt]), vertida_dT_max=filas[i_dt]['vertida_peor'],
        ubicacion_T_max=maximos['T'][1], ubicacion_dT_max=maximos['dT'][1],
        cumple=dict(T_max=bool(tm[i_tm] <= 70), dT=bool(dT[i_dt] <= 20)),
        criterios=dict(T_max_C=70, dT_C=20), residuo_balance_relativo_max=relativo, pasos=r.n_pasos,
        energia_bordes_J=r.balance[-1]['energia_bordes_J'], energia_hidratacion_J=r.balance[-1]['energia_hidratacion_J'],
        trazabilidad=dict(commit=commit(), sha256_config=sha256(args.config), sha256_grilla=sha256(args.grilla)),
        limitaciones=['Cielo a la temperatura del aire en la radiacion de onda larga', 'Caras verticales sin sol',
                      'h(V) por correlacion aproximada', 'Suelo con valores tipicos, sin ensayos',
                      'Diferencial nucleo-superficie aproximado por celdas, no a 25 mm',
                      'Hidratacion sin dependencia de la temperatura'])
    (out / 'resumen.json').write_text(json.dumps(resumen, ensure_ascii=False, indent=2), encoding='utf-8')
    otra = None
    if args.comparar_con:
        with (args.comparar_con / 'historia.csv').open(encoding='utf-8-sig') as f:
            otra = list(csv.DictReader(f))
        otra_res = json.loads((args.comparar_con / 'resumen.json').read_text(encoding='utf-8'))
    # Curvas nucleo / superficie de cada vertida (la tabla y su figura)
    hormigon_ids = sorted(v for v in g.vertidas if v != SUELO_ID)
    with (out / 'nucleo_superficie.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(['tiempo_h'] + [f'V{v:02d}_{p}_C' for v in hormigon_ids for p in ('nucleo', 'superficie')])
        for th, ns in ns_hist:
            w.writerow([th] + [x for v in hormigon_ids for x in (ns[v][:2] if v in ns else ('', ''))])
    graficar(out, horas, filas, clima, r, g, otra, otra_res['ambiente'] if otra else None)
    graficar_nucleo_superficie(out, ns_hist, hormigon_ids, clima, maximos['dT'][1])
    faltan = [hv for hv in horas_vtk if hv not in campos]
    if faltan:
        print(f'AVISO: no hay salida exacta en {faltan} h; use horas multiplo del intervalo de salida.')
    vtk = {}
    for hv, (T, activa) in campos.items():
        vtk[f'{hv:g}'] = campo_vtk(out / f'campo_{hv:g}h.vtk', g, T, activa, hv, f'snapshot {clima.snapshot}').name
    for clave, nombre in (('T', 'Tmax'), ('dT', 'dTmax')):
        if maximos[clave][2] is not None:
            T, activa = maximos[clave][2]; th = maximos[clave][1]['t_h']
            vtk[nombre] = campo_vtk(out / f'campo_{nombre}_{th:g}h.vtk', g, T, activa, th, f'{nombre} - snapshot {clima.snapshot}').name
    for k, ruta in vtk.items():
        print(f'Campo 3D ({k}): {out / ruta}')
    if vtk:
        resumen['vtk'] = vtk
        (out / 'resumen.json').write_text(json.dumps(resumen, ensure_ascii=False, indent=2), encoding='utf-8')
    for clave, u in (('T max', maximos['T'][1]), ('ΔT max', maximos['dT'][1])):
        print(f'{clave}: {json.dumps(u, ensure_ascii=False)}')
    print(f'T max hormigon {tm[i_tm]:.2f} C a las {horas[i_tm]:.1f} h; nucleo - superficie max {dT[i_dt]:.2f} C a las '
          f'{horas[i_dt]:.1f} h (V{filas[i_dt]["vertida_peor"]:02d}). Residuo {relativo:.1e}. {time.monotonic() - reloj:.0f} s.')
    print(f'Resultados: {out}')


def graficar(out, horas, filas, clima, r, g, otra, nombre_otra):
    fuente = np.array([f['fuente_clima'] for f in filas])
    def sombrear(ax):
        relleno = fuente == 'historico'
        if relleno.any():
            ax.fill_between(horas, 0, 1, where=relleno, transform=ax.get_xaxis_transform(), color='#9aa3a7', alpha=.18, lw=0,
                            label='Completado con histórico')
    fig, axs = plt.subplots(3, 1, figsize=(11, 8), sharex=True, layout='constrained')
    for ax, clave, et, col in zip(axs, ['T_aire_C', 'viento_m_s', 'ghi_W_m2'],
                                  ['Aire a 2 m [°C]', 'Viento a 10 m [m/s]', 'Radiación solar [W/m²]'], COLORES):
        ax.plot(horas, [f[clave] for f in filas], color=col, lw=1.4); ax.set_ylabel(et); ax.grid(alpha=.2); sombrear(ax)
    if (fuente == 'historico').any():
        axs[0].legend(loc='upper right', fontsize=8)
    axs[-1].set_xlabel('Tiempo desde la primera vertida [h]')
    fig.suptitle(f'Clima usado · snapshot {clima.snapshot} · cola: {clima.metodo_cola}')
    fig.savefig(out / 'clima.png', dpi=150); plt.close(fig)

    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    axs[0].plot(horas, [f['T_max_hormigon_C'] for f in filas], color='#b4531f', lw=2, label='Este caso')
    axs[1].plot(horas, [f['dT_nucleo_superficie_C'] for f in filas], color='#b4531f', lw=2, label='Este caso')
    if otra:
        th = [float(f['tiempo_h']) for f in otra]
        axs[0].plot(th, [float(f['T_max_hormigon_C'] or 'nan') for f in otra], color='#6b7280', ls='--', label=nombre_otra)
        axs[1].plot(th, [float(f['dT_nucleo_superficie_C'] or 'nan') for f in otra], color='#6b7280', ls='--', label=nombre_otra)
    axs[0].axhline(70, color='#a3281e', ls=':', lw=1.2, label='Límite 70 °C'); axs[1].axhline(20, color='#a3281e', ls=':', lw=1.2, label='Límite 20 °C')
    axs[0].set_ylabel('T máx del hormigón [°C]'); axs[1].set_ylabel('Núcleo − superficie [°C]')
    for ax in axs: ax.grid(alpha=.2); ax.legend(fontsize=8); sombrear(ax)
    axs[1].set_xlabel('Tiempo desde la primera vertida [h]')
    fig.suptitle('Temperatura máxima y diferencial por vertida (aproximación por celdas)')
    fig.savefig(out / 'temperaturas_y_diferencial.png', dpi=150); plt.close(fig)

    indices = [int(np.argmin(abs(horas - t))) for t in (24, 72, 168, horas[-1]) if t <= horas[-1]]
    fig, axs = plt.subplots(len(indices), 1, figsize=(11, 2.3 * len(indices)), layout='constrained', squeeze=False)
    ext = (g.origen[0], g.origen[0] + g.paso[0] * g.vertida_id.shape[0], g.origen[2], g.origen[2] + g.paso[2] * g.vertida_id.shape[2])
    mapa = plt.get_cmap('inferno').copy(); mapa.set_bad('#e5e7e9')
    vmax = max(np.nanmax(r.cortes_xz[j]) for j in indices); vmin = min(np.nanmin(r.cortes_xz[j]) for j in indices)
    for ax, j in zip(axs[:, 0], indices):
        im = ax.imshow(np.ma.masked_invalid(r.cortes_xz[j]), origin='lower', extent=ext, cmap=mapa, vmin=vmin, vmax=vmax, interpolation='nearest', aspect='auto')
        ax.axhline(0, color='white', lw=.6, ls='--'); ax.set(title=f't = {horas[j]:g} h', ylabel='z [m]')
    axs[-1, 0].set_xlabel('x [m]')
    fig.colorbar(im, ax=axs[:, 0].tolist(), label='Temperatura [°C]', shrink=.85)
    fig.suptitle(f'Cortes y = {r.y_corte_m:.2f} m · bajo la línea blanca, el suelo')
    fig.savefig(out / 'cortes_con_suelo.png', dpi=150); plt.close(fig)


def graficar_nucleo_superficie(out, ns_hist, ids, clima, u_dT):
    """Arriba: nucleo (linea llena) y superficie (punteada) de cada vertida, con el aire.
    Abajo: su diferencia, contra el criterio de 20 °C; el punto marca el maximo."""
    t = np.array([th for th, _ in ns_hist])
    fig, axs = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True, layout='constrained')
    axs[0].plot(t, [clima.en(x * 3600)[0] for x in t], color='#9aa3a7', lw=1, label='Aire')
    for v, col in zip(ids, COLORES):
        tn = np.array([ns[v][0] if v in ns else np.nan for _, ns in ns_hist])
        ts = np.array([ns[v][1] if v in ns else np.nan for _, ns in ns_hist])
        axs[0].plot(t, tn, color=col, lw=1.8, label=f'V{v:02d} núcleo')
        axs[0].plot(t, ts, color=col, lw=1.2, ls='--', label=f'V{v:02d} superficie')
        axs[1].plot(t, tn - ts, color=col, lw=1.6, label=f'V{v:02d}')
    axs[0].axhline(70, color='#a3281e', ls=':', lw=1.2, label='Límite 70 °C')
    axs[1].axhline(20, color='#a3281e', ls=':', lw=1.2, label='Límite 20 °C')
    if u_dT:
        axs[1].plot(u_dT['t_h'], u_dT['T_nucleo_C'] - u_dT['T_superficie_C'], 'o', color='#111', ms=6,
                    label=f"Máximo, V{u_dT['vertida']:02d} a las {u_dT['t_h']:g} h")
    axs[0].set_ylabel('Temperatura [°C]'); axs[1].set_ylabel('Núcleo − superficie [°C]')
    axs[1].set_xlabel('Tiempo desde la primera vertida [h]')
    axs[0].legend(fontsize=7, ncol=4); axs[1].legend(fontsize=8, ncol=3)
    for ax in axs: ax.grid(alpha=.2)
    fig.suptitle('Núcleo (celda más caliente) y superficie (celda expuesta más fría) de cada vertida')
    fig.savefig(out / 'nucleo_superficie.png', dpi=150); plt.close(fig)


if __name__ == '__main__':
    main()
