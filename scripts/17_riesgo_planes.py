"""E13: evalua cada plan de la familia en todos los escenarios (curva x clima) y su riesgo.

python -X utf8 scripts/17_riesgo_planes.py                       (malla de exploracion, todos los planes)
python -X utf8 scripts/17_riesgo_planes.py --malla final --planes cinco_e6_manta10,unica_manta10
python -X utf8 scripts/17_riesgo_planes.py --solo-resumir         (rehace tablas y figuras sin simular)

Cada corrida guarda una fila en escenarios.csv apenas termina; si se corta, al volver a
correr se saltan las que ya estan. Metricas por escenario s, como en el enunciado:
  MT = max en t y en el hormigon colocado de T;
  MD = max en t de [max T del hormigon colocado - min T a 25 mm bajo las caras exteriores].
Falla: MT > 70 o MD > 20. P_f = suma de los pesos w = p_m q_k de los escenarios que fallan.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import time
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
CAMPOS = ['plan', 'curva', 'clima', 'malla', 'D_h', 'horizonte_h', 'MT', 't_MT_h', 'xyz_MT', 'MD', 't_MD_h', 'xyz_nucleo_MD',
          'xyz_S25_MD', 'cara_S25', 'T_nucleo_MD', 'T_S25_MD', 'MD_celda', 'MT_336', 'MD_336', 'pico_al_final',
          'residuo', 'pasos', 'segundos']


# ---------------------------------------------------------------------------
# Lo que corre en cada proceso (se carga una vez por proceso)
# ---------------------------------------------------------------------------
_CACHE = {}


def _cargar(cfg, malla):
    if malla in _CACHE:
        return _CACHE[malla]
    from mcplan.voxelizacion import GrillaVoxeles
    from mcplan.hidratacion import leer_curvas_tabuladas
    from mcplan.clima import cargar_snapshot, snapshot_de_miembro, serie_clima
    tray = json.loads((RAIZ / 'resultados' / 'e13' / 'clima' / 'trayectorias.json').read_text(encoding='utf-8'))
    snap = cargar_snapshot(RAIZ / tray['carpeta'])
    H = tray['horizonte_h']
    climas = {e['nombre']: serie_clima(snapshot_de_miembro(snap, e['miembro']), cfg['t0'], H, cfg['clima']['metodo_cola'])
              for e in tray['trayectorias']}
    for clave in ('frio_extremo', 'calido_extremo'):
        climas[clave] = serie_clima(snapshot_de_miembro(snap, tray[clave]), cfg['t0'], H, cfg['clima']['metodo_cola'])
    _CACHE[malla] = (GrillaVoxeles.cargar(ruta_grilla(cfg, malla)), leer_curvas_tabuladas(RAIZ / cfg['curvas']['archivo']), climas)
    return _CACHE[malla]


def ruta_grilla(cfg, malla):
    m = cfg['mallas'][malla]
    if (m['dx'], m['dy'], m['dz']) == (0.2, 0.2, 0.05):
        return RAIZ / 'resultados' / 'e07_voxelizacion' / 'grilla_vertidas.npz'
    return RAIZ / 'resultados' / 'e13' / 'grillas' / f"grilla_{m['dx']:g}_{m['dy']:g}_{m['dz']:g}.npz"


def correr(cfg, plan, curva, clima, malla):
    """Una simulacion completa de un plan en un escenario. Devuelve una fila de CAMPOS."""
    from mcplan.construccion import ProgramaVertidas
    from mcplan.solver1d import Material
    from mcplan.solver_vertidas import resolver_vertidas
    from mcplan.bordes import caras_ambientales
    from mcplan.suelo import SUELO_ID, agregar_suelo, modelo_suelo
    from mcplan.riesgo import cronograma, grilla_con_plan, medir
    reloj = time.monotonic()
    grilla, curvas, climas = _cargar(cfg, malla)
    c = curvas[curva]
    vol = {v: float(grilla.vertidas[v]['volumen_ifc_m3']) for v in sorted(grilla.vertidas)}
    crono, D = cronograma(vol, float(cfg['colocacion']['capacidad_m3_h']), float(plan['espera_h']), plan['particion'] == 'unica')
    ultima = max(f['activacion_h'] for f in crono.values())
    H = max(float(cfg.get('horizonte_extendido_h', cfg['horizonte_h'])), ultima + float(cfg['extension_tras_ultima_h']))
    s = cfg['suelo']
    g = agregar_suelo(grilla_con_plan(grilla, crono, float(plan.get('temperatura_C', cfg['colocacion']['temperatura_C']))),
                      espesor_m=s['espesor_m'], margen_m=s['margen_m'], T_inicial_C=s['T_inicial_C'])
    k = float(cfg['hormigon']['conductividad_W_mK'])
    modelos = {'M1': c, 'SUELO': modelo_suelo(s['densidad_kg_m3'], s['calor_especifico_J_kgK'])}
    materiales = {'M1': Material('M1', k, c.rho_kg_m3, c.cp_J_kgK),
                  'SUELO': Material('SUELO', s['k_W_mK'], s['densidad_kg_m3'], s['calor_especifico_J_kgK'])}
    programa = ProgramaVertidas(g, modelos)
    prot = cfg['protecciones'][plan['proteccion']]
    (Rl, rl), (Ra, ra) = prot['lados'], prot['arriba']
    hormigon = [v for v in programa.vertidas if v != SUELO_ID]
    sup = cfg['superficie']
    aislacion = {v: (float(Ra), ra) for v in hormigon} if Ra else {}
    if 'terreno' in prot:                      # manta sobre el suelo que rodea la fundacion
        aislacion[SUELO_ID] = (float(prot['terreno'][0]), prot['terreno'][1])
    bordes = caras_ambientales(climas[clima], moldaje={v: (float(Rl), rl) for v in hormigon} if Rl else {},
                               aislacion=aislacion,
                               suelo=dict(id=SUELO_ID, T_fondo_C=s['T_fondo_C']),
                               absortividad=sup['absortividad_solar'], emisividad=sup['emisividad'])
    centro = lambda i: [round(float(x), 3) for x in np.asarray(g.origen) + (np.array(np.unravel_index(i, g.vertida_id.shape)) + .5) * np.asarray(g.paso)]
    hist = []
    def al_guardar(t, T, activa, sistema):
        m = medir(T, activa, sistema, t)
        if m:
            hist.append((t / 3600, m))
    salidas = np.unique(np.r_[np.arange(0, H, float(cfg['intervalo_salida_h'])), H]) * 3600
    r = resolver_vertidas(programa, materiales, bordes, H * 3600, max_dt_s=float(cfg['max_dt_s']),
                          tiempos_guardado_s=salidas, al_guardar=al_guardar)
    escala = np.array([max(1., abs(f['energia_colocacion_J']) + abs(f['energia_hidratacion_J']) + abs(f['energia_bordes_J'])) for f in r.balance])
    residuo = float(np.max(np.abs([f['residuo_balance_J'] for f in r.balance]) / escala))
    th = np.array([h for h, _ in hist]); mt = np.array([m['T_max'] for _, m in hist])
    md = np.array([m['T_max'] - m['T25_min'] for _, m in hist]); mdc = np.array([m['T_max'] - m['T_celda_min'] for _, m in hist])
    i_t, i_d = int(np.argmax(mt)), int(np.argmax(md)); dentro = th <= float(cfg['horizonte_h']) + 1e-9
    mT, mD = hist[i_t][1], hist[i_d][1]
    return dict(plan=plan['nombre'], curva=curva, clima=clima, malla=malla, D_h=round(D, 3), horizonte_h=H,
                MT=round(float(mt[i_t]), 3), t_MT_h=th[i_t], xyz_MT=centro(mT['i_max']), MD=round(float(md[i_d]), 3),
                t_MD_h=th[i_d], xyz_nucleo_MD=centro(mD['i_max']), xyz_S25_MD=centro(mD['i25']), cara_S25=mD['cara25'],
                T_nucleo_MD=round(mD['T_max'], 3), T_S25_MD=round(mD['T25_min'], 3), MD_celda=round(float(mdc.max()), 3),
                MT_336=round(float(mt[dentro].max()), 3), MD_336=round(float(md[dentro].max()), 3),
                pico_al_final=bool(th[i_t] > H - 24 or th[i_d] > H - 24), residuo=residuo, pasos=r.n_pasos,
                segundos=round(time.monotonic() - reloj, 1))


# ---------------------------------------------------------------------------
# Orquestacion
# ---------------------------------------------------------------------------
def preparar_grilla(cfg, malla):
    ruta = ruta_grilla(cfg, malla)
    if ruta.is_file():
        return ruta
    from mcplan.voxelizacion import extraer_mallas, voxelizar_mallas
    m = cfg['mallas'][malla]
    ruta.parent.mkdir(parents=True, exist_ok=True)
    voxelizar_mallas(extraer_mallas(RAIZ / 'fundaciones' / 'fundacion.ifc'), (m['dx'], m['dy'], m['dz']), max_celdas=12_000_000).guardar(ruta)
    return ruta


def leer_filas(ruta):
    if not ruta.is_file():
        return []
    with ruta.open(encoding='utf-8-sig') as f:
        filas = list(csv.DictReader(f))
    for x in filas:
        for c in ('D_h', 'horizonte_h', 'MT', 'MD', 'MD_celda', 'MT_336', 'MD_336', 't_MT_h', 't_MD_h'):
            x[c] = float(x[c])
    return filas


def resumir(cfg, out, filas, climas_peso):
    from mcplan.riesgo import probabilidades
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pesos = dict(base=cfg['curvas']['pesos']['base'], **cfg['curvas']['sensibilidad'])
    planes = [p['nombre'] for p in cfg['planes'] if any(f['plan'] == p['nombre'] for f in filas)]
    tabla = []
    for nombre in planes:
        F = [f for f in filas if f['plan'] == nombre]
        fila = dict(plan=nombre, D_h=F[0]['D_h'], escenarios=len(F), MT_max=max(f['MT'] for f in F), MD_max=max(f['MD'] for f in F),
                    MD_celda_max=max(f['MD_celda'] for f in F), pico_al_final=any(str(f['pico_al_final']) == 'True' for f in F))
        completo = len({(f['curva'], f['clima']) for f in F if f['clima'] in climas_peso}) == len(pesos['base']) * len(climas_peso)
        for clave, p in pesos.items():
            if completo:
                P = probabilidades([f for f in F if f['clima'] in climas_peso], p, climas_peso)
                fila[f'Pf_{clave}'] = round(P['conjunta'], 4)
                if clave == 'base':
                    fila.update(PT=round(P['T'], 4), PdT=round(P['dT'], 4), fallas=len(P['fallas']))
        fila['cumple'] = bool(completo and fila.get('Pf_base', 1) <= 0.05)
        # Sensibilidad climatica: los miembros extremos en vez de los representativos de sus grupos
        for extremo, reemplaza in (('frio_extremo', 'frio'), ('calido_extremo', 'calido')):
            G = [f for f in F if f['clima'] in climas_peso and f['clima'] != reemplaza] + \
                [dict(f, clima=reemplaza) for f in F if f['clima'] == extremo]
            if completo and len(G) == len(pesos['base']) * len(climas_peso):
                fila[f'Pf_con_{extremo}'] = round(probabilidades(G, pesos['base'], climas_peso)['conjunta'], 4)
        tabla.append(fila)
    tabla.sort(key=lambda x: x['D_h'])
    with (out / 'planes.csv').open('w', newline='', encoding='utf-8-sig') as f:
        claves = list(dict.fromkeys(k for t in tabla for k in t))
        w = csv.DictWriter(f, fieldnames=claves); w.writeheader(); w.writerows(tabla)
    factibles = [t for t in tabla if t['cumple']]
    resumen = dict(malla=out.name, criterio='P_f conjunta <= 0,05 con los pesos base', planes=tabla,
                   recomendado=factibles[0]['plan'] if factibles else None,
                   nota='Menor D entre los planes evaluados; no es una prueba de optimalidad global.')
    (out / 'resumen.json').write_text(json.dumps(resumen, ensure_ascii=False, indent=2), encoding='utf-8')

    # Figura 1: matriz de escenarios por plan (MD en cada celda; rojo = falla)
    curvas = list(pesos['base']); climas = list(climas_peso)
    fig, axs = plt.subplots(1, len(planes), figsize=(2.3 * len(planes) + 1, 4.2), layout='constrained', squeeze=False)
    for ax, nombre in zip(axs[0], planes):
        M = np.full((len(curvas), len(climas)), np.nan); falla = np.zeros_like(M, bool)
        for f in filas:
            if f['plan'] == nombre and f['clima'] in climas:
                i, j = curvas.index(f['curva']), climas.index(f['clima'])
                M[i, j] = f['MD']; falla[i, j] = f['MD'] > 20 or f['MT'] > 70
        ax.imshow(falla, cmap='RdYlGn_r', vmin=-.6, vmax=1.6, aspect='auto')
        for i in range(len(curvas)):
            for j in range(len(climas)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f'{M[i, j]:.1f}', ha='center', va='center', fontsize=7)
        t = next(x for x in tabla if x['plan'] == nombre)
        ax.set_title(f"{nombre}\nD = {t['D_h']:.1f} h · P_f = {100 * t.get('Pf_base', np.nan):.1f} %", fontsize=8)
        ax.set_xticks(range(len(climas)), climas, fontsize=7); ax.set_yticks(range(len(curvas)), curvas if ax is axs[0][0] else [], fontsize=7)
    fig.suptitle('MΔ [°C] por escenario (filas: curva; columnas: clima). Rojo: falla (MΔ > 20 °C o MT > 70 °C)')
    fig.savefig(out / 'matriz_escenarios.png', dpi=140); plt.close(fig)

    # Figura 2: D contra P_f
    fig, ax = plt.subplots(figsize=(8, 4.8), layout='constrained')
    for t in tabla:
        if 'Pf_base' in t:
            ax.plot(t['D_h'], 100 * t['Pf_base'], 'o', color='#258b62' if t['cumple'] else '#cd4c44')
            ax.annotate(t['plan'], (t['D_h'], 100 * t['Pf_base']), fontsize=7, xytext=(4, 4), textcoords='offset points')
    ax.axhline(5, color='#a3281e', ls=':', label='P_f = 5 %'); ax.set_xlabel('D = fin de la última vertida − t0 [h]')
    ax.set_ylabel('P_f conjunta [%]'); ax.grid(alpha=.3); ax.legend()
    fig.savefig(out / 'D_vs_Pf.png', dpi=140); plt.close(fig)
    return resumen


def figura_curvas(cfg, out):
    """Las 8 curvas juntas: Q(a) y su ΔT adiabatico acumulado, mas una tabla de resumen."""
    from mcplan.hidratacion import leer_curvas_tabuladas
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    curvas = leer_curvas_tabuladas(RAIZ / cfg['curvas']['archivo'])
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout='constrained')
    filas = []
    for n, c in curvas.items():
        h = c.a_s / 3600; ls = '-' if n.endswith('SS') else '--'
        ax[0].plot(h, c.q, ls, label=n); ax[1].plot(h, c.elevacion_adiabatica(c.a_s), ls, label=n)
        filas.append(dict(curva=n, pico_W_m3=round(float(c.q.max()), 1), t_pico_h=round(float(h[c.q.argmax()]), 2),
                          calor_total_MJ_m3=round(c.energia_total_J_m3 / 1e6, 2),
                          dT_adiabatico_K=round(c.energia_total_J_m3 / c.capacidad_volumetrica, 2), fin_registro_h=round(float(h[-1]), 1)))
    ax[0].set_ylabel('Q [W/m³ de hormigón]'); ax[1].set_ylabel('ΔT adiabático acumulado [K]'); ax[1].set_xlabel('Edad desde la activación [h]')
    ax[0].legend(ncol=4, fontsize=8); [a.grid(alpha=.3) for a in ax]
    fig.suptitle('Las 8 curvas del profesor (línea llena: campaña SS; punteada: 22A, menor confianza)')
    fig.savefig(out / 'curvas_8.png', dpi=140); plt.close(fig)
    with (out / 'curvas.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)


def main():
    from mcplan import bim
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=RAIZ / 'config' / 'riesgo_e13.yaml')
    p.add_argument('--malla', default='exploracion')
    p.add_argument('--planes', help='Nombres separados por coma; por defecto todos')
    p.add_argument('--curvas', help='Por defecto las 8')
    p.add_argument('--climas', help='Por defecto las 3 trayectorias y los 2 extremos')
    p.add_argument('--procesos', type=int, default=8)
    p.add_argument('--solo-resumir', action='store_true')
    p.add_argument('--suelo-T', help='Sensibilidad: T inicial y de fondo del suelo, p. ej. 12,12 (sale en otra carpeta)')
    args = p.parse_args()
    cfg = bim.cargar_yaml(args.config)
    sufijo = ''
    if args.suelo_T:
        ti, tf = (float(x) for x in args.suelo_T.split(','))
        cfg['suelo'].update(T_inicial_C=ti, T_fondo_C=tf); sufijo = f'_suelo{ti:g}_{tf:g}'
    tray = json.loads((RAIZ / 'resultados' / 'e13' / 'clima' / 'trayectorias.json').read_text(encoding='utf-8'))
    climas_peso = {e['nombre']: e['peso'] for e in tray['trayectorias']}
    out = RAIZ / 'resultados' / 'e13' / (args.malla + sufijo); out.mkdir(parents=True, exist_ok=True)
    ruta = out / 'escenarios.csv'
    if not args.solo_resumir:
        preparar_grilla(cfg, args.malla)
        planes = [x for x in cfg['planes'] if not args.planes or x['nombre'] in args.planes.split(',')]
        curvas = args.curvas.split(',') if args.curvas else list(cfg['curvas']['pesos']['base'])
        climas = args.climas.split(',') if args.climas else list(climas_peso) + ['frio_extremo', 'calido_extremo']
        hechas = {(f['plan'], f['curva'], f['clima']) for f in leer_filas(ruta)}
        tareas = [(pl, c, k) for pl in planes for c in curvas for k in climas if (pl['nombre'], c, k) not in hechas]
        print(f'{len(tareas)} corridas pendientes ({len(hechas)} ya hechas) con {args.procesos} procesos, malla {args.malla}.', flush=True)
        nuevo = not ruta.is_file()
        with ruta.open('a', newline='', encoding='utf-8-sig') as f, ProcessPoolExecutor(args.procesos) as pool:
            w = csv.DictWriter(f, fieldnames=CAMPOS)
            if nuevo:
                w.writeheader()
            futuros = {pool.submit(correr, cfg, pl, c, k, args.malla): (pl['nombre'], c, k) for pl, c, k in tareas}
            for n, fut in enumerate(as_completed(futuros), 1):
                fila = fut.result(); w.writerow(fila); f.flush()
                print(f"[{n}/{len(tareas)}] {fila['plan']} {fila['curva']} {fila['clima']}: MT {fila['MT']:.1f} MD {fila['MD']:.1f} "
                      f"(celda {fila['MD_celda']:.1f}) en {fila['segundos']:.0f} s", flush=True)
    figura_curvas(cfg, RAIZ / 'resultados' / 'e13' / 'curvas')
    filas = leer_filas(ruta)
    resumen = resumir(cfg, out, filas, climas_peso)
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    resumen.update(commit=commit, sha256_curvas=hashlib.sha256((RAIZ / cfg['curvas']['archivo']).read_bytes()).hexdigest(),
                   sha256_config=hashlib.sha256(args.config.read_bytes()).hexdigest(), clima=tray['snapshot'])
    (out / 'resumen.json').write_text(json.dumps(resumen, ensure_ascii=False, indent=2), encoding='utf-8')
    for t in resumen['planes']:
        print(f"{t['plan']:20s} D {t['D_h']:6.1f} h  MT {t['MT_max']:5.1f}  MΔ {t['MD_max']:5.1f}  "
              f"P_f {t.get('Pf_base', float('nan')):.3f}  {'CUMPLE' if t['cumple'] else ''}")
    print(f"Recomendado (menor D con P_f <= 5 %): {resumen['recomendado']}. Resultados: {out}")


if __name__ == '__main__':
    main()
