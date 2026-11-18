"""E13: descarga el ensemble de Open-Meteo y lo reduce a pocas historias climaticas con peso.

python -X utf8 scripts/16_clima_ensemble.py --descargar     (una vez: guarda data/weather_ensemble/<fecha>/)
python -X utf8 scripts/16_clima_ensemble.py                 (solo vuelve a elegir y graficar)

Cada miembro del ensemble es un futuro posible completo: aire, viento y sol coherentes
entre si y con su ciclo diario. Se ordenan por la temperatura media del aire en las horas
con pronostico, se parten en tres grupos iguales (frio, central, calido) y de cada grupo se
toma el miembro mas cercano a su mediana. Peso = fraccion de miembros que representa.
La cola despues del ultimo dia del ensemble es la misma para todos: promedio ERA5.
Grafica trayectorias.png y temperatura_min_max.png (minimas y maximas de los 21 dias).
"""
from datetime import datetime, timedelta
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mcplan import bim
from mcplan.clima import descargar_ensemble, cargar_snapshot, ultimo_snapshot, miembros_de, snapshot_de_miembro, serie_clima
from mcplan.riesgo import elegir_trayectorias

RAIZ = Path(__file__).resolve().parents[1]
COLORES = dict(frio='#287bb5', central='#258b62', calido='#cd4c44')
NOMBRE = dict(frio='fría', central='central', calido='cálida', frio_extremo='más frío', calido_extremo='más cálido')


def horizonte_maximo(cfg):
    return max(float(cfg.get('horizonte_extendido_h', 0)), float(cfg['horizonte_h']) + float(cfg['extension_tras_ultima_h']))


def series_del_ensemble(cfg, carpeta=None):
    """{miembro: SerieClima} desde t0 hasta el horizonte maximo de la familia de planes."""
    base = RAIZ / cfg['clima']['carpeta']
    carpeta = Path(carpeta) if carpeta else (ultimo_snapshot(base) if cfg['clima']['snapshot'] == 'ultimo' else base / cfg['clima']['snapshot'])
    snap = cargar_snapshot(carpeta)
    H = horizonte_maximo(cfg)
    return snap, {m: serie_clima(snapshot_de_miembro(snap, m), cfg['t0'], H, cfg['clima']['metodo_cola'])
                  for m in miembros_de(snap['pronostico'])}


def descargar_suelo(cfg, caso):
    """Temperatura del suelo de ERA5 en el sitio, ultimos dias disponibles antes de t0.

    El suelo no parte uniforme a una temperatura supuesta: su estado inicial sale de este
    archivo (capa 28-100 cm para el bloque y 100-255 cm para el fondo, ver config).
    """
    from datetime import timezone
    from mcplan.clima import _pedir, URL_HISTORICO
    sitio = caso['sitio']; ahora = datetime.now(timezone.utc)
    fin = ahora.date() - timedelta(days=4); inicio = fin - timedelta(days=11)
    capas = ['soil_temperature_0_to_7cm', 'soil_temperature_7_to_28cm', 'soil_temperature_28_to_100cm', 'soil_temperature_100_to_255cm']
    url, datos = _pedir(URL_HISTORICO, dict(latitude=sitio['latitud'], longitude=sitio['longitud'], hourly=','.join(capas),
                                            timezone=sitio['zona_horaria'], start_date=str(inicio), end_date=str(fin)))
    medias = {c: round(float(np.nanmean(np.array(datos['hourly'][c], dtype=float))), 2) for c in capas}
    carpeta = RAIZ / 'data' / 'suelo' / ahora.strftime('%Y-%m-%dT%H%M'); carpeta.mkdir(parents=True)
    meta = dict(source='Open-Meteo archive-api', model='ERA5 reanalysis (soil temperature)', retrieved_at=ahora.isoformat(timespec='seconds'),
                latitude=sitio['latitud'], longitude=sitio['longitud'], timezone=sitio['zona_horaria'], units='degC',
                periodo=[str(inicio), str(fin)], medias_C=medias, request=url)
    (carpeta / 'metadata.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    (carpeta / 'era5_suelo.json').write_text(json.dumps(datos, ensure_ascii=False), encoding='utf-8')
    print(f'Suelo ERA5 {inicio}..{fin}: ' + ', '.join(f"{c.replace('soil_temperature_', '')} {v} C" for c, v in medias.items()))
    print(f'Guardado en {carpeta}')


def figura_min_max(cfg, snap, series, elegidas, extra, ruta):
    """Temperatura del aire en los 21 dias que se simulan, con sus minimas y maximas.

    Arriba: hora a hora, el rango de los 51 miembros, las 3 trayectorias, los 2 miembros
    extremos y la minima y la maxima absolutas. Abajo: minima y maxima de cada dia
    (barra = rango del ensemble, punto = trayectoria central).
    """
    import matplotlib.dates as mdates
    H = float(cfg['horizonte_extendido_h'])
    t0 = datetime.fromisoformat(cfg['t0']).replace(tzinfo=None)          # hora de Chile
    s0 = next(iter(series.values()))
    dentro = s0.t_s / 3600 <= H + 1e-9
    fechas = np.array([t0 + timedelta(seconds=float(x)) for x in s0.t_s[dentro]])
    T = np.array([s.T_aire_C[dentro] for s in series.values()])          # miembro x hora
    fin = s0.t_s[s0.fuente == 'pronostico'].max() / 3600
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True, layout='constrained', height_ratios=[1.25, 1])

    a1.fill_between(fechas, T.min(0), T.max(0), color='#b8bec2', alpha=.45, lw=0, label=f'Rango de los {len(series)} miembros')
    for e in elegidas:
        a1.plot(fechas, series[e['miembro']].T_aire_C[dentro], color=COLORES.get(e['nombre'], 'k'), lw=1.3,
                label=f"Trayectoria {NOMBRE[e['nombre']]} ({e['miembro']})")
    for clave, color in (('frio_extremo', '#173f63'), ('calido_extremo', '#7a1f1a')):
        m = extra[clave]
        a1.plot(fechas, series[m].T_aire_C[dentro], color=color, lw=1, ls='--', label=f"Miembro {NOMBRE[clave]} ({m})")
    i_min, i_max = np.unravel_index(T.argmin(), T.shape), np.unravel_index(T.argmax(), T.shape)
    for (im, ih), txt, dy, c in ((i_min, 'Mínima', -5, '#173f63'), (i_max, 'Máxima', 3, '#7a1f1a')):
        v, f = T[im, ih], fechas[ih]
        a1.plot(f, v, 'o', ms=7, mfc='white', mec=c, mew=2, zorder=5)
        a1.annotate(f"{txt}: {v:.1f} °C\n{f:%d-%m %H:%M}".replace('.', ','), (f, v), xytext=(f, v + dy), ha='center',
                    fontsize=9, color=c, fontweight='bold', arrowprops=dict(arrowstyle='-', color=c, lw=.8))
    a1.set_ylabel('Temperatura del aire a 2 m [°C]')
    a1.set_ylim(T.min() - 7, T.max() + 6)
    a1.legend(fontsize=8, ncol=3, loc='lower left', bbox_to_anchor=(0, 1.01), frameon=False)

    dia = np.array([f.date() for f in fechas])
    dias = [d for d in sorted(set(dia)) if np.count_nonzero(dia == d) >= 24]   # solo dias completos
    tmin = np.array([[T[k, dia == d].min() for d in dias] for k in range(len(T))])   # miembro x dia
    tmax = np.array([[T[k, dia == d].max() for d in dias] for k in range(len(T))])
    cen = next(e['miembro'] for e in elegidas if e['nombre'] == 'central')
    Tc = series[cen].T_aire_C[dentro]
    x = np.array([datetime.combine(d, datetime.min.time()) + timedelta(hours=12) for d in dias])
    for datos, c, txt in ((tmax, '#cd4c44', 'Máxima del día'), (tmin, '#287bb5', 'Mínima del día')):
        a2.vlines(x, datos.min(0), datos.max(0), color=c, lw=6, alpha=.3)
        y = [(Tc[dia == d].max() if c == '#cd4c44' else Tc[dia == d].min()) for d in dias]
        a2.plot(x, y, 'o-', color=c, ms=4, lw=1.2, label=f'{txt} (trayectoria central)')
    a2.set_ylabel('Mínima y máxima diaria [°C]')
    a2.set_ylim(T.min() - 6, T.max() + 2)
    a2.legend(fontsize=8, loc='lower right', ncol=2, title='Barra = rango de los 51 miembros · solo días completos', title_fontsize=8)
    a2.xaxis.set_major_locator(mdates.DayLocator(interval=2)); a2.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))
    finf = t0 + timedelta(hours=fin)
    for ax in (a1, a2):
        ax.axvspan(finf, fechas[-1], color='#9aa3a7', alpha=.15, lw=0); ax.grid(alpha=.2)
    a1.text(finf + timedelta(hours=6), T.min() - 6, 'Sin pronóstico: promedio ERA5 2021-2025\n(igual para todos los miembros)',
            fontsize=8, color='#555', va='bottom')
    a2.set_xlabel(f"Fecha (hora de Chile) · desde el {t0:%d-%m-%Y %H:%M}, inicio de la vertida, durante {H / 24:.0f} días")
    bajado = datetime.strptime(snap['nombre'][:10], '%Y-%m-%d')
    fig.suptitle(f"Temperatura ambiente en Renaico · ensemble ECMWF de {len(series)} miembros descargado el {bajado:%d-%m-%Y}")
    fig.savefig(ruta, dpi=140); plt.close(fig)
    print(f'Aire: minima {T.min():.1f} C el {fechas[i_min[1]]:%d-%m %H:%M}, maxima {T.max():.1f} C el {fechas[i_max[1]]:%d-%m %H:%M}. Figura: {ruta}')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=RAIZ / 'config' / 'riesgo_e13.yaml')
    p.add_argument('--descargar', action='store_true', help='Baja un ensemble nuevo (requiere internet)')
    p.add_argument('--salida', type=Path, default=RAIZ / 'resultados' / 'e13' / 'clima')
    p.add_argument('--descargar-suelo', action='store_true',
                   help='Baja la temperatura del suelo de ERA5 (dias recientes) para el estado inicial del suelo')
    args = p.parse_args()
    cfg = bim.cargar_yaml(args.config)
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    if args.descargar_suelo:
        descargar_suelo(cfg, caso)
    if args.descargar:
        t0 = datetime.fromisoformat(cfg['t0'])
        historico = dict(desde=(t0 - timedelta(days=1)).strftime('%m-%d'),
                         dias=int(np.ceil(horizonte_maximo(cfg) / 24)) + 3, anos=list(cfg['clima']['historico_anos']))
        sitio = caso['sitio']
        carpeta = descargar_ensemble(RAIZ / cfg['clima']['carpeta'], sitio['latitud'], sitio['longitud'],
                                     zona=sitio['zona_horaria'], modelo=cfg['clima']['modelo'], dias=int(cfg['clima']['dias']),
                                     historico=historico)
        print(f'Ensemble guardado en {carpeta}')
    snap, series = series_del_ensemble(cfg)
    elegidas, extra = elegir_trayectorias(series, int(cfg['clima']['trayectorias']))
    s0 = next(iter(series.values()))
    salida = dict(snapshot=snap['nombre'], carpeta=str(Path(snap['carpeta']).relative_to(RAIZ)).replace('\\', '/'),
                  metadata=snap['metadata'], t0=cfg['t0'], horizonte_h=horizonte_maximo(cfg),
                  horas_pronostico=s0.horas_pronostico, horas_relleno=s0.horas_relleno, metodo_cola=s0.metodo_cola,
                  criterio='temperatura media del aire en las horas con pronostico; 3 grupos de igual tamano; '
                           'representante = miembro mas cercano a la mediana del grupo',
                  trayectorias=elegidas, frio_extremo=extra['frio_extremo'], calido_extremo=extra['calido_extremo'],
                  indice_por_miembro_C=extra['indice_C'])
    args.salida.mkdir(parents=True, exist_ok=True)
    (args.salida / 'trayectorias.json').write_text(json.dumps(salida, ensure_ascii=False, indent=2), encoding='utf-8')

    fig, axs = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True, layout='constrained')
    claves = [('T_aire_C', 'Aire a 2 m [°C]'), ('viento_m_s', 'Viento a 10 m [m/s]'), ('ghi_W_m2', 'Radiación [W/m²]')]
    for m, s in series.items():
        for ax, (c, _) in zip(axs, claves):
            ax.plot(s.t_s / 3600, getattr(s, c), color='#b8bec2', lw=.5, alpha=.5)
    for e in elegidas:
        s = series[e['miembro']]
        for ax, (c, et) in zip(axs, claves):
            ax.plot(s.t_s / 3600, getattr(s, c), color=COLORES.get(e['nombre'], 'k'), lw=1.5,
                    label=f"{e['nombre']}: {e['miembro']} (peso {e['peso']:.3f})"); ax.set_ylabel(et)
    fin = s0.t_s[s0.fuente == 'pronostico'].max() / 3600
    for ax in axs:
        ax.axvspan(fin, s0.t_s[-1] / 3600, color='#9aa3a7', alpha=.15, lw=0); ax.grid(alpha=.2)
    axs[0].legend(fontsize=8); axs[-1].set_xlabel('Horas desde t0')
    fig.suptitle(f"{len(series)} miembros del ensemble (gris) y las {len(elegidas)} trayectorias elegidas · "
                 f"sombreado: completado con ERA5 desde las {fin:.0f} h")
    fig.savefig(args.salida / 'trayectorias.png', dpi=140); plt.close(fig)
    figura_min_max(cfg, snap, series, elegidas, extra, args.salida / 'temperatura_min_max.png')
    for e in elegidas:
        print(f"{e['nombre']:8s} {e['miembro']}  peso {e['peso']:.3f}  T media {e['indice_C']:.2f} C  grupo {e['rango_C'][0]:.2f}-{e['rango_C'][1]:.2f} C")
    print(f'Pronostico: {s0.horas_pronostico} h; relleno: {s0.horas_relleno} h ({s0.metodo_cola}). Salida: {args.salida}')


if __name__ == '__main__':
    main()
