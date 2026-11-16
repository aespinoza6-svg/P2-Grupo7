"""E13: descarga el ensemble de Open-Meteo y lo reduce a pocas historias climaticas con peso.

python -X utf8 scripts/16_clima_ensemble.py --descargar     (una vez: guarda data/weather_ensemble/<fecha>/)
python -X utf8 scripts/16_clima_ensemble.py                 (solo vuelve a elegir y graficar)

Cada miembro del ensemble es un futuro posible completo: aire, viento y sol coherentes
entre si y con su ciclo diario. Se ordenan por la temperatura media del aire en las horas
con pronostico, se parten en tres grupos iguales (frio, central, calido) y de cada grupo se
toma el miembro mas cercano a su mediana. Peso = fraccion de miembros que representa.
La cola despues del ultimo dia del ensemble es la misma para todos: promedio ERA5.
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
    for e in elegidas:
        print(f"{e['nombre']:8s} {e['miembro']}  peso {e['peso']:.3f}  T media {e['indice_C']:.2f} C  grupo {e['rango_C'][0]:.2f}-{e['rango_C'][1]:.2f} C")
    print(f'Pronostico: {s0.horas_pronostico} h; relleno: {s0.horas_relleno} h ({s0.metodo_cola}). Salida: {args.salida}')


if __name__ == '__main__':
    main()
