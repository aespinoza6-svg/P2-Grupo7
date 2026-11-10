"""Clima del sitio: snapshots de pronostico, historico y serie para el solver.

Fuente: Open-Meteo (guia G07). Tres tipos de dato que NO son lo mismo:
- pronostico (forecast): lo que un modelo predice, depende de CUANDO se descargo;
- historico (reanalysis ERA5, archive-api): reconstruccion del pasado;
- observacion: medicion de una estacion (no se usa aun).

Cada descarga se guarda como snapshot inmutable en data/weather/<fecha de descarga>/
con metadata.json (fuente, modelo, fecha, coordenadas, unidades) y los JSON crudos.
Nunca se sobrescribe: un pronostico viejo es evidencia de lo que se sabia ese dia.

Variables (solo las que entran a una ecuacion):
- temperature_2m [degC]: temperatura del aire en Robin;
- wind_speed_10m [m/s]: coeficiente de conveccion h(V);
- shortwave_radiation [W/m2]: flujo solar absorbido en caras horizontales.
"""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import json
import numpy as np

VARIABLES = ('temperature_2m', 'wind_speed_10m', 'shortwave_radiation')
UNIDADES = {'temperature_2m': 'degC', 'wind_speed_10m': 'm/s', 'shortwave_radiation': 'W/m2'}
URL_PRONOSTICO = 'https://api.open-meteo.com/v1/forecast'
URL_HISTORICO = 'https://archive-api.open-meteo.com/v1/archive'


# ---------------------------------------------------------------------------
# Descarga (red) -> snapshot (disco). Separado de la fisica (guia G08, principio 1).
# ---------------------------------------------------------------------------

def _pedir(url, parametros, espera_s=60):
    direccion = url + '?' + urlencode(parametros)
    with urlopen(direccion, timeout=espera_s) as r:
        return direccion, json.loads(r.read().decode('utf-8'))


def descargar_snapshot(carpeta_base, latitud, longitud, *, zona='America/Santiago', dias=16,
                       historico=None, ahora=None, pedir=_pedir):
    """Descarga el pronostico (y opcionalmente historico) y lo guarda como snapshot nuevo.

    historico: dict(desde='MM-DD', dias=int, anos=[2021, ...]) para completar la cola.
    Devuelve la carpeta creada. Si ya existe una carpeta con la misma hora, falla:
    los snapshots no se sobrescriben.
    """
    ahora = ahora or datetime.now(timezone.utc)
    nombre = ahora.strftime('%Y-%m-%dT%H%M')
    carpeta = Path(carpeta_base) / nombre
    if carpeta.exists():
        raise FileExistsError(f'Ya existe el snapshot {carpeta}; no se sobrescribe.')
    comun = dict(latitude=latitud, longitude=longitud, hourly=','.join(VARIABLES),
                 wind_speed_unit='ms', timezone=zona)
    url, pronostico = pedir(URL_PRONOSTICO, {**comun, 'forecast_days': dias, 'models': 'best_match'})
    pedidos = [dict(tipo='pronostico', url=url)]
    historicos = {}
    if historico:
        for ano in historico['anos']:
            inicio = datetime.strptime(f"{ano}-{historico['desde']}", '%Y-%m-%d')
            fin = inicio + timedelta(days=int(historico['dias']))
            u, datos = pedir(URL_HISTORICO, {**comun, 'start_date': inicio.strftime('%Y-%m-%d'),
                                                         'end_date': fin.strftime('%Y-%m-%d')})
            historicos[str(ano)] = datos
            pedidos.append(dict(tipo='historico', ano=ano, url=u))
    metadata = dict(source='Open-Meteo', model='best_match (pronostico); ERA5 reanalysis via archive-api (historico)',
                    retrieved_at=ahora.isoformat(timespec='seconds'), latitude=latitud, longitude=longitud,
                    timezone=zona, variables=list(VARIABLES), units=dict(UNIDADES),
                    forecast_days=dias, requests=pedidos,
                    nota='Snapshot inmutable. Un pronostico descargado otro dia es otro dato.')
    carpeta.mkdir(parents=True)
    (carpeta / 'metadata.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
    (carpeta / 'forecast.json').write_text(json.dumps(pronostico, ensure_ascii=False), encoding='utf-8')
    if historicos:
        (carpeta / 'historico.json').write_text(json.dumps(historicos, ensure_ascii=False), encoding='utf-8')
    return carpeta


# ---------------------------------------------------------------------------
# Snapshot (disco) -> serie horaria para el solver
# ---------------------------------------------------------------------------

def cargar_snapshot(carpeta):
    carpeta = Path(carpeta)
    metadata = json.loads((carpeta / 'metadata.json').read_text(encoding='utf-8'))
    faltan = [c for c in ('source', 'model', 'retrieved_at', 'latitude', 'longitude', 'timezone', 'variables', 'units')
              if c not in metadata]
    if faltan:
        raise ValueError(f'metadata.json incompleta: faltan {faltan}.')
    pronostico = json.loads((carpeta / 'forecast.json').read_text(encoding='utf-8'))
    ruta = carpeta / 'historico.json'
    historico = json.loads(ruta.read_text(encoding='utf-8')) if ruta.is_file() else {}
    return dict(carpeta=str(carpeta), nombre=carpeta.name, metadata=metadata, pronostico=pronostico, historico=historico)


def ultimo_snapshot(carpeta_base):
    carpetas = sorted(p for p in Path(carpeta_base).iterdir() if (p / 'metadata.json').is_file())
    if not carpetas:
        raise FileNotFoundError(f'No hay snapshots en {carpeta_base}. Corra scripts/13_descargar_pronostico.py.')
    return carpetas[-1]


def _horas(datos):
    """Instantes UTC (epoch s) de una respuesta Open-Meteo con hora local + desfase."""
    desfase = int(datos.get('utc_offset_seconds', 0))
    t = [datetime.strptime(s, '%Y-%m-%dT%H:%M').replace(tzinfo=timezone.utc).timestamp() - desfase
         for s in datos['hourly']['time']]
    return np.array(t)


@dataclass
class SerieClima:
    """Serie horaria desde el inicio de la primera vertida (t = 0)."""
    t_s: np.ndarray
    T_aire_C: np.ndarray
    viento_m_s: np.ndarray
    ghi_W_m2: np.ndarray
    fuente: np.ndarray          # 'pronostico' o 'historico' por hora
    inicio_utc: float
    snapshot: str
    metodo_cola: str

    def en(self, t_s):
        """Valores interpolados (lineal) en el tiempo t del solver."""
        t = np.clip(t_s, self.t_s[0], self.t_s[-1])
        return (float(np.interp(t, self.t_s, self.T_aire_C)), float(np.interp(t, self.t_s, self.viento_m_s)),
                float(np.interp(t, self.t_s, self.ghi_W_m2)))

    @property
    def horas_pronostico(self):
        return int(np.count_nonzero(self.fuente == 'pronostico'))

    @property
    def horas_relleno(self):
        return int(np.count_nonzero(self.fuente != 'pronostico'))


def _sin_huecos(t, valores):
    v = np.array([np.nan if x is None else x for x in valores], dtype=float)
    ok = np.isfinite(v)
    if not ok.any():
        return v
    return np.interp(t, t[ok], v[ok])


def serie_clima(snapshot, inicio_iso, horizonte_h, metodo_cola='historico_promedio'):
    """Arma la serie horaria [inicio, inicio + horizonte] para el solver.

    Usa el pronostico donde alcanza. Donde no, completa la cola con el historico de los
    mismos dias y horas de anos anteriores:
    - 'historico_promedio': promedio de los anos guardados;
    - 'historico_frio': el ano con menor temperatura media en la ventana (escenario conservador).
    """
    if metodo_cola not in ('historico_promedio', 'historico_frio'):
        raise ValueError(f'Metodo de cola desconocido: {metodo_cola}.')
    inicio = datetime.fromisoformat(inicio_iso).astimezone(timezone.utc).timestamp()
    t_abs = inicio + 3600.0 * np.arange(0, int(np.ceil(horizonte_h)) + 1)
    P = snapshot['pronostico']
    tp = _horas(P)
    fuente = np.full(len(t_abs), 'historico', dtype=object)
    series = {}
    for var in VARIABLES:
        valores = _sin_huecos(tp, P['hourly'][var])
        dentro = (t_abs >= tp[0]) & (t_abs <= tp[-1])
        s = np.full(len(t_abs), np.nan)
        s[dentro] = np.interp(t_abs[dentro], tp, valores)
        series[var] = s
    fuente[(t_abs >= tp[0]) & (t_abs <= tp[-1])] = 'pronostico'
    falta = fuente != 'pronostico'
    if falta.any():
        if not snapshot['historico']:
            raise ValueError('El pronostico no cubre el horizonte y el snapshot no trae historico para completarlo.')
        anos = sorted(snapshot['historico'])
        cubos = {}
        for ano in anos:
            H = snapshot['historico'][ano]
            th = _horas(H)
            # Se lleva cada ano al ano del inicio para comparar el mismo dia y la misma hora.
            corrimiento = datetime.fromtimestamp(inicio, timezone.utc).year - int(ano)
            th = np.array([_mover_ano(x, corrimiento) for x in th])
            cubos[ano] = {var: np.interp(t_abs, th, _sin_huecos(th, H['hourly'][var]), left=np.nan, right=np.nan)
                          for var in VARIABLES}
        if metodo_cola == 'historico_frio':
            elegido = min(anos, key=lambda a: np.nanmean(cubos[a]['temperature_2m'][falta]))
            relleno = cubos[elegido]
            etiqueta = f'historico_frio:{elegido}'
        else:
            relleno = {var: np.nanmean([cubos[a][var] for a in anos], axis=0) for var in VARIABLES}
            etiqueta = 'historico_promedio:' + ','.join(anos)
        for var in VARIABLES:
            series[var][falta] = relleno[var][falta]
        fuente[falta] = 'historico'
    else:
        etiqueta = 'sin_relleno'
    for var in VARIABLES:
        if not np.isfinite(series[var]).all():
            raise ValueError(f'Quedan horas sin dato de {var}: amplie el historico descargado.')
    return SerieClima(t_s=t_abs - inicio, T_aire_C=series['temperature_2m'], viento_m_s=series['wind_speed_10m'],
                      ghi_W_m2=np.maximum(series['shortwave_radiation'], 0.0), fuente=fuente.astype(str),
                      inicio_utc=inicio, snapshot=snapshot['nombre'], metodo_cola=etiqueta)


def _mover_ano(epoch, anos):
    d = datetime.fromtimestamp(epoch, timezone.utc)
    try:
        return d.replace(year=d.year + anos).timestamp()
    except ValueError:      # 29 de febrero
        return (d + timedelta(days=365 * anos)).timestamp()


def serie_constante(T_aire_C, horizonte_h, viento_m_s=0.0, ghi_W_m2=0.0):
    """Clima constante (el caso de E07) con la misma interfaz que un snapshot."""
    t = 3600.0 * np.arange(0, int(np.ceil(horizonte_h)) + 1)
    uno = np.ones(len(t))
    return SerieClima(t, T_aire_C * uno, viento_m_s * uno, ghi_W_m2 * uno, np.full(len(t), 'constante'),
                      0.0, 'constante', 'sin_relleno')
