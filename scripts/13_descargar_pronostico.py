"""E10: descarga un snapshot meteorologico nuevo de Open-Meteo para el sitio del curso.

python -X utf8 scripts/13_descargar_pronostico.py
python -X utf8 scripts/13_descargar_pronostico.py --sin-historico

Guarda data/weather/<fecha UTC de descarga>/ con metadata.json, forecast.json e
historico.json (mismos dias de anos anteriores, para completar lo que el pronostico
no alcanza). Nunca sobrescribe: cada descarga queda como evidencia de ese dia.
Usa la latitud y longitud de config/caso_base.yaml (regla del curso).
"""
from datetime import datetime, timedelta
from pathlib import Path
import argparse
import numpy as np
from mcplan import bim
from mcplan.clima import descargar_snapshot, cargar_snapshot, serie_clima

RAIZ = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=RAIZ / 'config' / 'ambiente_e10.yaml')
    p.add_argument('--sin-historico', action='store_true', help='Solo el pronostico')
    args = p.parse_args()
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    amb = bim.cargar_yaml(args.config)
    sitio, clima = caso['sitio'], amb['clima']
    inicio = datetime.fromisoformat(amb['inicio_primera_vertida'])
    historico = None
    if not args.sin_historico:
        desde = inicio - timedelta(days=1)
        historico = dict(desde=desde.strftime('%m-%d'), dias=int(np.ceil(amb['horizonte_h'] / 24)) + 2,
                         anos=list(clima['historico_anos']))
    carpeta = descargar_snapshot(RAIZ / clima['carpeta'], sitio['latitud'], sitio['longitud'],
                                 zona=sitio['zona_horaria'], dias=int(clima['dias_pronostico']), historico=historico)
    snap = cargar_snapshot(carpeta)
    horas = snap['pronostico']['hourly']['time']
    print(f'Snapshot: {carpeta}')
    print(f"Descargado: {snap['metadata']['retrieved_at']} | pronostico de {horas[0]} a {horas[-1]} (hora local)")
    if historico:
        s = serie_clima(snap, amb['inicio_primera_vertida'], amb['horizonte_h'], clima['metodo_cola'])
        print(f"Horizonte desde {amb['inicio_primera_vertida']}: {s.horas_pronostico} h con pronostico, "
              f'{s.horas_relleno} h completadas con {s.metodo_cola}.')


if __name__ == '__main__':
    main()
