"""Copia la geometria del IFC y el ultimo snapshot de clima dentro de visor_laboratorio.html.

El visor calcula en el navegador, pero no puede leer el IFC ni la carpeta data/weather:
lleva una copia de los triangulos, del programa de cada vertida y de un snapshot de
Open-Meteo. Si el IFC cambia o se descarga un snapshot nuevo, corra:

python -X utf8 scripts/actualizar_visor_laboratorio.py
python -X utf8 scripts/actualizar_visor_laboratorio.py --comprobar

--comprobar no escribe nada: termina con error si alguna copia quedo vieja.
"""
from pathlib import Path
import argparse
import json
import re
import sys
import numpy as np
from mcplan.voxelizacion import extraer_mallas
from mcplan.clima import cargar_snapshot, ultimo_snapshot

RAIZ = Path(__file__).resolve().parents[1]
VISOR = RAIZ / 'visor_laboratorio.html'
BLOQUE = re.compile(r'(<script type="application/json" id="geometria">)(.*?)(</script>)', re.S)
BLOQUE_CLIMA = re.compile(r'(<script type="application/json" id="clima">)(.*?)(</script>)', re.S)


def geometria_desde_ifc(ruta_ifc):
    """Triangulos, volumen y programa de cada vertida, tal como los usa el visor."""
    datos = []
    for malla in extraer_mallas(ruta_ifc):
        p = malla.propiedades
        datos.append(dict(id=malla.vertida_id, nombre=malla.nombre,
                          v=np.round(malla.vertices, 6).ravel().tolist(), t=malla.triangulos.ravel().tolist(),
                          vol=malla.volumen_ifc_m3, inicio=float(p['Inicio_h']), temp=float(p['TemperaturaColocacion_C']),
                          mezcla=p['TipoHormigon'], elemento=p.get('Elemento')))
    return datos


def clima_desde_snapshot(carpeta):
    """Solo lo que usa el visor: metadata, horas y las tres variables."""
    s = cargar_snapshot(carpeta)
    corto = lambda d: dict(utc_offset_seconds=d.get('utc_offset_seconds', 0), hourly=d['hourly'])
    return dict(nombre=s['nombre'], metadata=s['metadata'], pronostico=corto(s['pronostico']),
                historico={a: corto(d) for a, d in s['historico'].items()})


def _bloque(regex, texto, nombre):
    encontrado = regex.search(texto)
    if not encontrado:
        raise ValueError(f'El visor no tiene el bloque de {nombre}.')
    return json.loads(encontrado.group(2))


def geometria_del_visor(texto):
    return _bloque(BLOQUE, texto, 'geometria')


def clima_del_visor(texto):
    return _bloque(BLOQUE_CLIMA, texto, 'clima')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ifc', type=Path, default=RAIZ / 'fundaciones' / 'fundacion.ifc')
    p.add_argument('--snapshot', type=Path, help='Carpeta del snapshot; por defecto el ultimo de data/weather')
    p.add_argument('--visor', type=Path, default=VISOR)
    p.add_argument('--comprobar', action='store_true', help='Solo revisa si las copias coinciden')
    args = p.parse_args()
    texto = args.visor.read_text(encoding='utf-8')
    geometria = geometria_desde_ifc(args.ifc)
    clima = clima_desde_snapshot(args.snapshot or ultimo_snapshot(RAIZ / 'data' / 'weather'))
    if args.comprobar:
        if geometria_del_visor(texto) != json.loads(json.dumps(geometria)):
            sys.exit('El visor tiene una geometria distinta al IFC. Corra este script sin --comprobar.')
        if clima_del_visor(texto) != json.loads(json.dumps(clima)):
            sys.exit('El visor tiene otro snapshot de clima. Corra este script sin --comprobar.')
        print(f"El visor coincide con el IFC y con el snapshot {clima['nombre']}.")
        return
    compacto = lambda x: json.dumps(x, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    texto = BLOQUE.sub(lambda m: m.group(1) + compacto(geometria) + m.group(3), texto, count=1)
    texto = BLOQUE_CLIMA.sub(lambda m: m.group(1) + compacto(clima) + m.group(3), texto, count=1)
    args.visor.write_text(texto, encoding='utf-8', newline='\n')
    print(f"Visor actualizado: {len(geometria)} vertidas y snapshot {clima['nombre']}.")


if __name__ == '__main__':
    main()
