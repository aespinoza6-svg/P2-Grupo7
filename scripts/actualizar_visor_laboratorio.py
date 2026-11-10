"""Copia la geometria del IFC dentro de visor_laboratorio.html.

El visor calcula en el navegador, pero no puede leer el IFC: lleva una copia
de los triangulos y del programa de cada vertida. Si el IFC cambia, corra:

python -X utf8 scripts/12_actualizar_visor_laboratorio.py
python -X utf8 scripts/12_actualizar_visor_laboratorio.py --comprobar

--comprobar no escribe nada: termina con error si la copia quedo vieja.
"""
from pathlib import Path
import argparse
import json
import re
import sys
import numpy as np
from mcplan.voxelizacion import extraer_mallas

RAIZ = Path(__file__).resolve().parents[1]
VISOR = RAIZ / 'visor_laboratorio.html'
BLOQUE = re.compile(r'(<script type="application/json" id="geometria">)(.*?)(</script>)', re.S)


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


def geometria_del_visor(texto):
    encontrado = BLOQUE.search(texto)
    if not encontrado:
        raise ValueError('El visor no tiene el bloque de geometria.')
    return json.loads(encontrado.group(2))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ifc', type=Path, default=RAIZ / 'fundaciones' / 'fundacion.ifc')
    p.add_argument('--visor', type=Path, default=VISOR)
    p.add_argument('--comprobar', action='store_true', help='Solo revisa si la copia coincide con el IFC')
    args = p.parse_args()
    texto = args.visor.read_text(encoding='utf-8')
    nueva = geometria_desde_ifc(args.ifc)
    if args.comprobar:
        if geometria_del_visor(texto) != json.loads(json.dumps(nueva)):
            sys.exit('El visor tiene una geometria distinta al IFC. Corra este script sin --comprobar.')
        print('El visor coincide con el IFC.')
        return
    contenido = json.dumps(nueva, ensure_ascii=False, separators=(',', ':'))
    args.visor.write_text(BLOQUE.sub(lambda m: m.group(1) + contenido + m.group(3), texto, count=1), encoding='utf-8')
    print(f'Geometria actualizada en {args.visor} ({len(nueva)} vertidas).')


if __name__ == '__main__':
    main()
