"""Copia dentro de visor_laboratorio.html lo que el navegador no puede leer del repo.

El visor calcula en el navegador, pero no puede abrir archivos. Lleva una copia de:
- la geometria y el programa de cada vertida (IFC),
- el ultimo snapshot de clima (data/weather),
- las 8 curvas de hidratacion del ensayo (data/hidratacion/source_terms.txt, E13),
- los resultados de riesgo de E13 (resultados/e13): se muestran, no se recalculan.
Si cambia alguno de ellos, corra:

python -X utf8 scripts/actualizar_visor_laboratorio.py
python -X utf8 scripts/actualizar_visor_laboratorio.py --comprobar

--comprobar no escribe nada: termina con error si alguna copia quedo vieja.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import sys
import numpy as np
from mcplan import bim
from mcplan.voxelizacion import extraer_mallas
from mcplan.clima import cargar_snapshot, ultimo_snapshot
from mcplan.hidratacion import leer_curvas_tabuladas

RAIZ = Path(__file__).resolve().parents[1]
VISOR = RAIZ / 'visor_laboratorio.html'
CURVAS = RAIZ / 'data' / 'hidratacion' / 'source_terms.txt'
E13 = RAIZ / 'resultados' / 'e13'
BLOQUE = re.compile(r'(<script type="application/json" id="geometria">)(.*?)(</script>)', re.S)
BLOQUE_CLIMA = re.compile(r'(<script type="application/json" id="clima">)(.*?)(</script>)', re.S)
BLOQUE_CURVAS = re.compile(r'(<script type="application/json" id="curvas">)(.*?)(</script>)', re.S)
BLOQUE_RIESGO = re.compile(r'(<script type="application/json" id="riesgo">)(.*?)(</script>)', re.S)


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
    corto = lambda d: dict(utc_offset_seconds=d.get('utc_offset_seconds', 0), hourly=d['hourly'],
                           **{k: d[k] for k in ('latitude', 'longitude', 'elevation', 'timezone', 'hourly_units') if k in d})
    return dict(nombre=s['nombre'], metadata=s['metadata'], pronostico=corto(s['pronostico']),
                historico={a: corto(d) for a, d in s['historico'].items()})


def curvas_desde_archivo(ruta=CURVAS):
    """Edad [h] y Q [W/m3] de cada curva, sin los NaN del final; rho y cp del encabezado."""
    curvas = leer_curvas_tabuladas(ruta)
    c0 = next(iter(curvas.values()))
    return dict(archivo=Path(ruta).relative_to(RAIZ).as_posix(), sha256=hashlib.sha256(Path(ruta).read_bytes()).hexdigest(),
                rho=c0.rho_kg_m3, cp=c0.cp_J_kgK,
                curvas=[dict(nombre=n, a=np.round(c.a_s / 3600, 6).tolist(), q=c.q.tolist()) for n, c in curvas.items()])


def _csv(ruta):
    with open(ruta, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def _num(x):
    return None if x in ('', None) else round(float(x), 4)


def riesgo_desde_resultados(carpeta=E13):
    """Lo justo para la pestana Riesgo: planes, escenarios (MT, MD), prueba de malla y suelo a 12 C."""
    cfg = bim.cargar_yaml(RAIZ / 'config' / 'riesgo_e13.yaml')
    mallas = {}
    for nombre in ('exploracion', 'final'):
        res = json.loads((carpeta / nombre / 'resumen.json').read_text(encoding='utf-8'))
        esc = {}
        for f in _csv(carpeta / nombre / 'escenarios.csv'):
            esc.setdefault(f['plan'], []).append([f['curva'], f['clima'], _num(f['MT']), _num(f['MD'])])
        planes = [dict(plan=f['plan'], D_h=_num(f['D_h']), escenarios=int(f['escenarios']), MT_max=_num(f['MT_max']),
                       MD_max=_num(f['MD_max']), Pf=_num(f.get('Pf_base')), Pf_solo_SS=_num(f.get('Pf_solo_SS')),
                       Pf_22A=_num(f.get('Pf_22A_a_la_mitad')), cumple=f['cumple'] == 'True')
                  for f in _csv(carpeta / nombre / 'planes.csv')]
        mallas[nombre] = dict(d=cfg['mallas'][nombre], planes=planes, escenarios=esc, recomendado=res.get('recomendado'),
                              commit=res.get('commit'), clima=res.get('clima'))
    # El mismo escenario en todas las mallas: asi se ve si MD converge al refinar.
    # Dos planes: el de la prueba original (10 C) y el propuesto (5 C).
    planes_conv, curva, clima = ['unica_manta15_T10', 'unica_manta15_T5'], 'M4_SS', 'frio'
    conv = []
    for nombre, d in cfg['mallas'].items():
        ruta = carpeta / nombre / 'escenarios.csv'
        md, mt = {}, {}
        for f in (_csv(ruta) if ruta.exists() else []):
            if f['plan'] in planes_conv and (f['curva'], f['clima']) == (curva, clima):
                md[f['plan']], mt[f['plan']] = _num(f['MD']), _num(f['MT'])
        if md:
            conv.append(dict(malla=nombre, dx=d['dx'], dz=d['dz'], MD=md, MT=mt))
    conv.sort(key=lambda f: (-f['dx'], -f['dz']))
    suelo = []
    for nombre in ('exploracion', 'final'):
        ruta = carpeta / f'{nombre}_suelo12_12' / 'escenarios.csv'
        for f in (_csv(ruta) if ruta.exists() else []):
            suelo.append(dict(malla=nombre, plan=f['plan'], curva=f['curva'], clima=f['clima'], MT=_num(f['MT']), MD=_num(f['MD'])))
    return dict(criterio='P_f = P(MT > 70 °C o MΔ > 20 °C) ≤ 5 %', curvas=list(cfg['curvas']['pesos']['base']),
                climas=['frio', 'central', 'calido'], mallas=mallas,
                convergencia=dict(planes=planes_conv, curva=curva, clima=clima, filas=conv), suelo12=suelo)


def _bloque(regex, texto, nombre):
    encontrado = regex.search(texto)
    if not encontrado:
        raise ValueError(f'El visor no tiene el bloque de {nombre}.')
    return json.loads(encontrado.group(2))


def geometria_del_visor(texto):
    return _bloque(BLOQUE, texto, 'geometria')


def clima_del_visor(texto):
    return _bloque(BLOQUE_CLIMA, texto, 'clima')


def curvas_del_visor(texto):
    return _bloque(BLOQUE_CURVAS, texto, 'curvas')


def riesgo_del_visor(texto):
    return _bloque(BLOQUE_RIESGO, texto, 'riesgo')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ifc', type=Path, default=RAIZ / 'fundaciones' / 'fundacion.ifc')
    p.add_argument('--snapshot', type=Path, help='Carpeta del snapshot; por defecto el ultimo de data/weather')
    p.add_argument('--visor', type=Path, default=VISOR)
    p.add_argument('--comprobar', action='store_true', help='Solo revisa si las copias coinciden')
    args = p.parse_args()
    texto = args.visor.read_text(encoding='utf-8')
    copias = [(BLOQUE, 'geometria', geometria_desde_ifc(args.ifc), 'una geometria distinta al IFC'),
              (BLOQUE_CLIMA, 'clima', clima_desde_snapshot(args.snapshot or ultimo_snapshot(RAIZ / 'data' / 'weather')),
               'otro snapshot de clima'),
              (BLOQUE_CURVAS, 'curvas', curvas_desde_archivo(), 'otras curvas de hidratacion'),
              (BLOQUE_RIESGO, 'riesgo', riesgo_desde_resultados(), 'otros resultados de E13')]
    if args.comprobar:
        for regex, nombre, datos, problema in copias:
            if _bloque(regex, texto, nombre) != json.loads(json.dumps(datos)):
                sys.exit(f'El visor tiene {problema}. Corra este script sin --comprobar.')
        print('El visor coincide con el IFC, el snapshot de clima, las curvas y los resultados de E13.')
        return
    compacto = lambda x: json.dumps(x, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    for regex, nombre, datos, _ in copias:
        if not regex.search(texto):
            raise ValueError(f'El visor no tiene el bloque de {nombre}.')
        texto = regex.sub(lambda m: m.group(1) + compacto(datos) + m.group(3), texto, count=1)
    args.visor.write_text(texto, encoding='utf-8', newline='\n')
    print(f"Visor actualizado: {len(copias[0][2])} vertidas, snapshot {copias[1][2]['nombre']}, "
          f"{len(copias[2][2]['curvas'])} curvas y los planes de E13.")


if __name__ == '__main__':
    main()
