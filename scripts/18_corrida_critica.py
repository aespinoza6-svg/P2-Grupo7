"""F01: repite un escenario de E13 y guarda su historia hora a hora y el campo 3D.

python -X utf8 scripts/18_corrida_critica.py                              (propuesta, M4_SS frio, malla final)
python -X utf8 scripts/18_corrida_critica.py --plan unica_moldaje         (otro plan, mismo escenario)
python -X utf8 scripts/18_corrida_critica.py --cambiar max_dt_s=300       (sensibilidad al paso temporal)
python -X utf8 scripts/18_corrida_critica.py --cambiar intervalo_salida_h=0.25

No cambia el calculo: usa la misma funcion correr() de scripts/17_riesgo_planes.py y solo
envuelve el solver para registrar, en cada salida, la T maxima del hormigon, la T minima en
S25 y MΔ, y para guardar el campo 3D en el instante del MΔ maximo. --cambiar modifica un
valor numerico de config/riesgo_e13.yaml solo en esta corrida.

Salidas en resultados/e13/criticos/:
  <nombre>.json            la misma fila que escenarios.csv (MT, MD, instantes, ubicaciones, residuo)
  <nombre>_historia.csv    t_h, T_max_C, T25_min_C, MD_C en cada salida
  <nombre>.npz             campo 3D en el instante del MΔ maximo (no se versiona, unos 2 MB)
  resumen.csv              una fila por corrida guardada en la carpeta
"""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import sys
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))
SALIDA = RAIZ / 'resultados' / 'e13' / 'criticos'


def cargar_script_17():
    spec = importlib.util.spec_from_file_location('riesgo17', RAIZ / 'scripts' / '17_riesgo_planes.py')
    modulo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(modulo)
    return modulo


def correr_con_historia(r17, cfg, plan, curva, clima, malla):
    """Llama a r17.correr() y devuelve (fila, historia, campo en el instante del MΔ maximo)."""
    import mcplan.solver_vertidas as sv
    from mcplan.riesgo import medir
    historia, campo = [], {}
    original = sv.resolver_vertidas

    def envoltura(programa, materiales, bordes, t_final, **kw):
        previo = kw.get('al_guardar')

        def registrar(t, T, activa, sistema):
            if previo:
                previo(t, T, activa, sistema)
            m = medir(T, activa, sistema, t)
            if not m:
                return
            md = m['T_max'] - m['T25_min']
            historia.append((t / 3600, m['T_max'], m['T25_min'], md))
            if md > campo.get('md', -np.inf):
                g = sistema.programa.grilla
                campo.update(md=md, t_h=t / 3600, T=np.array(T, dtype=np.float32), activa=np.array(activa),
                             vertida_id=np.array(g.vertida_id), origen=np.asarray(g.origen),
                             paso=np.asarray(g.paso), i_max=m['i_max'], i25=m['i25'])

        kw['al_guardar'] = registrar
        return original(programa, materiales, bordes, t_final, **kw)

    sv.resolver_vertidas = envoltura
    try:
        fila = r17.correr(cfg, plan, curva, clima, malla)
    finally:
        sv.resolver_vertidas = original
    return fila, np.array(historia), campo


def guardar(nombre, fila, historia, campo):
    SALIDA.mkdir(parents=True, exist_ok=True)
    (SALIDA / f'{nombre}.json').write_text(json.dumps(fila, ensure_ascii=False, indent=2), encoding='utf-8')
    with (SALIDA / f'{nombre}_historia.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['t_h', 'T_max_C', 'T25_min_C', 'MD_C'])
        for t, tm, t25, md in historia:
            w.writerow([f'{t:.4f}', f'{tm:.4f}', f'{t25:.4f}', f'{md:.4f}'])
    np.savez_compressed(SALIDA / f'{nombre}.npz', campo_T=campo['T'], campo_activa=campo['activa'],
                        campo_vertida=campo['vertida_id'], origen=campo['origen'], paso=campo['paso'],
                        t_campo_h=campo['t_h'], i_max_campo=campo['i_max'], i25_campo=campo['i25'])


def resumir():
    filas = []
    for ruta in sorted(SALIDA.glob('critico_*.json')):
        j = json.loads(ruta.read_text(encoding='utf-8'))
        filas.append(dict(corrida=ruta.stem, plan=j['plan'], curva=j['curva'], clima=j['clima'], malla=j['malla'],
                          MT=j['MT'], t_MT_h=j['t_MT_h'], MD=j['MD'], t_MD_h=j['t_MD_h'], MD_336=j['MD_336'],
                          xyz_S25_MD=j['xyz_S25_MD'], pasos=j['pasos'], residuo=j['residuo']))
    if filas:
        with (SALIDA / 'resumen.csv').open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0]))
            w.writeheader()
            w.writerows(filas)
    return filas


def main():
    from mcplan import bim
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', type=Path, default=RAIZ / 'config' / 'riesgo_e13.yaml')
    p.add_argument('--malla', default='final')
    p.add_argument('--plan', default='unica_manta15_T5')
    p.add_argument('--curva', default='M4_SS')
    p.add_argument('--clima', default='frio')
    p.add_argument('--cambiar', nargs='*', default=[], metavar='CLAVE=VALOR',
                   help='Cambia un valor numerico de la configuracion solo en esta corrida')
    args = p.parse_args()
    cfg = bim.cargar_yaml(args.config)
    cambios = dict(c.split('=') for c in args.cambiar)
    for clave, valor in cambios.items():
        if clave not in cfg:
            raise SystemExit(f'La configuracion no tiene la clave {clave!r}.')
        cfg[clave] = float(valor)
    plan = next((x for x in cfg['planes'] if x['nombre'] == args.plan), None)
    if plan is None:
        raise SystemExit(f'No existe el plan {args.plan!r} en {args.config.name}.')
    r17 = cargar_script_17()
    r17.preparar_grilla(cfg, args.malla)
    fila, historia, campo = correr_con_historia(r17, cfg, plan, args.curva, args.clima, args.malla)
    nombre = f'critico_{args.plan}_{args.curva}_{args.clima}_{args.malla}' + ''.join(f'_{k}{v}' for k, v in cambios.items())
    guardar(nombre, fila, historia, campo)
    resumir()
    print(f"{nombre}: MT {fila['MT']:.3f} °C a las {fila['t_MT_h']:.0f} h; MΔ {fila['MD']:.3f} °C a las "
          f"{fila['t_MD_h']:.1f} h; {fila['pasos']} pasos en {fila['segundos']:.0f} s. Resultados: {SALIDA}")


if __name__ == '__main__':
    main()
