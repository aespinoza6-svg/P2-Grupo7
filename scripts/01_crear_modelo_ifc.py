"""
Paso 1 — Crear el modelo IFC de la fundación a partir de config/.

    python scripts/01_crear_modelo_ifc.py
    python scripts/01_crear_modelo_ifc.py --plan config/plan_v01.yaml --salida modelos/otro.ifc

Lee:   config/caso_base.yaml (sitio, geometría, hormigón) y el plan de vertidas.
Crea:  modelos/fundacion_G7_<nombre del plan>.ifc  (se abre en Bonsai).
"""

import argparse
from pathlib import Path

from mcplan import bim

RAIZ = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan", default=RAIZ / "config" / "plan_v01.yaml", type=Path)
    parser.add_argument("--salida", default=None, type=Path)
    args = parser.parse_args()

    caso = bim.cargar_yaml(RAIZ / "config" / "caso_base.yaml")
    plan = bim.cargar_yaml(args.plan)
    salida = args.salida or RAIZ / "modelos" / f"fundacion_G7_{plan['nombre']}.ifc"

    modelo = bim.crear_modelo_fundacion(caso, plan)
    salida.parent.mkdir(parents=True, exist_ok=True)
    modelo.write(str(salida))

    print(f"Modelo IFC creado: {salida}")
    sitio = caso["sitio"]
    print(f"  esquema: {modelo.schema}   sitio: {sitio['localidad']} ({sitio['latitud']}, {sitio['longitud']})")
    for v in bim.leer_vertidas(modelo):
        print(f"  {v['VertidaID']}  {v['Nombre']:<16} inicio {v['InicioPlanificado']}  "
              f"T_col = {v['TemperaturaColocacion_C']} °C")


if __name__ == "__main__":
    main()
