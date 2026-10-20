"""
Paso 2 — Leer un IFC desde Python (sin Blender) y mostrar lo que contiene.

    python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc
    python scripts/02_inspeccionar_ifc.py modelos/fundacion_G7_plan_v01.ifc --comparar modelos/editado_en_bonsai.ifc

Muestra: jerarquía espacial, vertidas con sus propiedades, volumen y cotas
calculados desde la geometría. Con --comparar, lista qué cambió entre dos
versiones (p. ej. después de editar el modelo en Bonsai).
"""

import argparse
from pathlib import Path

from mcplan import bim


def imprimir_jerarquia(elemento, nivel=0):
    """Recorre IfcProject > IfcSite > IfcBuilding > ... y sus elementos contenidos."""
    print("  " * nivel + f"- {elemento.is_a()}  '{elemento.Name}'")
    for rel in getattr(elemento, "IsDecomposedBy", []):        # agregación (partes)
        for hijo in rel.RelatedObjects:
            imprimir_jerarquia(hijo, nivel + 1)
    for rel in getattr(elemento, "ContainsElements", []):      # contención espacial
        for hijo in sorted(rel.RelatedElements, key=lambda e: e.Name or ""):
            print("  " * (nivel + 1) + f"- {hijo.is_a()}  '{hijo.Name}'  GlobalId={hijo.GlobalId}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ifc", type=Path)
    parser.add_argument("--comparar", type=Path, default=None, help="otra versión del IFC")
    args = parser.parse_args()

    modelo = bim.abrir(args.ifc)
    print(f"Archivo: {args.ifc}   esquema: {modelo.schema}\n")

    print("JERARQUÍA")
    for proyecto in modelo.by_type("IfcProject"):
        imprimir_jerarquia(proyecto)

    print("\nVERTIDAS (property set MC_Vertida)")
    geometria = bim.geometria_vertidas(modelo)
    for v in bim.leer_vertidas(modelo):
        g = geometria[v["VertidaID"]]
        print(f"\n  {v['VertidaID']}  ({v['ClaseIFC']}, material: {v['Material']})")
        for campo in bim.PROPIEDADES_VERTIDA[1:]:
            print(f"    {campo:<24} {v.get(campo, '¡FALTA!')}")
        print(f"    {'Volumen [m³]':<24} {g['volumen_m3']:.2f}")
        print(f"    {'Cota z [m]':<24} {g['min_xyz_m'][2]:.2f} -> {g['max_xyz_m'][2]:.2f}")
        faltan = bim.propiedades_faltantes(v)
        if faltan:
            print(f"    ATENCIÓN: faltan {faltan}")
    total = sum(g["volumen_m3"] for g in geometria.values())
    print(f"\n  Volumen total de hormigón: {total:.1f} m³")

    if args.comparar:
        print(f"\nDIFERENCIAS  {args.ifc.name}  ->  {args.comparar.name}")
        cambios = bim.comparar_vertidas(modelo, args.comparar)
        print("\n".join("  " + c for c in cambios) if cambios else "  (sin cambios en las vertidas)")


if __name__ == "__main__":
    main()
