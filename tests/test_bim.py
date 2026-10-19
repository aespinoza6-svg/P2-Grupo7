"""
Tests del modelo BIM/IFC (mcplan.bim).

Comprueban que el IFC generado desde config/ contiene lo que exige el
enunciado (vertidas identificables con sus propiedades), que la geometría
tiene los volúmenes esperados y que un cambio en el modelo se detecta
desde Python (round-trip: escribir -> leer -> editar -> comparar).
"""

from pathlib import Path

import numpy as np
import pytest

from mcplan import bim

RAIZ = Path(__file__).resolve().parents[1]
CONFIG = RAIZ / "config"


@pytest.fixture(scope="module")
def entradas():
    return (
        bim.cargar_yaml(CONFIG / "sitio.yaml"),
        bim.cargar_yaml(CONFIG / "fundacion.yaml"),
        bim.cargar_yaml(CONFIG / "plan_v01.yaml"),
    )


@pytest.fixture()
def ruta_ifc(entradas, tmp_path):
    """Genera el IFC en una carpeta temporal y devuelve su ruta."""
    ruta = tmp_path / "fundacion.ifc"
    bim.crear_modelo_fundacion(*entradas).write(str(ruta))
    return ruta


def test_tres_vertidas_con_todas_las_propiedades(ruta_ifc):
    vertidas = bim.leer_vertidas(ruta_ifc)
    assert [v["VertidaID"] for v in vertidas] == ["V01", "V02", "V03"]
    for v in vertidas:
        assert bim.propiedades_faltantes(v) == []
        assert v["ClaseIFC"] == "IfcFooting"
        assert v["Material"] == "Hormigon G30"


def test_inicio_planificado_se_calcula_desde_el_plan(ruta_ifc):
    """V02 parte 72 h después de fecha_inicio (2026-12-01 08:00)."""
    v02 = bim.leer_vertidas(ruta_ifc)[1]
    assert v02["InicioPlanificado"] == "2026-12-04T08:00:00-03:00"


def test_volumenes_coinciden_con_la_geometria(ruta_ifc, entradas):
    """Volumen calculado por IfcOpenShell vs fórmula exacta.
    El cilindro se aproxima con triángulos: se acepta 1 % de diferencia."""
    _, fundacion, plan = entradas
    g = fundacion["geometria"]
    geo = bim.geometria_vertidas(ruta_ifc)
    for v in plan["vertidas"]:
        espesor = v["z_superior_m"] - v["z_inferior_m"]
        if v["elemento"] == "losa":
            exacto = g["losa_largo_x_m"] * g["losa_ancho_y_m"] * espesor
        else:
            exacto = np.pi * (g["pedestal_diametro_m"] / 2) ** 2 * espesor
        assert geo[v["id"]]["volumen_m3"] == pytest.approx(exacto, rel=0.01)
        assert geo[v["id"]]["min_xyz_m"][2] == pytest.approx(v["z_inferior_m"])
        assert geo[v["id"]]["max_xyz_m"][2] == pytest.approx(v["z_superior_m"])


def test_altura_total_de_la_fundacion(ruta_ifc, entradas):
    geo = bim.geometria_vertidas(ruta_ifc)
    z_max = max(d["max_xyz_m"][2] for d in geo.values())
    assert z_max == pytest.approx(entradas[1]["geometria"]["altura_total_m"])


def test_coordenadas_del_sitio(ruta_ifc):
    modelo = bim.abrir(ruta_ifc)
    sitio = modelo.by_type("IfcSite")[0]
    assert tuple(sitio.RefLatitude)[:2] == (-37, -41)
    assert tuple(sitio.RefLongitude)[:2] == (-72, -37)


def test_cambio_de_propiedad_se_detecta(ruta_ifc, tmp_path):
    """Simula lo que haríamos en Bonsai: abrir, cambiar la temperatura de
    colocación de V02, guardar como otra versión y comparar desde Python."""
    modelo = bim.abrir(ruta_ifc)
    bim.editar_propiedad(modelo, "V02", "TemperaturaColocacion_C", 15.0)
    ruta_v2 = tmp_path / "fundacion_v2.ifc"
    modelo.write(str(ruta_v2))

    cambios = bim.comparar_vertidas(ruta_ifc, ruta_v2)
    assert cambios == ["MODIFICADA V02.TemperaturaColocacion_C: 20.0 -> 15.0"]


def test_globalid_estable_al_regenerar(entradas):
    """Regenerar el IFC desde config/ entrega las mismas identidades
    (requisito para poder comparar versiones y para la trazabilidad)."""
    ids_1 = [v["GlobalId"] for v in bim.leer_vertidas(bim.crear_modelo_fundacion(*entradas))]
    ids_2 = [v["GlobalId"] for v in bim.leer_vertidas(bim.crear_modelo_fundacion(*entradas))]
    assert ids_1 == ids_2
    assert len(set(ids_1)) == 3  # y son distintas entre sí


def test_plan_con_hueco_es_rechazado(entradas):
    """Si las vertidas no llenan la fundación, se avisa (no se inventa nada)."""
    sitio, fundacion, plan = entradas
    plan_malo = {**plan, "vertidas": [dict(v) for v in plan["vertidas"]]}
    plan_malo["vertidas"][1]["z_inferior_m"] = 1.20  # deja un hueco 1.15-1.20
    with pytest.raises(ValueError, match="hueco o traslape"):
        bim.crear_modelo_fundacion(sitio, fundacion, plan_malo)
