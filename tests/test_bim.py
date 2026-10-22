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
    """(caso, plan): los mismos archivos que usa scripts/01_crear_modelo_ifc.py."""
    return (
        bim.cargar_yaml(CONFIG / "caso_base.yaml"),
        bim.cargar_yaml(CONFIG / "plan_v02.yaml"),
    )


@pytest.fixture()
def ruta_ifc(entradas, tmp_path):
    """Genera el IFC en una carpeta temporal y devuelve su ruta."""
    ruta = tmp_path / "fundacion.ifc"
    bim.crear_modelo_fundacion(*entradas).write(str(ruta))
    return ruta


def test_cinco_vertidas_con_todas_las_propiedades(ruta_ifc):
    vertidas = bim.leer_vertidas(ruta_ifc)
    assert [v["VertidaID"] for v in vertidas] == [1, 2, 3, 4, 5]
    for v in vertidas:
        assert bim.propiedades_faltantes(v) == []
        assert v["ClaseIFC"] == "IfcFooting"
        assert v["Material"] == "Hormigon M1"
        assert v["Inicio_h"] >= 0


def test_inicio_planificado_se_calcula_desde_el_plan(ruta_ifc):
    """V02 parte 6 h después de fecha_inicio (2026-12-01 08:00)."""
    v02 = bim.leer_vertidas(ruta_ifc)[1]
    assert v02["Inicio_h"] == 6.0
    assert v02["InicioPlanificado"] == "2026-12-01T14:00:00-03:00"


def test_volumenes_coinciden_con_la_geometria(ruta_ifc, entradas):
    """Volumen IFC frente al de prismas/frustos octagonales y cilindro faceteado."""
    caso, plan = entradas
    g = caso["geometria"]
    geo = bim.geometria_vertidas(ruta_ifc)
    n_oct = g["lados_planta"]
    n_cil = g["lados_cilindro"]
    r_ext = (g["distancia_caras_exteriores_m"] / 2) / np.cos(np.pi / n_oct)
    area_oct = lambda r: n_oct / 2 * r**2 * np.sin(2 * np.pi / n_oct)
    area_cil = lambda r: n_cil / 2 * r**2 * np.sin(2 * np.pi / n_cil)
    for v in plan["vertidas"]:
        espesor = v["z_superior_m"] - v["z_inferior_m"]
        if v["forma"] == "prisma_octagonal":
            exacto = area_oct(r_ext) * espesor
        elif v["forma"] == "tronco_octagonal":
            r0 = bim._radio_octagono_en_z(g, v["z_inferior_m"])
            r1 = bim._radio_octagono_en_z(g, v["z_superior_m"])
            a0, a1 = area_oct(r0), area_oct(r1)
            exacto = espesor / 3 * (a0 + a1 + np.sqrt(a0 * a1))
        else:
            exacto = area_cil(g["pedestal_diametro_m"] / 2) * espesor
        vid = int(v["id"][1:])
        assert geo[vid]["volumen_m3"] == pytest.approx(exacto, rel=0.001)
        assert geo[vid]["min_xyz_m"][2] == pytest.approx(v["z_inferior_m"])
        assert geo[vid]["max_xyz_m"][2] == pytest.approx(v["z_superior_m"])
    assert sum(d["volumen_m3"] for d in geo.values()) == pytest.approx(636.78, abs=0.15)


def test_altura_total_de_la_fundacion(ruta_ifc, entradas):
    geo = bim.geometria_vertidas(ruta_ifc)
    z_max = max(d["max_xyz_m"][2] for d in geo.values())
    assert z_max == pytest.approx(entradas[0]["geometria"]["altura_total_m"])


def test_planta_octagonal_mide_18_m_entre_caras(ruta_ifc):
    geo = bim.geometria_vertidas(ruta_ifc)
    lo = np.min([d["min_xyz_m"] for d in geo.values()], axis=0)
    hi = np.max([d["max_xyz_m"] for d in geo.values()], axis=0)
    assert hi - lo == pytest.approx([18.0, 18.0, 4.30], abs=0.01)


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
    assert cambios == ["MODIFICADA V02.TemperaturaColocacion_C: 18.0 -> 15.0"]



def test_cambio_de_geometria_se_detecta(ruta_ifc, tmp_path):
    """Simula mover el pedestal (V05) 0,5 m hacia arriba en Bonsai:
    la comparación debe informar el cambio de cotas aunque ninguna
    propiedad haya cambiado."""
    import ifcopenshell.api.geometry
    import ifcopenshell.util.placement

    modelo = bim.abrir(ruta_ifc)
    pedestal = bim.elementos_vertida(modelo)[-1]
    matriz = ifcopenshell.util.placement.get_local_placement(pedestal.ObjectPlacement)
    matriz[2, 3] += 0.5
    ifcopenshell.api.geometry.edit_object_placement(modelo, product=pedestal, matrix=matriz)
    ruta_v2 = tmp_path / "fundacion_movida.ifc"
    modelo.write(str(ruta_v2))

    cambios = bim.comparar_vertidas(ruta_ifc, ruta_v2)
    assert len(cambios) == 2
    assert all(c.startswith("MODIFICADA V05.Geometria.") for c in cambios)


def test_globalid_estable_al_regenerar(entradas):
    """Regenerar el IFC desde config/ entrega las mismas identidades
    (requisito para poder comparar versiones y para la trazabilidad)."""
    ids_1 = [v["GlobalId"] for v in bim.leer_vertidas(bim.crear_modelo_fundacion(*entradas))]
    ids_2 = [v["GlobalId"] for v in bim.leer_vertidas(bim.crear_modelo_fundacion(*entradas))]
    assert ids_1 == ids_2
    assert len(set(ids_1)) == 5  # y son distintas entre sí


def test_plan_con_hueco_es_rechazado(entradas):
    """Si las vertidas no llenan la fundación, se avisa (no se inventa nada)."""
    caso, plan = entradas
    plan_malo = {**plan, "vertidas": [dict(v) for v in plan["vertidas"]]}
    plan_malo["vertidas"][2]["z_inferior_m"] = 1.20  # deja un hueco 1.10-1.20
    with pytest.raises(ValueError, match="hueco o traslape"):
        bim.crear_modelo_fundacion(caso, plan_malo)
