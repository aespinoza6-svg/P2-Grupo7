"""
Modelo BIM de la fundación en formato IFC (con IfcOpenShell).

Arquitectura del curso:

    Bonsai / Blender  ->  IFC  ->  IfcOpenShell  ->  geometría + propiedades  ->  ...  ->  solver

Este módulo hace dos cosas:

1) CREAR el IFC de la fundación a partir de los archivos de config/
   (sitio, geometría, plan de vertidas). Así el modelo es reproducible:
   si cambia una cota o una fecha, se cambia el YAML y se regenera.
   El archivo resultante se abre y edita en Bonsai como cualquier IFC.

2) LEER cualquier IFC (generado aquí o editado en Bonsai) y extraer las
   vertidas con sus propiedades, su geometría y las diferencias entre dos
   versiones del modelo.

Nunca se importa Blender/Bonsai: el IFC es el único puente. Por eso el
cálculo se puede correr en un computador sin Blender (p. ej. el cluster).

Estructura del IFC creado (IFC4):

    IfcProject  "P2 Grupo 7 - Fundación aerogenerador"
    └── IfcSite  "Renaico - P.E. Vergara"          (lat/lon + pset MC_Sitio)
        └── IfcBuilding  "Aerogenerador"
            └── IfcBuildingStorey  "Nivel fundación"
                ├── IfcFooting  "Vertida_01"         (pset UANDES_HormigonMasivo)
                ├── ...
                └── IfcFooting  "V05 - pedestal"

Cada vertida es un IfcFooting propio: tiene identidad (GlobalId) y su
property set UANDES_HormigonMasivo con los datos de planificación y construcción.
"""

from __future__ import annotations

import uuid
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import yaml

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.api.aggregate
import ifcopenshell.api.context
import ifcopenshell.api.geometry
import ifcopenshell.api.material
import ifcopenshell.api.project
import ifcopenshell.api.pset
import ifcopenshell.api.root
import ifcopenshell.api.spatial
import ifcopenshell.api.unit
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.geolocation
import ifcopenshell.util.shape

# Property Set docente indicado en GUIA_BONSAI_FUNDACION.md. No usa el prefijo
# Pset_, que queda reservado para conjuntos estandarizados por buildingSMART.
PSET_VERTIDA = "UANDES_HormigonMasivo"
PSET_SITIO = "MC_Sitio"

# Campos mínimos que cada vertida debe tener (enunciado, sección 7).
PROPIEDADES_VERTIDA = (
    "VertidaID",
    "Elemento",
    "Forma",
    "TipoHormigon",
    "Inicio_h",
    "InicioPlanificado",
    "TemperaturaColocacion_C",
    "Moldaje",
    "Aislacion",
    "RetiroMoldaje_h",
)


# ---------------------------------------------------------------------------
# Lectura de la configuración
# ---------------------------------------------------------------------------

def cargar_yaml(ruta: str | Path) -> dict:
    """Lee un archivo YAML de config/ y lo devuelve como diccionario."""
    with open(ruta, encoding="utf-8") as f:
        return yaml.safe_load(f)


def validar_plan(geometria: dict, plan: dict) -> None:
    """Revisa continuidad vertical, formas permitidas y cotas de cada zona."""
    vertidas = sorted(plan["vertidas"], key=lambda v: v["z_inferior_m"])
    if not vertidas:
        raise ValueError("El plan no contiene vertidas.")
    ids = [v["id"] for v in vertidas]
    if len(ids) != len(set(ids)):
        raise ValueError("El plan contiene IDs de vertida duplicados.")

    z = 0.0
    for v in vertidas:
        if not np.isclose(v["z_inferior_m"], z):
            raise ValueError(
                f"{v['id']}: empieza en z={v['z_inferior_m']} m pero debería "
                f"empezar en z={z:.3f} m (hueco o traslape)."
            )
        if v["z_superior_m"] <= v["z_inferior_m"]:
            raise ValueError(f"{v['id']}: z_superior debe ser mayor que z_inferior.")
        if v.get("forma") not in {"prisma_octagonal", "tronco_octagonal", "cilindro"}:
            raise ValueError(f"{v['id']}: forma IFC desconocida {v.get('forma')!r}.")
        z0, z1, forma = v["z_inferior_m"], v["z_superior_m"], v["forma"]
        h_borde, h_centro, h_total = (
            geometria["espesor_borde_m"], geometria["espesor_centro_m"], geometria["altura_total_m"]
        )
        if forma == "prisma_octagonal" and (z0 < -1e-9 or z1 > h_borde + 1e-9):
            raise ValueError(f"{v['id']}: el prisma octagonal debe estar dentro de z=0..{h_borde} m.")
        if forma == "tronco_octagonal" and (z0 < h_borde - 1e-9 or z1 > h_centro + 1e-9):
            raise ValueError(f"{v['id']}: la transición octagonal debe estar dentro de z={h_borde}..{h_centro} m.")
        if forma == "cilindro" and (z0 < h_centro - 1e-9 or z1 > h_total + 1e-9):
            raise ValueError(f"{v['id']}: el cilindro debe estar dentro de z={h_centro}..{h_total} m.")
        z = v["z_superior_m"]

    if not np.isclose(z, geometria["altura_total_m"]):
        raise ValueError(
            f"Las vertidas terminan en z={z} m, no en {geometria['altura_total_m']:.3f} m."
        )


# ---------------------------------------------------------------------------
# Creación del modelo IFC
# ---------------------------------------------------------------------------

# Espacio de nombres para generar GlobalId deterministas (ver _guid_estable).
_ESPACIO_GUID = uuid.uuid5(uuid.NAMESPACE_URL, "https://github.com/aespinoza6-svg/P2-Grupo7")


def _guid_estable(clave: str) -> str:
    """GlobalId IFC que depende solo de 'clave' (p. ej. "vertida/V01").

    IfcOpenShell inventa un GlobalId aleatorio para cada objeto nuevo. Si lo
    dejáramos así, cada vez que se regenera el IFC la "misma" vertida tendría
    otra identidad y no se podrían comparar versiones. Con uuid5 la misma
    clave produce siempre el mismo identificador.
    """
    return ifcopenshell.guid.compress(uuid.uuid5(_ESPACIO_GUID, clave).hex)


def _vertida_id_numero(identificador) -> int:
    """Convierte V01, Vertida_01 o 1 al VertidaID entero usado por la guía."""
    digitos = "".join(c for c in str(identificador) if c.isdigit())
    if not digitos:
        raise ValueError(f"ID de vertida inválido: {identificador!r}")
    return int(digitos)


def _etiqueta_vertida(identificador) -> str:
    """Etiqueta legible V01 para mostrar IDs IFC enteros."""
    return f"V{_vertida_id_numero(identificador):02d}"


def _fijar_guids(modelo: ifcopenshell.file, claves: dict) -> None:
    """Asigna GlobalId deterministas a todo objeto IfcRoot del modelo.

    claves: {id de la entidad: clave legible} para los objetos importantes; el resto
    (relaciones, property sets) usa su clase y orden de creación, que también
    es siempre el mismo porque el script es determinista.
    """
    for entidad in modelo.by_type("IfcRoot"):
        clave = claves.get(entidad.id()) or f"{entidad.is_a()}#{entidad.id()}"
        entidad.GlobalId = _guid_estable(clave)


def _angulo_compuesto(grados_decimales: float) -> list[int]:
    """Convierte grados decimales al formato IFC (grados, minutos, segundos,
    millonésimas de segundo), p. ej. -37.6921 -> [-37, -41, -31, -560000]."""
    return list(ifcopenshell.util.geolocation.dd2dms(grados_decimales, use_us=True))


def _radio_octagono_en_z(geometria: dict, z: float) -> float:
    """Radio circunscrito del octágono: constante en el borde y lineal en la transición."""
    n = int(geometria["lados_planta"])
    apotema_ext = float(geometria["distancia_caras_exteriores_m"]) / 2
    apotema_int = float(geometria["plataforma_distancia_caras_m"]) / 2
    r_ext = apotema_ext / np.cos(np.pi / n)
    r_int = apotema_int / np.cos(np.pi / n)
    h_borde, h_centro = geometria["espesor_borde_m"], geometria["espesor_centro_m"]
    if z <= h_borde:
        return float(r_ext)
    fraccion = (z - h_borde) / (h_centro - h_borde)
    return float(r_ext + fraccion * (r_int - r_ext))


def _crear_representacion_vertida(modelo, contexto, v: dict, geometria: dict):
    """Crea un sólido cerrado faceteado para prisma, tronco octagonal o cilindro."""
    z0, z1 = float(v["z_inferior_m"]), float(v["z_superior_m"])
    h = z1 - z0
    forma = v["forma"]
    if forma == "cilindro":
        n = int(geometria["lados_cilindro"])
        radio0 = radio1 = float(geometria["pedestal_diametro_m"]) / 2
        offset = 0.0
    else:
        n = int(geometria["lados_planta"])
        radio0 = _radio_octagono_en_z(geometria, z0)
        radio1 = _radio_octagono_en_z(geometria, z1)
        offset = np.pi / n  # caras alineadas con X/Y, como indica la guía

    vertices = []
    for radio, z_local in ((radio0, 0.0), (radio1, h)):
        vertices.extend(
            (radio * np.cos(2 * np.pi * i / n + offset),
             radio * np.sin(2 * np.pi * i / n + offset), z_local)
            for i in range(n)
        )
    caras = [list(reversed(range(n))), list(range(n, 2 * n))]
    caras.extend([[i, (i + 1) % n, (i + 1) % n + n, i + n] for i in range(n)])
    return ifcopenshell.api.geometry.add_mesh_representation(
        modelo, context=contexto, vertices=[vertices], faces=[caras], force_faceted_brep=True
    )


def crear_modelo_fundacion(caso: dict, plan: dict) -> ifcopenshell.file:
    """Construye el IFC de la fundación con una IfcFooting por vertida.

    caso: diccionario de config/caso_base.yaml (secciones sitio, geometria, hormigon).
    plan: diccionario de config/plan_vXX.yaml (las vertidas).
    """
    sitio, geo, h = caso["sitio"], caso["geometria"], caso["hormigon"]
    validar_plan(geo, plan)

    # 1) Archivo, proyecto y unidades (metros: así las cotas del YAML se
    #    escriben tal cual; IfcOpenShell por defecto usaría milímetros).
    modelo = ifcopenshell.api.project.create_file(version="IFC4")
    modelo.header.file_name.author = ("Grupo 7",)
    modelo.header.file_name.organization = ("UANDES - Metodos Computacionales",)
    proyecto = ifcopenshell.api.root.create_entity(
        modelo, ifc_class="IfcProject", name="P2 Grupo 7 - Fundación aerogenerador"
    )
    unidades = [
        ifcopenshell.api.unit.add_si_unit(modelo, unit_type="LENGTHUNIT"),
        ifcopenshell.api.unit.add_si_unit(modelo, unit_type="AREAUNIT"),
        ifcopenshell.api.unit.add_si_unit(modelo, unit_type="VOLUMEUNIT"),
    ]
    ifcopenshell.api.unit.assign_unit(modelo, units=unidades)

    # 2) Contexto geométrico: "Body" = la forma 3D del objeto.
    modelo3d = ifcopenshell.api.context.add_context(modelo, context_type="Model")
    cuerpo = ifcopenshell.api.context.add_context(
        modelo, context_type="Model", context_identifier="Body",
        target_view="MODEL_VIEW", parent=modelo3d,
    )

    # 3) Estructura espacial: proyecto > sitio > edificio > nivel.
    ifc_sitio = ifcopenshell.api.root.create_entity(
        modelo, ifc_class="IfcSite", name=f"{sitio['localidad']} - {sitio['parque_eolico']}"
    )
    ifc_sitio.RefLatitude = _angulo_compuesto(sitio["latitud"])
    ifc_sitio.RefLongitude = _angulo_compuesto(sitio["longitud"])
    pset = ifcopenshell.api.pset.add_pset(modelo, product=ifc_sitio, name=PSET_SITIO)
    ifcopenshell.api.pset.edit_pset(modelo, pset=pset, properties={
        "Grupo": sitio["grupo"],
        "Localidad": sitio["localidad"],
        "Region": sitio["region"],
        "ParqueEolico": sitio["parque_eolico"],
        "Latitud_deg": float(sitio["latitud"]),
        "Longitud_deg": float(sitio["longitud"]),
    })
    edificio = ifcopenshell.api.root.create_entity(modelo, ifc_class="IfcBuilding", name="Aerogenerador")
    nivel = ifcopenshell.api.root.create_entity(modelo, ifc_class="IfcBuildingStorey", name="Nivel fundación")
    ifcopenshell.api.aggregate.assign_object(modelo, products=[ifc_sitio], relating_object=proyecto)
    ifcopenshell.api.aggregate.assign_object(modelo, products=[edificio], relating_object=ifc_sitio)
    ifcopenshell.api.aggregate.assign_object(modelo, products=[nivel], relating_object=edificio)
    claves = {proyecto.id(): "proyecto", ifc_sitio.id(): "sitio",
              edificio.id(): "edificio", nivel.id(): "nivel"}

    # 4) Material con sus propiedades térmicas en los property sets ESTÁNDAR
    #    de IFC (Pset_MaterialCommon y Pset_MaterialThermal).
    hormigon =ifcopenshell.api.material.add_material(modelo, name=h["nombre"], category="concrete")
    ps = ifcopenshell.api.pset.add_pset(modelo, product=hormigon, name="Pset_MaterialCommon")
    ifcopenshell.api.pset.edit_pset(modelo, pset=ps, properties={"MassDensity": float(h["densidad_kg_m3"])})
    ps = ifcopenshell.api.pset.add_pset(modelo, product=hormigon, name="Pset_MaterialThermal")
    ifcopenshell.api.pset.edit_pset(modelo, pset=ps, properties={
        "ThermalConductivity": float(h["conductividad_W_mK"]),
        "SpecificHeatCapacity": float(h["calor_especifico_J_kgK"]),
    })

    # 5) Una IfcFooting por vertida: geometría + material + property set.
    inicio = datetime.fromisoformat(plan["fecha_inicio"])
    for v in plan["vertidas"]:
        zapata = ifcopenshell.api.root.create_entity(
            modelo, ifc_class="IfcFooting", predefined_type="PAD_FOOTING",
            name=f"Vertida_{_vertida_id_numero(v['id']):02d}",
        )
        ifcopenshell.api.spatial.assign_container(modelo, products=[zapata], relating_structure=nivel)
        claves[zapata.id()] = f"vertida/{v['id']}"

        # Cada sólido se construye en coordenadas locales y se posiciona en z_inferior.
        forma = _crear_representacion_vertida(modelo, cuerpo, v, geo)
        ifcopenshell.api.geometry.assign_representation(modelo, product=zapata, representation=forma)
        ubicacion = np.eye(4)
        ubicacion[2, 3] = v["z_inferior_m"]
        ifcopenshell.api.geometry.edit_object_placement(modelo, product=zapata, matrix=ubicacion)

        ifcopenshell.api.material.assign_material(modelo, products=[zapata], type="IfcMaterial", material=hormigon)

        pset = ifcopenshell.api.pset.add_pset(modelo, product=zapata, name=PSET_VERTIDA)
        ifcopenshell.api.pset.edit_pset(modelo, pset=pset, properties={
            "VertidaID": _vertida_id_numero(v["id"]),
            "Elemento": v["elemento"],
            "Forma": v["forma"],
            "TipoHormigon": v["tipo_hormigon"],
            "Inicio_h": float(v["inicio_h"]),
            "InicioPlanificado": (inicio + timedelta(hours=v["inicio_h"])).isoformat(),
            "TemperaturaColocacion_C": float(v["temperatura_colocacion_C"]),
            "Moldaje": v["moldaje"],
            "Aislacion": v["aislacion"],
            "RetiroMoldaje_h": float(v["retiro_moldaje_h"]),
        })

    # 6) Identidades estables: regenerar el IFC no cambia los GlobalId.
    _fijar_guids(modelo, claves)
    return modelo


# ---------------------------------------------------------------------------
# Lectura del modelo IFC (sirve para archivos generados aquí o en Bonsai)
# ---------------------------------------------------------------------------

def abrir(modelo_o_ruta: ifcopenshell.file | str | Path) -> ifcopenshell.file:
    """Acepta un modelo ya abierto o la ruta a un .ifc."""
    if isinstance(modelo_o_ruta, ifcopenshell.file):
        return modelo_o_ruta
    return ifcopenshell.open(str(modelo_o_ruta))


def elementos_vertida(modelo: ifcopenshell.file) -> list:
    """Todos los elementos que tienen el property set docente de vertidas
    (sin importar su clase IFC), ordenados por VertidaID.

    Recibe un modelo YA ABIERTO: los elementos devueltos solo son válidos
    mientras ese modelo exista.
    """
    encontrados = []
    for elemento in modelo.by_type("IfcElement"):
        psets = ifcopenshell.util.element.get_psets(elemento)
        if PSET_VERTIDA in psets:
            encontrados.append((psets[PSET_VERTIDA].get("VertidaID", ""), elemento))
    return [e for _, e in sorted(encontrados, key=lambda par: _vertida_id_numero(par[0]))]


def leer_vertidas(modelo) -> list[dict]:
    """Lista de vertidas como diccionarios: identidad IFC + propiedades.

    El GlobalId es la identidad permanente del objeto: no cambia aunque en
    Bonsai se edite su nombre, su forma o sus propiedades.
    """
    # Se abre aquí (y no solo dentro de elementos_vertida) para que el archivo
    # siga vivo mientras se leen sus elementos: si Python lo libera antes,
    # los elementos apuntan a memoria inválida y el programa se cae.
    modelo = abrir(modelo)
    vertidas = []
    for elemento in elementos_vertida(modelo):
        props = ifcopenshell.util.element.get_psets(elemento)[PSET_VERTIDA]
        props.pop("id", None)  # id interno de IfcOpenShell, no es un dato nuestro
        material = ifcopenshell.util.element.get_material(elemento)
        vertidas.append({
            "GlobalId": elemento.GlobalId,
            "ClaseIFC": elemento.is_a(),
            "Nombre": elemento.Name,
            "Material": material.Name if material else None,
            **props,
        })
    return vertidas


def propiedades_faltantes(vertida: dict) -> list[str]:
    """Campos del enunciado que le faltan a una vertida leída del IFC."""
    return [p for p in PROPIEDADES_VERTIDA if p not in vertida]


def geometria_vertidas(modelo) -> dict[str, dict]:
    """Volumen [m³] y caja envolvente [m] de cada vertida, calculados por
    IfcOpenShell a partir de la geometría del IFC (triangulada).

    Sirve para comprobar que lo que se modeló es lo que se quería
    (y, desde la semana 3, para voxelizar).
    """
    modelo = abrir(modelo)
    ajustes = ifcopenshell.geom.settings()
    ajustes.set("use-world-coords", True)  # coordenadas globales, no locales
    resultado = {}
    for elemento in elementos_vertida(modelo):
        vid = ifcopenshell.util.element.get_psets(elemento)[PSET_VERTIDA]["VertidaID"]
        forma = ifcopenshell.geom.create_shape(ajustes, elemento)
        vertices = np.array(forma.geometry.verts).reshape(-1, 3)
        resultado[vid] = {
            "GlobalId": elemento.GlobalId,
            "volumen_m3": float(ifcopenshell.util.shape.get_volume(forma.geometry)),
            "min_xyz_m": vertices.min(axis=0),
            "max_xyz_m": vertices.max(axis=0),
        }
    return resultado


def editar_propiedad(modelo, vertida_id: str, propiedad: str, valor) -> None:
    """Cambia una propiedad del pset de vertida (lo mismo que hacerlo a mano
    en el panel de propiedades de Bonsai). No guarda el archivo."""
    for elemento in elementos_vertida(modelo):
        psets = ifcopenshell.util.element.get_psets(elemento)
        if _vertida_id_numero(psets[PSET_VERTIDA].get("VertidaID")) == _vertida_id_numero(vertida_id):
            pset = modelo.by_id(psets[PSET_VERTIDA]["id"])
            ifcopenshell.api.pset.edit_pset(modelo, pset=pset, properties={propiedad: valor})
            return
    raise KeyError(f"No existe la vertida {vertida_id!r} en el modelo.")


def _resumen_geometria(g: dict) -> dict:
    """Volumen y caja envolvente redondeados a 1 mm / 0,01 m³, para que
    diferencias de redondeo numérico no se informen como cambios."""
    return {
        "Geometria.Volumen_m3": round(g["volumen_m3"], 2),
        "Geometria.Min_xyz_m": tuple(round(float(c), 3) for c in g["min_xyz_m"]),
        "Geometria.Max_xyz_m": tuple(round(float(c), 3) for c in g["max_xyz_m"]),
    }


def comparar_vertidas(antes, despues) -> list[str]:
    """Diferencias entre dos versiones del modelo, vertida por vertida.

    Las vertidas se emparejan por GlobalId (identidad), no por nombre:
    así un cambio de nombre se detecta como cambio y no como "borrada + nueva".

    Se comparan las PROPIEDADES (pset de vertida, nombre, material) y también
    la GEOMETRÍA: volumen y caja envolvente calculados desde el IFC. Así se
    detecta tanto "cambié la temperatura" como "moví o estiré la pieza en Bonsai".
    """
    antes, despues = abrir(antes), abrir(despues)
    a = {v["GlobalId"]: v for v in leer_vertidas(antes)}
    b = {v["GlobalId"]: v for v in leer_vertidas(despues)}
    geo_a = {g["GlobalId"]: g for g in geometria_vertidas(antes).values()}
    geo_b = {g["GlobalId"]: g for g in geometria_vertidas(despues).values()}
    for gid in a.keys() & b.keys():
        a[gid] = {**a[gid], **_resumen_geometria(geo_a[gid])}
        b[gid] = {**b[gid], **_resumen_geometria(geo_b[gid])}
    cambios = []
    for gid in sorted(a.keys() - b.keys()):
        cambios.append(f"ELIMINADA  {_etiqueta_vertida(a[gid].get('VertidaID'))} ({gid})")
    for gid in sorted(b.keys() - a.keys()):
        cambios.append(f"NUEVA      {_etiqueta_vertida(b[gid].get('VertidaID'))} ({gid})")
    for gid in sorted(a.keys() & b.keys()):
        for campo in sorted(a[gid].keys() | b[gid].keys()):
            va, vb = a[gid].get(campo), b[gid].get(campo)
            if va != vb:
                cambios.append(
                    f"MODIFICADA {_etiqueta_vertida(a[gid].get('VertidaID'))}.{campo}: {va!r} -> {vb!r}"
                )
    return cambios
