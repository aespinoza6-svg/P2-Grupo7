"""Suelo bajo la fundacion: un bloque de celdas con capacidad termica.

Por que no basta una base adiabatica: el suelo absorbe calor del hormigon caliente y
tiene inercia. Se modela como una "vertida" mas, activa desde t = 0, sin hidratacion:
asi el mismo solver lo resuelve sin codigo especial.

Bordes del bloque (ver mcplan.bordes.caras_ambientales):
- fondo: temperatura fija (a algunos metros casi no cambia en 14 dias);
- costados: adiabaticos (campo lejano, el bloque es lo bastante ancho);
- cara superior alrededor de la fundacion: aire y sol, como el hormigon.
Valores por defecto: tipicos de un suelo limoso humedo, sin ensayos del sitio.
"""
import numpy as np
from mcplan.hidratacion import HidratacionExponencial
from mcplan.voxelizacion import GrillaVoxeles

SUELO_ID = 90


def agregar_suelo(grilla: GrillaVoxeles, *, espesor_m, margen_m, T_inicial_C, mezcla='SUELO'):
    """Nueva grilla con capas de suelo bajo z = 0 y un margen lateral alrededor."""
    if espesor_m <= 0 or margen_m < 0:
        raise ValueError('Espesor de suelo > 0 y margen >= 0.')
    if SUELO_ID in grilla.vertidas:
        raise ValueError('La grilla ya tiene suelo.')
    dx, dy, dz = grilla.paso
    mx, my, nzs = int(np.ceil(margen_m / dx - 1e-9)), int(np.ceil(margen_m / dy - 1e-9)), int(np.ceil(espesor_m / dz - 1e-9))
    nx, ny, nz = grilla.vertida_id.shape
    ids = np.zeros((nx + 2 * mx, ny + 2 * my, nz + nzs), dtype=np.int32)
    ids[mx:mx + nx, my:my + ny, nzs:] = grilla.vertida_id
    ids[:, :, :nzs] = SUELO_ID
    origen = np.asarray(grilla.origen, dtype=float) - np.array([mx * dx, my * dy, nzs * dz])
    volumen = float(np.count_nonzero(ids == SUELO_ID) * grilla.volumen_celda)
    vertidas = {**grilla.vertidas, SUELO_ID: dict(Nombre='Suelo', VertidaID=SUELO_ID, Inicio_h=0.0,
                                                  TemperaturaColocacion_C=float(T_inicial_C), TipoHormigon=mezcla,
                                                  volumen_ifc_m3=volumen)}
    return GrillaVoxeles(origen, grilla.paso, ids, vertidas, grilla.centros_compartidos)


def modelo_suelo(rho_kg_m3, cp_J_kgK):
    """'Hidratacion' nula con la capacidad del suelo (el solver exige que coincidan)."""
    return HidratacionExponencial(elevacion_final_C=0.0, tau_s=1.0, rho_kg_m3=rho_kg_m3, cp_J_kgK=cp_J_kgK)


def nucleo_superficie_por_vertida(T, activa, sistema, excluir=(SUELO_ID,)):
    """Para cada vertida activa: su nucleo (celda mas caliente) y su superficie (celda
    expuesta mas fria). Devuelve {vertida: (T_nucleo, T_superficie, i_nucleo, i_superficie)},
    con i = indice plano de la celda en la grilla (np.unravel_index lo vuelve (i, j, k)).

    Superficie: celda activa con al menos una cara al aire (x±, y±, z+) en ese instante.
    Es una aproximacion por celdas: el criterio del curso mide a 25 mm de la cara.
    """
    ids = sistema.programa.grilla.vertida_id
    superficie = np.zeros(ids.shape, dtype=bool)
    for cara in ('x-', 'x+', 'y-', 'y+', 'z+'):
        superficie |= sistema.expuestas[cara]
    plano = np.ravel(T)
    salida = {}
    for vid in sistema.programa.vertidas:
        if vid in excluir:
            continue
        m = np.flatnonzero(activa & (ids == vid))
        s = np.flatnonzero(activa & (ids == vid) & superficie)
        if m.size and s.size:
            i_n, i_s = m[np.argmax(plano[m])], s[np.argmin(plano[s])]
            salida[vid] = (float(plano[i_n]), float(plano[i_s]), int(i_n), int(i_s))
    return salida


def diferencial_por_vertida(T, activa, sistema, excluir=(SUELO_ID,)):
    """Nucleo - superficie dentro de cada vertida activa; devuelve (peor ΔT, vertida)."""
    peor, cual = np.nan, None
    for vid, (tn, ts, _, _) in nucleo_superficie_por_vertida(T, activa, sistema, excluir).items():
        if not tn - ts <= peor:
            peor, cual = tn - ts, vid
    return peor, cual
