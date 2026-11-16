"""E13: riesgo termico de un plan de vertidas bajo escenarios ponderados.

Piezas, en el orden en que se usan:
1. cronograma: con una capacidad de colocacion r [m3/h], cada vertida dura d = V/r.
   D(u) = fin de la ultima vertida - t0. La vertida se ACTIVA de una vez en el punto
   medio de su colocacion (aproximacion termica declarada); el cronograma cuenta su
   duracion completa.
2. temperatura_25mm: estimacion de T a 25 mm bajo cada cara exterior existente (S25).
3. metricas MT y MΔ del enunciado, sobre el hormigon colocado (sin suelo ni celdas futuras).
4. elegir_trayectorias: reduce el ensemble a pocas historias completas con pesos.
5. probabilidades: enumeracion ponderada exacta del modelo discreto.
"""
from __future__ import annotations
import copy
import numpy as np
from mcplan.solver3d import CARAS
from mcplan.solver_vertidas import _es_campo
from mcplan.solver1d import Neumann
from mcplan.voxelizacion import GrillaVoxeles

SUELO_ID = 90
PROFUNDIDAD_S25_M = 0.025
LIMITE_T_C, LIMITE_DT_C = 70.0, 20.0


# ---------------------------------------------------------------------------
# 1. Cronograma y plan
# ---------------------------------------------------------------------------
def cronograma(volumenes_m3, capacidad_m3_h, espera_h=0.0, unica=False):
    """Inicio, fin y activacion [h desde t0] de cada vertida, y D [h].

    volumenes_m3: {vertida_id: V} en el orden de colocacion.
    espera_h: espera entre el FIN de una vertida y el inicio de la siguiente (una bomba,
        sin simultaneidad). unica=True: todo el volumen es una sola colocacion continua.
    """
    if capacidad_m3_h <= 0 or espera_h < 0:
        raise ValueError('Capacidad > 0 y espera >= 0.')
    ids = list(volumenes_m3)
    if unica:
        d = sum(volumenes_m3.values()) / capacidad_m3_h
        filas = {v: dict(inicio_h=0.0, fin_h=d, duracion_h=d, activacion_h=d / 2) for v in ids}
        return filas, d
    filas, t = {}, 0.0
    for j, v in enumerate(ids):
        d = volumenes_m3[v] / capacidad_m3_h
        filas[v] = dict(inicio_h=t, fin_h=t + d, duracion_h=d, activacion_h=t + d / 2)
        t += d + (espera_h if j < len(ids) - 1 else 0.0)
    return filas, filas[ids[-1]]['fin_h']


def grilla_con_plan(grilla, crono, temperatura_C, mezcla='M1'):
    """Copia de la grilla del IFC con los horarios y la temperatura del plan."""
    vertidas = copy.deepcopy(grilla.vertidas)
    for v, f in crono.items():
        vertidas[v].update(Inicio_h=float(f['activacion_h']), TemperaturaColocacion_C=float(temperatura_C), TipoHormigon=mezcla)
    return GrillaVoxeles(grilla.origen, grilla.paso, grilla.vertida_id, vertidas, grilla.centros_compartidos)


# ---------------------------------------------------------------------------
# 2-3. Superficie a 25 mm y metricas del enunciado
# ---------------------------------------------------------------------------
def temperatura_25mm(T, sistema, t_s, profundidad_m=PROFUNDIDAD_S25_M, excluir=(SUELO_ID,)):
    """T estimada a 25 mm bajo cada cara exterior del hormigon colocado.

    La celda solo conoce su centro (a d/2 de la cara). Con el flujo que sale por esa
    cara, q, se reconstruye la temperatura de la cara: T_cara = T_c - q (d/2)/k, y se
    interpola en linea recta entre la cara y el centro. Si d/2 <= 25 mm se usa T_c.
    Devuelve (T_minima, indice plano de la celda, cara).
    """
    ids = sistema.programa.grilla.vertida_id
    mejor = (np.inf, None, None)
    for cara in ('x-', 'x+', 'y-', 'y+', 'z+'):
        e = CARAS[cara][0]
        mask = sistema.expuestas[cara] & ~np.isin(ids, excluir)
        if not mask.any():
            continue
        Tc = T[mask]
        borde = sistema.bordes[cara]
        if _es_campo(borde):
            G, Tref, Q = borde.coeficientes(t_s, cara, mask, sistema)
            q_sale = (G * (Tc - Tref) - Q) / sistema.areas[e]        # W/m2 que salen
        elif isinstance(borde, Neumann):
            q_sale = np.full(Tc.shape, -float(borde.q) if np.isscalar(borde.q) else 0.0)
        else:
            raise ValueError(f'S25: tipo de borde no previsto en {cara}.')
        medio = sistema.d[e] / 2
        T_cara = Tc - q_sale * medio / sistema.k[mask]
        T25 = T_cara + (Tc - T_cara) * min(profundidad_m / medio, 1.0)
        i = int(np.argmin(T25))
        if T25[i] < mejor[0]:
            mejor = (float(T25[i]), int(np.flatnonzero(mask)[i]), cara)
    return mejor


def medir(T, activa, sistema, t_s, excluir=(SUELO_ID,)):
    """MΔ(t) del enunciado: max T en el hormigon colocado - min T en S25 (y version por celdas)."""
    ids = sistema.programa.grilla.vertida_id
    hormigon = activa & ~np.isin(ids, excluir)
    if not hormigon.any():
        return None
    plano = np.ravel(T)
    i_max = int(np.flatnonzero(hormigon)[np.argmax(T[hormigon])])
    t25, i25, cara = temperatura_25mm(T, sistema, t_s, excluir=excluir)
    sup = np.zeros(T.shape, bool)
    for c in ('x-', 'x+', 'y-', 'y+', 'z+'):
        sup |= sistema.expuestas[c]
    sup &= hormigon
    return dict(T_max=float(plano[i_max]), i_max=i_max, T25_min=t25, i25=i25, cara25=cara,
                T_celda_min=float(T[sup].min()) if sup.any() else np.nan)


# ---------------------------------------------------------------------------
# 4. Pocas historias climaticas completas desde el ensemble
# ---------------------------------------------------------------------------
def elegir_trayectorias(series, n=3):
    """series: {miembro: SerieClima}. Ordena por la temperatura media del aire en las
    horas con pronostico, parte en n grupos de igual tamano (frio, central, calido) y
    toma de cada grupo el miembro mas cercano a su mediana. Peso = fraccion de miembros
    del grupo. Devuelve tambien el miembro mas frio y el mas calido (sensibilidad)."""
    indice = {m: float(np.mean(s.T_aire_C[s.fuente == 'pronostico'])) for m, s in series.items()}
    orden = sorted(indice, key=indice.get)
    grupos = np.array_split(np.array(orden), n)
    nombres = ['frio', 'central', 'calido'] if n == 3 else [f'g{j + 1}' for j in range(n)]
    elegidas = []
    for nombre, grupo in zip(nombres, grupos):
        mediana = float(np.median([indice[m] for m in grupo]))
        rep = min(grupo, key=lambda m: abs(indice[m] - mediana))
        elegidas.append(dict(nombre=nombre, miembro=str(rep), peso=len(grupo) / len(orden),
                             indice_C=indice[rep], grupo=[str(m) for m in grupo],
                             rango_C=[indice[grupo[0]], indice[grupo[-1]]]))
    return elegidas, dict(frio_extremo=orden[0], calido_extremo=orden[-1], indice_C=indice)


# ---------------------------------------------------------------------------
# 5. Enumeracion ponderada
# ---------------------------------------------------------------------------
def probabilidades(filas, p_material, q_clima):
    """filas: dicts con curva, clima, MT, MD. Pesos w = p_m q_k (independencia material-clima).

    Escenarios con peso cero no cuentan. Devuelve la conjunta y las marginales.
    """
    P = dict(conjunta=0.0, T=0.0, dT=0.0, peso_total=0.0, fallas=[])
    for f in filas:
        w = p_material.get(f['curva'], 0.0) * q_clima.get(f['clima'], 0.0)
        if w == 0:
            continue
        falla_T, falla_dT = f['MT'] > LIMITE_T_C, f['MD'] > LIMITE_DT_C
        P['peso_total'] += w
        P['T'] += w * falla_T
        P['dT'] += w * falla_dT
        if falla_T or falla_dT:
            P['conjunta'] += w
            P['fallas'].append(f"{f['curva']}|{f['clima']}")
    if not np.isclose(P['peso_total'], 1.0):
        raise ValueError(f"Los pesos suman {P['peso_total']:.4f}, no 1: faltan escenarios.")
    return P
