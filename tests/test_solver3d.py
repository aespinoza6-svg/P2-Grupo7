"""
Tests de verificación del solver 2D/3D (mcplan.solver3d).

Siguen la "batería mínima de tests" de la guía del curso
(Antecedentes/guias/DIFUSION_DE_CALOR_1D_2D_3D.md, sección 33):

    2D exacto, 3D exacto, energía adiabática, fuente uniforme,
    Neumann, Robin  (+ estabilidad, interfaz entre materiales y simetría).

Correr con:   pytest -v
"""

import numpy as np
import pytest

from mcplan import analiticas as exacta
from mcplan import solver1d
from mcplan.solver3d import (
    Dirichlet,
    Malla3D,
    Material,
    Neumann,
    Robin,
    paso_tiempo_maximo,
    resolver,
)

DIDACTICO = Material("didactico", k=1.0, rho=1000.0, cp=10.0)   # alpha = 1e-4 m²/s
HORMIGON = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)
SUELO = Material("suelo", k=1.0, rho=1800.0, cp=1000.0)

TODAS_DIRICHLET_20 = {cara: Dirichlet(20.0) for cara in ("x-", "x+", "y-", "y+", "z-", "z+")}


def malla_cubo(n, dims=3, material=DIDACTICO):
    """Cuadrado (dims=2, nz=1) o cubo de lado 1 m con n celdas por lado."""
    forma = (n, n, n) if dims == 3 else (n, n, 1)
    return Malla3D.uniforme((1.0, 1.0, 1.0), forma, material)


def campo_seno(malla, T_borde=20.0, amplitud=40.0, dims=3):
    """Condición inicial del test exacto: seno en cada dirección activa."""
    x, y, z = malla.centros
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    coords, largos = ((X, Y, Z), (1, 1, 1)) if dims == 3 else ((X, Y), (1, 1))
    return coords, largos, exacta.seno_dirichlet_nd(coords, 0.0, largos, 1e-4, T_borde, amplitud)


def error_seno(n, dims, t_final=300.0):
    """Error máximo contra la solución exacta, con Fo fijo (dt ∝ dx²)."""
    malla = malla_cubo(n, dims)
    coords, largos, T0 = campo_seno(malla, dims=dims)
    bordes = {c: b for c, b in TODAS_DIRICHLET_20.items() if dims == 3 or c[0] != "z"}
    # Fo = 0.1: las esquinas con 3 caras Dirichlet exigen Fo <= 1/9 en 3D.
    dt = 0.1 * malla.d[0] ** 2 / DIDACTICO.alpha
    res = resolver(malla, T0, bordes, t_final, dt=dt, n_guardados=1)
    T_ex = exacta.seno_dirichlet_nd(coords, res.t[-1], largos, DIDACTICO.alpha, 20.0, 40.0)
    return np.max(np.abs(res.T[-1] - T_ex))


# ---------------------------------------------------------------------------
# Estabilidad
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims, limite", [(2, 1 / 4), (3, 1 / 6)])
def test_paso_maximo_es_fourier_un_cuarto_en_2d_y_un_sexto_en_3d(dims, limite):
    """Celda interior: 4 vecinas en 2D, 6 en 3D -> Fo <= 1/4 y 1/6."""
    malla = malla_cubo(10, dims)
    dt_max = paso_tiempo_maximo(malla, {})  # todo aislado
    assert DIDACTICO.alpha * dt_max / malla.d[0] ** 2 == pytest.approx(limite)


def test_dt_inestable_lanza_error():
    malla = malla_cubo(6)
    dt_max = paso_tiempo_maximo(malla, TODAS_DIRICHLET_20)
    with pytest.raises(ValueError, match="supera el paso máximo estable"):
        resolver(malla, 20.0, TODAS_DIRICHLET_20, 100 * dt_max, dt=1.5 * dt_max)


def test_cara_mal_escrita_lanza_error():
    """Un error de tipeo ("z_sup") no debe quedar en silencio como cara aislada."""
    with pytest.raises(KeyError, match="Caras desconocidas"):
        resolver(malla_cubo(4), 20.0, {"z_sup": Dirichlet(0.0)}, 10.0)


# ---------------------------------------------------------------------------
# Tests 2 y 3 de la guía: soluciones exactas
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims", [2, 3])
def test_solucion_exacta_seno(dims):
    """El seno decae como exp(-alpha·pi²·(1/Lx² + 1/Ly² [+ 1/Lz²])·t):
    comprueba que el stencil 2D/3D está bien armado."""
    assert error_seno(20, dims) < 0.1  # en una amplitud de 40 °C


@pytest.mark.parametrize("dims", [2, 3])
def test_convergencia_orden_2(dims):
    """Al duplicar las celdas (y dt ∝ dx²) el error baja ~4 veces: orden 2."""
    e_grueso, e_fino = error_seno(16, dims), error_seno(32, dims)
    orden = np.log2(e_grueso / e_fino)
    assert orden == pytest.approx(2.0, abs=0.1)


# ---------------------------------------------------------------------------
# Tests 4, 5 y 6 de la guía: energía, fuente y Neumann
# ---------------------------------------------------------------------------

def test_conservacion_energia_bordes_adiabaticos():
    """Todo aislado y sin fuente: la energía no cambia y la temperatura
    tiende al promedio inicial."""
    malla = malla_cubo(8)
    T0 = 20.0 + 30.0 * np.random.default_rng(0).random(malla.forma)
    res = resolver(malla, T0, {}, t_final=5000.0)
    E = res.energia()
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-12
    assert np.ptp(res.T[-1]) < 0.1
    assert res.T[-1].mean() == pytest.approx(T0.mean())


def test_fuente_uniforme_adiabatica_aumento_lineal():
    """dT/dt = Q/(rho·cp) en TODAS las celdas: la versión 3D del aumento
    adiabático del hormigón."""
    malla = malla_cubo(5, material=HORMIGON)
    Q = 1000.0  # W/m³
    res = resolver(malla, 15.0, {}, t_final=3 * 86400.0, fuente=Q, dt=3600.0)
    T_ex = exacta.adiabatico_fuente_uniforme(res.t, 15.0, Q, HORMIGON.rho, HORMIGON.cp)
    assert np.allclose(res.T, T_ex[:, None, None, None])


def test_neumann_flujo_entrante_suma_energia():
    """Entran q = 50 W/m² por la cara superior (z+) y el resto está aislado:
    la energía crece exactamente q·A·t."""
    malla = malla_cubo(6)
    q, A = 50.0, 1.0 * 1.0
    res = resolver(malla, 20.0, {"z+": Neumann(q)}, t_final=2000.0)
    assert np.allclose(res.energia() - res.energia()[0], q * A * res.t)
    # y el calor se acumula arriba: la capa superior está más caliente que la inferior
    assert res.T[-1][:, :, -1].mean() > res.T[-1][:, :, 0].mean()


# ---------------------------------------------------------------------------
# Test 7 de la guía: Robin
# ---------------------------------------------------------------------------

def test_robin_bloque_caliente_se_enfria_hacia_el_aire():
    """Bloque a 60 °C con aire a 10 °C en todas las caras: el calor sale
    (flujo con el signo de T_aire - T_superficie), ninguna celda baja de la
    temperatura del aire y el núcleo queda más caliente que las esquinas."""
    malla = malla_cubo(8)
    bordes = {c: Robin(h=10.0, T_inf=10.0) for c in TODAS_DIRICHLET_20}
    res = resolver(malla, 60.0, bordes, t_final=3000.0)
    assert np.all(np.diff(res.energia()) < 0)
    assert res.T.min() >= 10.0
    assert res.T[-1][4, 4, 4] > res.T[-1][0, 0, 0]


def test_3d_con_caras_laterales_aisladas_es_igual_al_1d():
    """Una barra 3D con las caras y, z aisladas no tiene por qué fluir el calor
    de lado: debe dar EXACTAMENTE lo mismo que solver1d (Dirichlet + Robin con
    moldaje). Conecta este solver con todo lo verificado en la semana 1."""
    n, L = 20, 1.0
    izq, der = Dirichlet(60.0), Robin(h=15.0, T_inf=5.0, R_extra=0.05)
    res1 = solver1d.resolver(solver1d.Malla1D.uniforme(L, n, DIDACTICO), 20.0, izq, der,
                             t_final=3000.0, dt=5.0, n_guardados=3)
    malla = Malla3D.uniforme((L, 0.3, 0.2), (n, 3, 2), DIDACTICO)
    res3 = resolver(malla, 20.0, {"x-": izq, "x+": der}, t_final=3000.0, dt=5.0, n_guardados=3)
    for j in range(3):
        for l in range(2):
            assert np.allclose(res3.T[:, :, j, l], res1.T, rtol=0, atol=1e-10)


def test_estacionario_robin_con_moldaje():
    """Estado estacionario de la misma barra contra la fórmula de resistencias
    en serie (estacionario_dirichlet_robin)."""
    n, L = 20, 1.0
    malla = Malla3D.uniforme((L, 0.2, 0.2), (n, 2, 2), DIDACTICO)
    bordes = {"x-": Dirichlet(60.0), "x+": Robin(h=15.0, T_inf=5.0, R_extra=0.05)}
    res = resolver(malla, 20.0, bordes, t_final=40000.0, n_guardados=1)
    T_ex = exacta.estacionario_dirichlet_robin(malla.centros[0], L, 1.0, 60.0, 15.0, 5.0, 0.05)
    assert np.allclose(res.T[-1][:, 0, 0], T_ex, atol=1e-3)


# ---------------------------------------------------------------------------
# Interfaces y simetría
# ---------------------------------------------------------------------------

def test_interfaz_hormigon_sobre_suelo():
    """Hormigón sobre suelo (dos materiales en z), arriba 40 °C y abajo 10 °C:
    el estacionario es lineal por tramos con el quiebre en la interfaz
    (resistencias en serie). Verifica la conductancia entre materiales."""
    n = 20
    indice = np.zeros((2, 2, n), dtype=int)
    indice[:, :, n // 2:] = 1  # mitad inferior suelo (0), superior hormigón (1)
    malla = Malla3D.por_materiales((1.0, 1.0, 2.0), indice, [SUELO, HORMIGON])
    bordes = {"z-": Dirichlet(10.0), "z+": Dirichlet(40.0)}
    res = resolver(malla, 20.0, bordes, t_final=200 * 86400.0, n_guardados=1)
    T_ex = exacta.estacionario_dos_capas(malla.centros[2], 1.0, SUELO.k, 1.0, HORMIGON.k,
                                         10.0, 40.0)
    assert np.allclose(res.T[-1][0, 0, :], T_ex, atol=1e-3)


def test_simetria_del_campo():
    """Problema simétrico (cubo, mismas condiciones en todas las caras) ->
    campo simétrico al intercambiar x, y, z. Detecta errores de índices."""
    malla = malla_cubo(7)
    bordes = {c: Robin(h=8.0, T_inf=0.0) for c in TODAS_DIRICHLET_20}
    T = resolver(malla, 50.0, bordes, t_final=2000.0, n_guardados=1).T[-1]
    assert np.allclose(T, T.transpose(1, 0, 2))
    assert np.allclose(T, T.transpose(2, 1, 0))
    assert np.allclose(T, T[::-1, :, :])


# ---------------------------------------------------------------------------
# Método de Schmidt (Anexo 4 del ICH)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dims", [2, 3])
def test_en_el_limite_es_el_metodo_de_schmidt(dims):
    """Con Fo justo en el límite (1/4 en 2D, 1/6 en 3D) el término de la propia
    celda desaparece: la nueva temperatura es el PROMEDIO de sus vecinas.
    Esa es la regla del método de Schmidt, así que nuestro solver la contiene."""
    malla = malla_cubo(6, dims)
    T0 = 20.0 + 30.0 * np.random.default_rng(1).random(malla.forma)
    dt = paso_tiempo_maximo(malla, {})  # Fo = 1/4 (2D) o 1/6 (3D)
    T1 = resolver(malla, T0, {}, t_final=dt, dt=dt, n_guardados=1).T[-1]

    i = j = 2
    l = 2 if dims == 3 else 0
    vecinas = [T0[i - 1, j, l], T0[i + 1, j, l], T0[i, j - 1, l], T0[i, j + 1, l]]
    if dims == 3:
        vecinas += [T0[i, j, l - 1], T0[i, j, l + 1]]
    assert T1[i, j, l] == pytest.approx(np.mean(vecinas))
