"""
Batería de verificación del lab E05 (jueves 29 de octubre).

Un test por cada caso del enunciado, en el mismo orden:

    1  solución 1D exacta            pytest -k exacta_1d
    2  solución 2D sinusoidal        pytest -k seno_2d
    3  solución 3D sinusoidal        pytest -k seno_3d
    4  conservación de energía       pytest -k energia
    5  fuente uniforme               pytest -k fuente
    6  Neumann                       pytest -k neumann
    7  Robin                         pytest -k robin

    Todos:              pytest tests/test_e05_bateria.py -v
    Viendo los números: pytest tests/test_e05_bateria.py -v -s

Para CAMBIAR PARÁMETROS: editar el bloque "PARÁMETROS DEL LAB", guardar y
volver a correr. Las tolerancias son relativas (a |A| o a la diferencia de
temperaturas), así que cambiar A, T0, los largos o el material no rompe nada.

La batería completa (más casos) está en test_solver1d.py y test_solver3d.py.
Explicación de cada caso: docs/2_guia_oral.md, sección E05.
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
    resolver,
)

# =============================================================================
# PARÁMETROS DEL LAB: cambiar aquí
# =============================================================================

Lx, Ly, Lz = 1.0, 1.0, 1.0       # largos del dominio [m] (en 2D no se usa Lz)
N = 20                           # celdas en CADA dirección: dx = Lx/N, dy = Ly/N, dz = Lz/N

k, rho, cp = 1.0, 1000.0, 10.0   # didáctico: alpha = 1e-4 m²/s (pasa en minutos)
# k, rho, cp = 2.0, 2400.0, 900.0  # hormigón: alpha = 9.3e-7 m²/s (pasa en horas/días)

T0 = 20.0                        # temperatura de las caras [°C]       (casos 1, 2, 3)
A = 40.0                         # amplitud inicial del seno [°C]      (casos 1, 2, 3, 4)

T_FINAL = None                   # [s]. None = lo que tarda el seno en bajar a ~60 %: 0.5/(alpha·λ)
DT = None                        # [s]. None = 90 % del máximo estable. Si es inestable, el solver avisa.

q0 = 1000.0                      # fuente volumétrica uniforme [W/m³]  (caso 5)
q_SUP = 50.0                     # flujo que ENTRA por la cara superior z+ [W/m²] (caso 6)
h, T_INF, R_EXTRA = 10.0, 5.0, 0.0   # convección [W/m²K], aire [°C], moldaje [m²K/W] (caso 7)

TOLERANCIA = 0.01                # error máximo aceptado: 1 % de |A|

# =============================================================================

MATERIAL = Material("lab", k=k, rho=rho, cp=cp)
CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")


def tiempo_final(largos):
    """T_FINAL, o el tiempo en que exp(-alpha·λ·t) = exp(-0.5) ≈ 0.61."""
    if T_FINAL is not None:
        return T_FINAL
    return 0.5 / (MATERIAL.alpha * exacta.lambda_seno(largos))


def amplitud(T, modo):
    """Cuánto del seno queda en el campo T: proyección de (T - T0) sobre el modo."""
    return float(np.sum((T - T0) * modo) / np.sum(modo * modo))


def comprobar_seno(nombre, centros, largos, T_solver, t):
    """Compara con la fórmula del enunciado y mide λ a partir del decaimiento."""
    mallas = np.meshgrid(*centros, indexing="ij")
    T_ex = exacta.seno_dirichlet_nd(mallas, t, largos, MATERIAL.alpha, T0, A)
    T_ex = T_ex.reshape(T_solver.shape)
    error = np.max(np.abs(T_solver - T_ex))
    lam = exacta.lambda_seno(largos)
    queda = np.exp(-MATERIAL.alpha * lam * t)          # fracción del seno que queda
    print(f"  [{nombre}] t = {t:.4g} s | queda {100 * queda:.1f} % del seno | "
          f"error máx = {error:.3e} °C ({100 * error / max(abs(A), 1e-30):.3f} % de A)")
    assert error <= TOLERANCIA * abs(A)

    # λ medido: cuánto se apagó el seno en el solver -> -ln(amplitud/A)/(alpha·t)
    if A == 0 or queda < 1e-3:
        print(f"  [{nombre}] λ teórico = {lam:.4f} 1/m² (el seno ya se apagó: no se mide λ)")
        return
    modo = (T_ex - T0) / (A * queda)
    lam_medido = -np.log(amplitud(T_solver, modo) / A) / (MATERIAL.alpha * t)
    print(f"  [{nombre}] λ teórico = {lam:.4f} 1/m², λ medido = {lam_medido:.4f} 1/m² "
          f"({100 * (lam_medido / lam - 1):+.2f} %)")
    assert lam_medido == pytest.approx(lam, rel=0.03)


def energia_en_grados(dE, malla):
    """Una diferencia de energía [J] expresada como cambio de T promedio [°C]."""
    return dE / malla.capacidad.sum()


# -----------------------------------------------------------------------------
# 1, 2, 3: soluciones exactas
# -----------------------------------------------------------------------------

def test_1_exacta_1d():
    """T(x,t) = T0 + A·sin(pi·x/Lx)·exp(-alpha·λ·t), λ = (pi/Lx)². Solver 1D."""
    malla = solver1d.Malla1D.uniforme(Lx, N, MATERIAL)
    T_ini = exacta.seno_dirichlet(malla.x, 0.0, Lx, MATERIAL.alpha, T0, A)
    res = solver1d.resolver(malla, T_ini, solver1d.Dirichlet(T0), solver1d.Dirichlet(T0),
                            t_final=tiempo_final([Lx]), dt=DT, n_guardados=1)
    print(f"\n  [1D] N = {N}, dx = {malla.dx:.4g} m, dt = {res.dt:.4g} s, Fo = {res.fourier:.3f}")
    comprobar_seno("1D", [malla.x], [Lx], res.T[-1], res.t[-1])


def test_2_seno_2d():
    """T = T0 + A·sin(pi·x/Lx)·sin(pi·y/Ly)·exp(-alpha·λ·t), λ = (pi/Lx)² + (pi/Ly)²."""
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, 1), MATERIAL)
    x, y, _ = malla.centros
    X, Y = np.meshgrid(x, y, indexing="ij")
    T_ini = exacta.seno_dirichlet_nd((X, Y), 0.0, (Lx, Ly), MATERIAL.alpha, T0, A)[:, :, None]
    bordes = {c: Dirichlet(T0) for c in CARAS[:4]}       # z± aisladas: 2D
    res = resolver(malla, T_ini, bordes, tiempo_final((Lx, Ly)), dt=DT, n_guardados=1)
    print(f"\n  [2D] N = {N}x{N}, dt = {res.dt:.4g} s, Fo = {res.fourier:.3f}")
    comprobar_seno("2D", [x, y], (Lx, Ly), res.T[-1][:, :, 0], res.t[-1])


def test_3_seno_3d():
    """T = T0 + A·sin·sin·sin·exp(-alpha·λ·t), λ = (pi/Lx)² + (pi/Ly)² + (pi/Lz)²."""
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, N), MATERIAL)
    centros = malla.centros
    T_ini = exacta.seno_dirichlet_nd(np.meshgrid(*centros, indexing="ij"), 0.0,
                                     (Lx, Ly, Lz), MATERIAL.alpha, T0, A)
    bordes = {c: Dirichlet(T0) for c in CARAS}
    res = resolver(malla, T_ini, bordes, tiempo_final((Lx, Ly, Lz)), dt=DT, n_guardados=1)
    print(f"\n  [3D] N = {N}x{N}x{N}, dt = {res.dt:.4g} s, Fo = {res.fourier:.3f}")
    comprobar_seno("3D", centros, (Lx, Ly, Lz), res.T[-1], res.t[-1])


# -----------------------------------------------------------------------------
# 4, 5, 6: energía
# -----------------------------------------------------------------------------

def test_4_energia_adiabatica():
    """Todas las caras aisladas, sin fuente: la energía E = sum(rho·cp·V·T) no cambia
    y la temperatura tiende al promedio inicial."""
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, N), MATERIAL)
    T_ini = T0 + A * np.random.default_rng(0).random(malla.forma)  # campo cualquiera
    res = resolver(malla, T_ini, {}, tiempo_final((Lx, Ly, Lz)), dt=DT, n_guardados=10)
    E = res.energia()
    deriva = energia_en_grados(np.max(np.abs(E - E[0])), malla)
    print(f"\n  [energía] E0 = {E[0]:.6e} J | máxima variación = {np.max(np.abs(E - E[0])):.3e} J "
          f"(= {deriva:.2e} °C de T promedio) | dispersión de T: {np.ptp(T_ini):.2f} -> "
          f"{np.ptp(res.T[-1]):.2f} °C")
    assert deriva < 1e-9
    assert res.T[-1].mean() == pytest.approx(T_ini.mean())


def test_5_fuente_uniforme():
    """Bloque aislado con fuente q0: dT/dt = q0/(rho·cp) en TODAS las celdas,
    y la energía crece exactamente q0·V·t."""
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, N), MATERIAL)
    res = resolver(malla, T0, {}, tiempo_final((Lx, Ly, Lz)), dt=DT, fuente=q0, n_guardados=10)
    T_ex = exacta.adiabatico_fuente_uniforme(res.t, T0, q0, rho, cp)
    error = np.max(np.abs(res.T - T_ex[:, None, None, None]))
    ganada = res.energia() - res.energia()[0]
    esperada = q0 * Lx * Ly * Lz * res.t
    print(f"\n  [fuente] dT/dt = q0/(rho·cp) = {q0 / (rho * cp):.4g} °C/s | T final = "
          f"{res.T[-1].mean():.4f} °C (exacta {T_ex[-1]:.4f}) | error máx = {error:.2e} °C | "
          f"energía ganada = {ganada[-1]:.4e} J (q0·V·t = {esperada[-1]:.4e} J)")
    assert error < 1e-9 * max(1.0, abs(T_ex[-1]))
    assert ganada == pytest.approx(esperada, rel=1e-9, abs=1e-6)


def test_6_neumann():
    """Entra q_SUP [W/m²] por la cara superior y el resto está aislado:
    la energía crece exactamente q·A·t y la parte de arriba se calienta primero."""
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, N), MATERIAL)
    res = resolver(malla, T0, {"z+": Neumann(q_SUP)}, tiempo_final((Lx, Ly, Lz)),
                   dt=DT, n_guardados=10)
    ganada = res.energia() - res.energia()[0]
    esperada = q_SUP * Lx * Ly * res.t
    arriba, abajo = res.T[-1][:, :, -1].mean(), res.T[-1][:, :, 0].mean()
    print(f"\n  [Neumann] energía ganada = {ganada[-1]:.4e} J (q·A·t = {esperada[-1]:.4e} J) | "
          f"T capa superior = {arriba:.3f} °C, capa inferior = {abajo:.3f} °C")
    assert ganada == pytest.approx(esperada, rel=1e-9, abs=1e-6)
    if q_SUP != 0:
        assert np.sign(arriba - abajo) == np.sign(q_SUP)  # entra calor -> arriba más caliente


# -----------------------------------------------------------------------------
# 7: Robin
# -----------------------------------------------------------------------------

def test_7_robin():
    """(a) Estacionario: x- a T0 + A, x+ con aire (h, T_INF, R_EXTRA), resto aislado.
        Debe dar el perfil lineal de resistencias en serie
        q = (T_izq - T_INF)/(Lx/k + R_EXTRA + 1/h).
    (b) Signo: bloque a T0 con aire a T_INF en todas las caras: el calor va del
        más caliente al más frío y ninguna celda sale del rango [T0, T_INF]."""
    # (a) estado estacionario (dt automático: se simulan 10·Lx²/alpha)
    T_izq = T0 + A
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, 2, 2), MATERIAL)
    bordes = {"x-": Dirichlet(T_izq), "x+": Robin(h, T_INF, R_EXTRA)}
    res = resolver(malla, T0, bordes, 10 * Lx**2 / MATERIAL.alpha, n_guardados=1)
    T_ex = exacta.estacionario_dirichlet_robin(malla.centros[0], Lx, k, T_izq, h, T_INF, R_EXTRA)
    error = np.max(np.abs(res.T[-1][:, 0, 0] - T_ex))
    q = (T_izq - T_INF) / (Lx / k + R_EXTRA + 1 / h)
    print(f"\n  [Robin a] q = (T_izq - T_INF)/(L/k + R + 1/h) = {q:.4g} W/m² | "
          f"T superficie = T_izq - q·L/k = {T_izq - q * Lx / k:.3f} °C | "
          f"error del perfil = {error:.2e} °C")
    assert error < 1e-6 * max(1.0, abs(T_izq - T_INF))

    # (b) signo del flujo
    malla = Malla3D.uniforme((Lx, Ly, Lz), (N, N, N), MATERIAL)
    res = resolver(malla, T0, {c: Robin(h, T_INF, R_EXTRA) for c in CARAS},
                   tiempo_final((Lx, Ly, Lz)), dt=DT, n_guardados=10)
    dE = np.diff(res.energia())
    print(f"  [Robin b] bloque a {T0} °C, aire a {T_INF} °C: la energía "
          f"{'baja' if dE.sum() < 0 else 'sube' if dE.sum() > 0 else 'no cambia'} "
          f"({dE.sum():.3e} J) | T entre {res.T.min():.2f} y {res.T.max():.2f} °C")
    redondeo = 1e-12 * malla.capacidad.sum() * max(abs(T_INF - T0), 1.0)
    assert np.all(np.sign(T_INF - T0) * dE >= -redondeo)    # nunca va al revés
    assert np.sign(round(dE.sum() / redondeo)) == np.sign(T_INF - T0)
    assert res.T.min() >= min(T0, T_INF) - 1e-9
    assert res.T.max() <= max(T0, T_INF) + 1e-9
