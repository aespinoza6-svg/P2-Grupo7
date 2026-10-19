"""
Tests de verificación del solver 1D (mcplan.solver1d).

Cada test compara el solver con algo que sabemos que es cierto:
una solución exacta, una ley de conservación o un límite de estabilidad.
Correr con:   pytest -v
"""

import numpy as np
import pytest

from mcplan import soluciones_analiticas as exacta
from mcplan.solver1d import (
    Dirichlet,
    Malla1D,
    Material,
    Neumann,
    Robin,
    paso_tiempo_maximo,
    resolver,
)

# Material "didáctico" con alpha = 1e-4 m²/s (como en los ejemplos del curso):
# los procesos ocurren en minutos/horas y los tests corren rápido.
DIDACTICO = Material("didactico", k=1.0, rho=1000.0, cp=10.0)
# Hormigón con los valores provisionales de config/fundacion.yaml.
HORMIGON = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)


# ---------------------------------------------------------------------------
# Propiedades y estabilidad
# ---------------------------------------------------------------------------

def test_difusividad_hormigon():
    """alpha = k/(rho·cp) = 2/(2400·900) ≈ 9.26e-7 m²/s."""
    assert HORMIGON.alpha == pytest.approx(9.259259e-7, rel=1e-6)


def test_paso_maximo_celda_interior_es_fourier_un_medio():
    """Con bordes adiabáticos todas las celdas cumplen Fo <= 1/2."""
    malla = Malla1D.uniforme(L=1.0, n=50, material=DIDACTICO)
    dt_max = paso_tiempo_maximo(malla, Neumann(0.0), Neumann(0.0))
    assert DIDACTICO.alpha * dt_max / malla.dx**2 == pytest.approx(0.5)


def test_paso_maximo_borde_dirichlet_es_fourier_un_tercio():
    """La celda junto a un Dirichlet tiene la superficie a dx/2
    (conductancia 2k/dx), así que su límite es más estricto: Fo <= 1/3."""
    malla = Malla1D.uniforme(L=1.0, n=50, material=DIDACTICO)
    dt_max = paso_tiempo_maximo(malla, Dirichlet(0.0), Neumann(0.0))
    assert DIDACTICO.alpha * dt_max / malla.dx**2 == pytest.approx(1.0 / 3.0)


def test_dt_inestable_lanza_error():
    """El solver se niega a correr con un dt mayor al estable."""
    malla = Malla1D.uniforme(L=1.0, n=20, material=DIDACTICO)
    dt_max = paso_tiempo_maximo(malla, Neumann(), Neumann())
    with pytest.raises(ValueError, match="supera el paso máximo estable"):
        resolver(malla, 20.0, Neumann(), Neumann(), t_final=100 * dt_max, dt=1.5 * dt_max)


def _escalon(malla):
    """Condición inicial discontinua: 20 °C a la izquierda, 60 °C a la derecha."""
    return np.where(malla.x < 0.5 * malla.L, 20.0, 60.0)


def test_esquema_estable_respeta_principio_del_maximo():
    """Con Fo <= 1/2 y sin fuentes, T nunca sale del rango inicial [20, 60]
    (el calor no puede crear temperaturas más extremas que las existentes)."""
    malla = Malla1D.uniforme(L=1.0, n=40, material=DIDACTICO)
    res = resolver(malla, _escalon(malla), Neumann(), Neumann(), t_final=2000.0)
    assert res.T.min() >= 20.0 - 1e-9
    assert res.T.max() <= 60.0 + 1e-9


def test_esquema_inestable_explota():
    """Con Fo = 0.6 > 1/2 el modo de mayor frecuencia se amplifica por
    |1 - 4·Fo| = 1.4 en cada paso: la solución explota (no es física)."""
    malla = Malla1D.uniforme(L=1.0, n=40, material=DIDACTICO)
    dt = 0.6 * malla.dx**2 / DIDACTICO.alpha
    res = resolver(malla, _escalon(malla), Neumann(), Neumann(),
                   t_final=200 * dt, dt=dt, permitir_inestable=True)
    assert np.abs(res.T[-1]).max() > 1e6


# ---------------------------------------------------------------------------
# Solución exacta y convergencia
# ---------------------------------------------------------------------------

def _error_seno(n, fourier=0.25, t_final=1500.0):
    """Error máximo |T_num - T_exacta| del caso seno con n celdas y Fo fijo.
    (Fo = 0.25 < 1/3: las celdas junto a los bordes Dirichlet lo exigen.)"""
    L, T_b, A = 1.0, 20.0, 50.0
    malla = Malla1D.uniforme(L, n, DIDACTICO)
    T0 = exacta.seno_dirichlet(malla.x, 0.0, L, DIDACTICO.alpha, T_b, A)
    dt = fourier * malla.dx**2 / DIDACTICO.alpha
    res = resolver(malla, T0, Dirichlet(T_b), Dirichlet(T_b), t_final=t_final, dt=dt)
    T_ex = exacta.seno_dirichlet(res.x, res.t[-1], L, DIDACTICO.alpha, T_b, A)
    return np.max(np.abs(res.T[-1] - T_ex))


def test_solucion_exacta_seno_dirichlet():
    """Con 50 celdas el error es menor a 0.05 °C sobre una amplitud de 50 °C."""
    assert _error_seno(50) < 0.05


def test_convergencia_espacial_orden_2():
    """Al duplicar el número de celdas (con Fo fijo, dt ∝ dx²) el error debe
    caer ~4 veces: el esquema es de orden 2 en espacio."""
    errores = [_error_seno(n) for n in (10, 20, 40, 80)]
    ordenes = np.log2(np.array(errores[:-1]) / np.array(errores[1:]))
    assert ordenes[-1] == pytest.approx(2.0, abs=0.2)


# ---------------------------------------------------------------------------
# Conservación de energía y fuentes
# ---------------------------------------------------------------------------

def test_conservacion_energia_bordes_adiabaticos():
    """Sin flujos por los bordes ni fuente, la energía total no cambia
    (a precisión de máquina) y T tiende al promedio inicial."""
    rng = np.random.default_rng(0)
    malla = Malla1D.uniforme(L=1.0, n=30, material=DIDACTICO)
    T0 = 20.0 + 30.0 * rng.random(malla.n)
    res = resolver(malla, T0, Neumann(), Neumann(), t_final=50000.0)
    E = res.energia()
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-12
    assert res.T[-1] == pytest.approx(np.full(malla.n, T0.mean()), abs=1e-3)


def test_balance_con_flujo_neumann_y_fuente():
    """Energía ganada = calor que entra por los bordes + calor generado:
    E(t) - E(0) = (q_izq + q_der)·t + Q·L·t   [J/m²]."""
    L, q_izq, Q = 2.0, 500.0, 100.0
    malla = Malla1D.uniforme(L, 40, HORMIGON)
    res = resolver(malla, 20.0, Neumann(q_izq), Neumann(0.0),
                   t_final=5 * 86400.0, fuente=Q)
    ganada = res.energia() - res.energia()[0]
    esperada = (q_izq + Q * L) * res.t
    assert ganada == pytest.approx(esperada, rel=1e-10)


def test_fuente_uniforme_adiabatica_aumento_lineal():
    """Bloque aislado con Q uniforme: T sube linealmente, igual en todo x."""
    Q = 1000.0  # W/m³, orden de magnitud del calor de hidratación temprano
    malla = Malla1D.uniforme(L=3.3, n=33, material=HORMIGON)
    res = resolver(malla, 20.0, Neumann(), Neumann(), t_final=3 * 86400.0, fuente=Q)
    T_ex = exacta.adiabatico_fuente_uniforme(res.t[-1], 20.0, Q, HORMIGON.rho, HORMIGON.cp)
    assert res.T[-1] == pytest.approx(np.full(malla.n, T_ex), rel=1e-12)


# ---------------------------------------------------------------------------
# Condición Robin y materiales distintos (estado estacionario)
# ---------------------------------------------------------------------------

def test_robin_tiende_a_temperatura_del_aire():
    """Expuesto a aire a 15 °C por ambos lados, todo termina a 15 °C."""
    malla = Malla1D.uniforme(L=1.0, n=20, material=DIDACTICO)
    res = resolver(malla, 60.0, Robin(10.0, 15.0), Robin(10.0, 15.0), t_final=2e5)
    assert res.T[-1] == pytest.approx(np.full(malla.n, 15.0), abs=1e-6)


def test_estacionario_dirichlet_robin_con_moldaje():
    """Régimen permanente con resistencias en serie
    (conducción + moldaje R_extra + película de aire 1/h): perfil lineal exacto."""
    L, k, T_izq, h, T_inf, R = 1.0, DIDACTICO.k, 60.0, 10.0, 20.0, 0.14
    malla = Malla1D.uniforme(L, 20, DIDACTICO)
    res = resolver(malla, 20.0, Dirichlet(T_izq), Robin(h, T_inf, R_extra=R), t_final=2e5)
    T_ex = exacta.estacionario_dirichlet_robin(res.x, L, k, T_izq, h, T_inf, R)
    assert res.T[-1] == pytest.approx(T_ex, abs=1e-6)


def test_estacionario_dos_materiales_en_serie():
    """Dos capas de distinta conductividad: la interfaz se trata con
    resistencias en serie y el perfil estacionario es exacto."""
    B = Material("B", k=0.5, rho=1000.0, cp=10.0)
    malla = Malla1D.por_capas([(0.6, DIDACTICO), (0.4, B)], n=50)
    res = resolver(malla, 0.0, Dirichlet(80.0), Dirichlet(10.0), t_final=2e5)
    T_ex = exacta.estacionario_dos_capas(res.x, 0.6, DIDACTICO.k, 0.4, B.k, 80.0, 10.0)
    assert res.T[-1] == pytest.approx(T_ex, abs=1e-6)


def test_capas_deben_caer_en_caras_de_celda():
    """Si una interfaz cae dentro de una celda, se pide otra malla."""
    with pytest.raises(ValueError, match="número entero"):
        Malla1D.por_capas([(0.55, DIDACTICO), (0.45, HORMIGON)], n=10)
