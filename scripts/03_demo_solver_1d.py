"""
Paso 3 — Demostraciones del solver 1D (figuras para la revisión oral).

    python scripts/03_demo_solver_1d.py            # guarda las figuras
    python scripts/03_demo_solver_1d.py --mostrar  # además abre las ventanas

Genera en resultados/figuras/:
  fig01  Dirichlet: solución numérica vs exacta (verificación)
  fig02  Dirichlet + Neumann adiabático (calentamiento desde un borde)
  fig03  Robin con aire que oscila cada 24 h (la onda diaria penetra poco)
  fig04  Fuente de calor ILUSTRATIVA: núcleo vs superficie, con y sin aislación
  fig05  Estabilidad: Fo = 0.45 (estable) vs Fo = 0.55 (inestable)

Las figuras 03 y 04 usan una losa de hormigón de 3,30 m. Son ilustraciones
de la física: el clima es una onda ficticia y la fuente NO es todavía el
modelo de hidratación (semana 3).
"""

import argparse
import sys

import matplotlib

if "--mostrar" not in sys.argv:
    matplotlib.use("Agg")  # sin ventanas: solo guardar archivos
import matplotlib.pyplot as plt
import numpy as np

from estilo_graficos import AZUL, GRIS, NARANJO, VERDE_AGUA, colores_tiempo, guardar
from mcplan import soluciones_analiticas as exacta
from mcplan.solver1d import Dirichlet, Malla1D, Material, Neumann, Robin, resolver

HORA = 3600.0
DIA = 24 * HORA

# Material didáctico (alpha = 1e-4 m²/s): procesos de minutos, fácil de ver.
DIDACTICO = Material("didactico", k=1.0, rho=1000.0, cp=10.0)
# Hormigón: valores provisionales de config/fundacion.yaml (alpha ≈ 9.3e-7 m²/s).
HORMIGON = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)


def fig01_dirichlet_vs_exacta():
    L, T_b, A = 1.0, 20.0, 50.0
    malla = Malla1D.uniforme(L, 40, DIDACTICO)
    T0 = exacta.seno_dirichlet(malla.x, 0.0, L, DIDACTICO.alpha, T_b, A)
    minutos = [0, 5, 10, 20, 40, 80]
    res = resolver(malla, T0, Dirichlet(T_b), Dirichlet(T_b), t_final=80 * 60.0,
                   tiempos_guardado=[m * 60.0 for m in minutos])

    fig, ax = plt.subplots(figsize=(8, 4.8))
    x_fino = np.linspace(0, L, 300)
    for T, t, c in zip(res.T, res.t, colores_tiempo(len(res.t))):
        ax.plot(x_fino, exacta.seno_dirichlet(x_fino, t, L, DIDACTICO.alpha, T_b, A), color=c,
                label=f"t = {t / 60:.0f} min")
        ax.plot(res.x[::2], T[::2], "o", color=c, mfc="white", ms=5, mew=1.3)
    error = np.max(np.abs(res.T[-1] - exacta.seno_dirichlet(res.x, res.t[-1], L, DIDACTICO.alpha, T_b, A)))
    ax.set(xlabel="x [m]", ylabel="T [°C]",
           title=f"Dirichlet en ambos bordes: línea = exacta, círculos = solver "
                 f"(N = 40, Fo = {res.fourier:.2f}, error máx. final = {error:.3f} °C)")
    ax.legend(title="instante", fontsize=9)
    guardar(fig, "fig01_dirichlet_vs_exacta.png")


def fig02_dirichlet_neumann():
    malla = Malla1D.uniforme(1.0, 50, DIDACTICO)
    minutos = [0, 10, 30, 60, 120, 240]
    res = resolver(malla, 20.0, Dirichlet(60.0), Neumann(0.0), t_final=240 * 60.0,
                   tiempos_guardado=[m * 60.0 for m in minutos])

    fig, ax = plt.subplots(figsize=(8, 4.8))
    for T, t, c in zip(res.T, res.t, colores_tiempo(len(res.t))):
        ax.plot(res.x, T, color=c, label=f"t = {t / 60:.0f} min")
    ax.axhline(60.0, color=GRIS, ls="--", lw=1)
    ax.set(xlabel="x [m]", ylabel="T [°C]",
           title="x = 0: Dirichlet T = 60 °C  |  x = L: Neumann adiabático (dT/dx = 0)")
    ax.annotate("pendiente nula en el borde aislado", xy=(1.0, res.T[-1][-1]), xytext=(0.55, 30),
                arrowprops=dict(arrowstyle="->", color=GRIS), color=GRIS)
    ax.legend(title="instante", fontsize=9)
    guardar(fig, "fig02_dirichlet_neumann.png")


def fig03_robin_ambiente_diario():
    """Losa de 3,30 m. x = profundidad desde la superficie superior.
    Arriba: aire 15 ± 7 °C (onda ficticia de 24 h), h = 10 W/m²K.
    Abajo: adiabático (simplificación provisional del suelo)."""
    L, n = 3.30, 66          # dx = 5 cm -> el centro de la 1a celda está a 25 mm
    malla = Malla1D.uniforme(L, n, HORMIGON)

    def T_aire(t):
        return 15.0 + 7.0 * np.sin(2 * np.pi * t / DIA)

    res = resolver(malla, 15.0, Robin(10.0, T_aire), Neumann(0.0), t_final=6 * DIA, n_guardados=6 * 48)
    dias = res.t / DIA
    profundidades = [0.025, 0.175, 0.475]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    a1.plot(dias, T_aire(res.t), color=GRIS, ls="--", lw=1.2, label="aire")
    for z, c in zip(profundidades, [NARANJO, VERDE_AGUA, AZUL]):
        i = int(np.argmin(np.abs(res.x - z)))
        a1.plot(dias, res.T[:, i], color=c, label=f"{100 * res.x[i]:.1f} cm bajo la superficie")
    a1.set(xlabel="tiempo [días]", ylabel="T [°C]", title="Historia temporal")
    a1.legend(fontsize=9)

    # Envolvente: rango de T en el último día, en función de la profundidad.
    ultimo_dia = res.t >= 5 * DIA
    rango = res.T[ultimo_dia].max(axis=0) - res.T[ultimo_dia].min(axis=0)
    delta = np.sqrt(2 * HORMIGON.alpha / (2 * np.pi / DIA))  # profundidad de penetración
    a2.plot(100 * res.x[:20], rango[:20], "o-", color=AZUL, ms=4, label="solver")
    a2.plot(100 * res.x[:20], rango[0] * np.exp(-(res.x[:20] - res.x[0]) / delta), color=GRIS, ls="--",
            lw=1.2, label=f"decaimiento teórico exp(-z/δ), δ = {100 * delta:.0f} cm")
    a2.set(xlabel="profundidad [cm]", ylabel="oscilación diaria (máx - mín) [°C]",
           title="La onda diaria se amortigua en pocos centímetros")
    a2.legend(fontsize=9)
    fig.suptitle("Robin con aire oscilante (15 ± 7 °C, h = 10 W/m²K) — onda ficticia, no es clima real")
    guardar(fig, "fig03_robin_ambiente_diario.png")


def fig04_fuente_nucleo_superficie():
    """Losa de 3,30 m con una fuente ILUSTRATIVA Q(t) = Q0·exp(-t/tau)
    (aumento adiabático total Q0·tau/(rho·cp) = 40 °C). Compara superficie
    desnuda vs 5 cm de poliestireno (R = 0.05/0.035 ≈ 1.43 m²K/W)."""
    L, n = 3.30, 66
    malla = Malla1D.uniforme(L, n, HORMIGON)
    Q0, tau = 1000.0, 1.0 * DIA

    def fuente(t, x):
        return Q0 * np.exp(-t / tau)

    casos = [("sin aislación", 0.0, NARANJO), ("con 5 cm de poliestireno", 0.05 / 0.035, AZUL)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for nombre, R, c in casos:
        res = resolver(malla, 20.0, Robin(10.0, 15.0, R_extra=R), Neumann(0.0),
                       t_final=14 * DIA, fuente=fuente, n_guardados=14 * 24)
        dias = res.t / DIA
        T_max = res.T.max(axis=1)        # máximo en el hormigón en cada instante
        T_25 = res.T[:, 0]               # celda cuyo centro está a 25 mm de la superficie
        a1.plot(dias, T_max, color=c, label=f"T máx ({nombre})")
        a1.plot(dias, T_25, color=c, ls=":", label=f"T a 25 mm ({nombre})")
        a2.plot(dias, T_max - T_25, color=c, label=nombre)
        print(f"  {nombre:<26} T_max = {T_max.max():5.1f} °C   max(T_max - T_25mm) = {(T_max - T_25).max():5.1f} °C")
    a1.axhline(70, color=GRIS, ls="--", lw=1)
    a1.text(0.2, 70.8, "límite 70 °C", color=GRIS)
    a2.axhline(20, color=GRIS, ls="--", lw=1)
    a2.text(0.2, 20.6, "límite 20 °C", color=GRIS)
    a1.set(xlabel="tiempo [días]", ylabel="T [°C]", title="Núcleo vs superficie", ylim=(10, 75))
    a2.set(xlabel="tiempo [días]", ylabel="T máx - T a 25 mm [°C]", title="Diferencial núcleo-superficie")
    a1.legend(fontsize=8, loc="center right")
    a2.legend(fontsize=9)
    fig.suptitle("Fuente ILUSTRATIVA (no es el modelo de hidratación) — aire 15 °C, base adiabática")
    guardar(fig, "fig04_fuente_nucleo_superficie.png")


def fig05_estabilidad():
    malla = Malla1D.uniforme(1.0, 40, DIDACTICO)
    T0 = np.where(malla.x < 0.5, 20.0, 60.0)
    fig, ejes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=False)
    for ax, fo in zip(ejes, [0.45, 0.55]):
        dt = fo * malla.dx**2 / DIDACTICO.alpha
        pasos = [0, 5, 10, 20]
        res = resolver(malla, T0, Neumann(), Neumann(), t_final=20 * dt, dt=dt,
                       tiempos_guardado=[p * dt for p in pasos], permitir_inestable=True)
        for T, p, c in zip(res.T, pasos, colores_tiempo(len(pasos))):
            ax.plot(res.x, T, "o-", color=c, ms=3, lw=1.5, label=f"{p} pasos")
        estado = "estable" if fo <= 0.5 else "INESTABLE"
        ax.set(xlabel="x [m]", ylabel="T [°C]", title=f"Fo = {fo} ({estado})")
        ax.legend(fontsize=9)
    fig.suptitle("Esquema explícito: si Fo > 1/2 aparece un serrucho que crece en cada paso")
    guardar(fig, "fig05_estabilidad.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mostrar", action="store_true", help="abrir las figuras en ventanas")
    args = parser.parse_args()

    for f in (fig01_dirichlet_vs_exacta, fig02_dirichlet_neumann, fig03_robin_ambiente_diario,
              fig04_fuente_nucleo_superficie, fig05_estabilidad):
        print(f.__name__)
        f()
    if args.mostrar:
        plt.show()


if __name__ == "__main__":
    main()
