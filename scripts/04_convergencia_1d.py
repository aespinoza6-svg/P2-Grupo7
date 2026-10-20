"""
Paso 4 — Estudio de convergencia del solver 1D (verificación).

    python scripts/04_convergencia_1d.py

Caso: decaimiento del seno con bordes Dirichlet (tiene solución exacta).

1) Convergencia ESPACIAL: N = 10, 20, ..., 320 celdas con Fo fijo (dt ∝ dx²).
   Error vs solución exacta. Se espera orden 2: duplicar N divide el error por 4.

2) Convergencia TEMPORAL: malla fija (N = 50) y dt = dt_max, dt_max/2, ...
   Para aislar el error del tiempo se compara contra una referencia en la
   MISMA malla (así el error espacial es idéntico y se cancela). La
   referencia se obtiene por extrapolación de Richardson con dos pasos muy
   finos h y 2h: si T(dt) ≈ T_exacta_malla + C·dt, entonces
   2·T(h) - T(2h) elimina el término C·dt.
   Se espera orden 1 (Euler explícito): dividir dt por 2 divide el error por 2.

Salida: resultados/figuras/fig06_convergencia.png y resultados/convergencia_1d.csv
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mcplan import analiticas as exacta
from mcplan.solver1d import Dirichlet, Malla1D, Material, paso_tiempo_maximo, resolver

# --- Estilo y salida de las figuras -----------------------------------------
RAIZ = Path(__file__).resolve().parents[1]
CARPETA_FIGURAS = RAIZ / "resultados" / "figuras"
# Colores para curvas que representan COSAS distintas; GRIS para límites y referencias.
AZUL, NARANJO, VERDE_AGUA, GRIS = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "savefig.dpi": 160})


def guardar(fig, nombre):
    CARPETA_FIGURAS.mkdir(parents=True, exist_ok=True)
    ruta = CARPETA_FIGURAS / nombre
    fig.savefig(ruta, bbox_inches="tight")
    print(f"  figura -> {ruta.relative_to(RAIZ)}")


DIDACTICO = Material("didactico", k=1.0, rho=1000.0, cp=10.0)
L, T_B, A, T_FINAL = 1.0, 20.0, 50.0, 1500.0


def simular(n, dt):
    malla = Malla1D.uniforme(L, n, DIDACTICO)
    T0 = exacta.seno_dirichlet(malla.x, 0.0, L, DIDACTICO.alpha, T_B, A)
    return resolver(malla, T0, Dirichlet(T_B), Dirichlet(T_B), t_final=T_FINAL, dt=dt, n_guardados=1)


def orden_observado(h, error):
    """p = log(e1/e2)/log(h1/h2) entre resoluciones consecutivas."""
    return np.log(error[:-1] / error[1:]) / np.log(h[:-1] / h[1:])


def main():
    # --- 1) Espacio --------------------------------------------------------
    fo = 0.25
    ns = np.array([10, 20, 40, 80, 160, 320])
    dxs, err_x = [], []
    for n in ns:
        dx = L / n
        res = simular(n, fo * dx**2 / DIDACTICO.alpha)
        T_ex = exacta.seno_dirichlet(res.x, res.t[-1], L, DIDACTICO.alpha, T_B, A)
        dxs.append(dx)
        err_x.append(np.max(np.abs(res.T[-1] - T_ex)))
    dxs, err_x = np.array(dxs), np.array(err_x)
    p_x = orden_observado(dxs, err_x)

    # --- 2) Tiempo ---------------------------------------------------------
    n = 50
    malla = Malla1D.uniforme(L, n, DIDACTICO)
    dt_max = paso_tiempo_maximo(malla, Dirichlet(T_B), Dirichlet(T_B))
    h = dt_max / 256
    referencia = 2 * simular(n, h).T[-1] - simular(n, 2 * h).T[-1]  # Richardson
    dts = dt_max / 2.0 ** np.arange(0, 6)
    err_t = np.array([np.max(np.abs(simular(n, dt).T[-1] - referencia)) for dt in dts])
    p_t = orden_observado(dts, err_t)

    # --- Tabla -------------------------------------------------------------
    print("Convergencia espacial (Fo = 0.25)")
    print("     N        dx [m]     error [°C]   orden")
    for i, (nn, dx, e) in enumerate(zip(ns, dxs, err_x)):
        print(f"  {nn:4d}   {dx:.5f}   {e:.3e}   {'' if i == 0 else f'{p_x[i - 1]:.3f}'}")
    print("\nConvergencia temporal (N = 50, referencia de Richardson)")
    print("     dt [s]     error [°C]   orden")
    for i, (dt, e) in enumerate(zip(dts, err_t)):
        print(f"  {dt:8.4f}   {e:.3e}   {'' if i == 0 else f'{p_t[i - 1]:.3f}'}")

    csv = RAIZ / "resultados" / "convergencia_1d.csv"
    with open(csv, "w", encoding="utf-8") as f:
        f.write("estudio,paso,error_C\n")
        f.writelines(f"espacial,{dx},{e}\n" for dx, e in zip(dxs, err_x))
        f.writelines(f"temporal,{dt},{e}\n" for dt, e in zip(dts, err_t))
    print(f"\n  tabla  -> {csv.relative_to(RAIZ)}")

    # --- Figura ------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    a1.loglog(dxs, err_x, "o-", color=AZUL, ms=6, label="error del solver")
    a1.loglog(dxs, err_x[0] * (dxs / dxs[0]) ** 2, "--", color=GRIS, lw=1.2, label="pendiente 2 (referencia)")
    a1.set(xlabel="dx [m]", ylabel="error máximo [°C]", title=f"Espacio: orden observado ≈ {p_x[-1]:.2f}")
    a1.legend(fontsize=9)
    a2.loglog(dts, err_t, "o-", color=NARANJO, ms=6, label="error del solver")
    a2.loglog(dts, err_t[0] * (dts / dts[0]), "--", color=GRIS, lw=1.2, label="pendiente 1 (referencia)")
    a2.set(xlabel="dt [s]", ylabel="error máximo [°C]", title=f"Tiempo: orden observado ≈ {p_t[-1]:.2f}")
    a2.legend(fontsize=9)
    fig.suptitle("Convergencia del solver 1D (caso seno con bordes Dirichlet)")
    guardar(fig, "fig06_convergencia.png")


if __name__ == "__main__":
    main()
