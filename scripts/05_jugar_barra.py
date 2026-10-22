"""
Paso 5 — Laboratorio: jugar con la ecuación del calor en una barra.

    python scripts/05_jugar_barra.py               # calcula y abre el gráfico
    python scripts/05_jugar_barra.py --sin-ventana # solo guarda el gráfico

Ecuación que se resuelve (en mcplan/solver1d.py):

    rho·cp · dT/dt = d/dx( k · dT/dx ) + Q

Cambia los valores de la sección PERILLAS, guarda el archivo y vuelve a correr.
El programa muestra la temperatura del CENTRO de la barra en el tiempo y los
perfiles T(x). Gráfico en resultados/tmp/jugar_barra.png (no se sube a GitHub).
"""

import sys
from pathlib import Path

import matplotlib

if "--sin-ventana" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mcplan.solver1d import Dirichlet, Malla1D, Material, Neumann, Robin, resolver

# =============================================================================
# PERILLAS: cambia estos valores
# =============================================================================

L = 1.0                    # largo de la barra [m]
N = 50                     # número de celdas -> dx = L/N (más celdas = más preciso y más lento)
DT = None                  # paso de tiempo [s]. None = automático (90 % del máximo estable).
                           # Si pones uno mayor al máximo, el programa avisa con un error.

# Material. Didáctico: alpha = k/(rho·cp) = 1e-4 m²/s (procesos de minutos).
# Hormigón: k=2.0, rho=2400, cp=900 (alpha ≈ 9.3e-7 m²/s, procesos de días).
k = 1.0                    # conductividad [W/(m K)]: conduce
rho = 1000.0               # densidad [kg/m³]
cp = 10.0                  # calor específico [J/(kg K)]: almacena

T_INICIAL = 20.0           # condición inicial [°C] (toda la barra igual)

# Condiciones de borde. Opciones:
#   Dirichlet(T)                  temperatura fija en el extremo [°C]
#   Neumann(q)                    flujo que ENTRA [W/m²]; Neumann(0) = aislado
#   Robin(h, T_aire)              aire a T_aire con convección h [W/(m² K)]
#   Robin(h, T_aire, R_extra=R)   + moldaje o aislación (R = espesor/k) [m² K/W]
BORDE_IZQ = Dirichlet(60.0)
BORDE_DER = Dirichlet(20.0)

Q = 0.0                    # fuente de calor interna [W/m³] (0 = sin fuente)

T_FINAL_MIN = 120.0        # cuánto tiempo simular [minutos]

OBJETIVO_CENTRO = None     # p. ej. 35.0: avisa cuándo el centro llega a esa T [°C]

# =============================================================================


def main():
    material = Material("barra", k=k, rho=rho, cp=cp)
    malla = Malla1D.uniforme(L, N, material)
    res = resolver(malla, T_INICIAL, BORDE_IZQ, BORDE_DER, t_final=T_FINAL_MIN * 60.0,
                   fuente=Q, dt=DT, n_guardados=200)

    # Temperatura en x = L/2 (interpolada entre los centros de celda).
    T_centro = np.array([np.interp(L / 2, res.x, perfil) for perfil in res.T])
    minutos = res.t / 60.0

    print(f"alpha = {material.alpha:.3e} m²/s   dx = {malla.dx:.4f} m   "
          f"dt = {res.dt:.4g} s (máximo estable {res.dt_max:.4g} s)   "
          f"Fo = {res.fourier:.3f}   pasos = {res.n_pasos}")
    print(f"T centro: inicial {T_centro[0]:.2f} °C -> final {T_centro[-1]:.2f} °C "
          f"(máx {T_centro.max():.2f} °C)")
    print(f"Tiempo característico L²/alpha = {L**2 / material.alpha / 60:.1f} min")

    # Estado estacionario teórico, solo si ambos bordes son Dirichlet:
    # T(x) = recta entre T_izq y T_der  +  parábola de la fuente  Q·x·(L-x)/(2k)
    if (isinstance(BORDE_IZQ, Dirichlet) and isinstance(BORDE_DER, Dirichlet)
            and not any(callable(v) for v in (BORDE_IZQ.T, BORDE_DER.T, Q))):
        T_teo = (BORDE_IZQ.T + BORDE_DER.T) / 2 + Q * L**2 / (8 * k)
        print(f"T centro estacionaria teórica = (T_izq + T_der)/2 + Q·L²/(8k) = {T_teo:.2f} °C")

    if OBJETIVO_CENTRO is not None:
        alcanzado = np.nonzero(T_centro >= OBJETIVO_CENTRO)[0]
        if len(alcanzado):
            print(f"El centro llega a {OBJETIVO_CENTRO} °C a los {minutos[alcanzado[0]]:.1f} min")
        else:
            print(f"El centro NO llega a {OBJETIVO_CENTRO} °C en {T_FINAL_MIN} min")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    a1.plot(minutos, T_centro, color="#2a78d6", lw=2)
    if OBJETIVO_CENTRO is not None:
        a1.axhline(OBJETIVO_CENTRO, color="#52514e", ls="--", lw=1)
    a1.set(xlabel="tiempo [min]", ylabel="T [°C]", title="Temperatura del centro de la barra")
    a1.grid(alpha=0.3)
    indices = np.linspace(0, len(res.t) - 1, 6).astype(int)
    for i, c in zip(indices, plt.cm.Blues(np.linspace(0.35, 0.95, len(indices)))):
        a2.plot(res.x, res.T[i], color=c, lw=2, label=f"t = {minutos[i]:.0f} min")
    a2.axvline(L / 2, color="#52514e", ls=":", lw=1)
    a2.set(xlabel="x [m]", ylabel="T [°C]", title="Perfiles de temperatura")
    a2.legend(fontsize=8)
    a2.grid(alpha=0.3)
    fig.tight_layout()

    salida = Path(__file__).resolve().parents[1] / "resultados" / "tmp" / "jugar_barra.png"
    salida.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(salida, dpi=130)
    print(f"Gráfico: {salida}")
    if "--sin-ventana" not in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
