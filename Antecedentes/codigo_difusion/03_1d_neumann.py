"""
1D: Dirichlet a la izquierda + Neumann adiabático a la derecha.

Muestra perfiles temporales para experimentar de forma interactiva.
Figura: images/02_1d_neumann.png
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "images"   # figuras que se insertan en esta guía
OUT = ROOT / "output"   # campos y datos guardados
FIG.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

L = 1.0
alpha = 1.0e-4
Tini = 20.0       # °C: temperatura inicial de la barra
T_left = 60.0     # °C: temperatura impuesta en x=0

nx = 101
x = np.linspace(0.0, L, nx)
dx = x[1] - x[0]

r_target = 0.45
dt = r_target*dx**2/alpha
t_end = 3*6000.0
nsteps = int(np.ceil(t_end/dt))
dt = t_end/nsteps
r = alpha*dt/dx**2

T = np.full(nx, Tini)

# x=0: Dirichlet T=T_left
# x=L: Neumann adiabático dT/dx=0
snapshot_steps = np.rint(np.linspace(0, nsteps, 20)).astype(int)
snapshots = {0: T.copy()}

for n in range(1, nsteps+1):
    Told = T.copy()

    T[1:-1] = Told[1:-1] + r*(
        Told[2:] - 2*Told[1:-1] + Told[:-2]
    )

    T[0] = T_left

    # Ghost node:
    # dT/dx = 0 => T_ghost = T[-2]
    # Laplaciano en el borde => 2(T[-2]-T[-1])/dx^2
    T[-1] = Told[-1] + 2*r*(Told[-2] - Told[-1])

    if n in snapshot_steps:
        snapshots[n] = T.copy()

fig, ax = plt.subplots(figsize=(8.5, 5.0))
colors = plt.cm.viridis(np.linspace(0.05, 0.90, len(snapshot_steps)))

for step, color in zip(snapshot_steps, colors):
    t = step * dt
    ax.plot(x, snapshots[step], "o-", color=color, lw=2.0, ms=3.5,
            mfc="white", mew=0.9, markevery=5,
            label=f"t = {t/60:.0f} min")

ax.axhline(T_left, color="black", ls="--", lw=1.0,
           label=fr"$T_{{izq}}={T_left:.0f}$ °C")
ax.set_xlabel("x [m]")
ax.set_ylabel("T [°C]")
ax.set_title("1D: borde izquierdo fijo y borde derecho adiabático")
ax.legend(title="Perfil", loc="best", ncol=2, fontsize=9, title_fontsize=10)
ax.grid(True, alpha=0.35)
fig.tight_layout()
# plt.savefig(FIG/"02_1d_neumann.png", dpi=180)
plt.show()
plt.close(fig)

print(f"T en extremo adiabático al final: {T[-1]:.3f} °C")
print(f"El estado estacionario esperado es T={T_left:.0f} °C en toda la barra.")
