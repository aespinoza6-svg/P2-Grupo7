"""
Verificación 1D y evolución temporal (FTCS, Dirichlet homogéneo).

La figura superpone perfiles a varios instantes. Línea continua: solución
exacta. Marcadores: solución numérica.

Figura: images/01_1d_dirichlet.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "images"   # figuras que se insertan en esta guía
OUT = ROOT / "output"   # campos y datos guardados
FIG.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

L = 1.0                 # m
alpha = 1.0e-4          # m^2/s
Tizq = 0.0               # °C
Tder = 10.0               # °C
Tini = 20.
A = 0*50.0                # °C

nx = 101
x = np.linspace(0.0, L, nx)
dx = x[1] - x[0]

# FTCS explícito. En 1D: r <= 1/2
r_target = 0.45
dt = r_target * dx**2 / alpha
t_end = 3*1500.0
nsteps = int(np.ceil(t_end / dt))
dt = t_end / nsteps
r = alpha * dt / dx**2

T = Tini + A * np.sin(np.pi * x / L)
T[0] = Tizq
T[-1] = Tder

# Seis perfiles equiespaciados en el tiempo para mostrar la difusión.
snapshot_steps = np.rint(np.linspace(0, nsteps, 20)).astype(int)
snapshots = {0: T.copy()}

for n in range(1, nsteps + 1):
    Told = T.copy()
    T[1:-1] = Told[1:-1] + r * (
        Told[2:] - 2.0*Told[1:-1] + Told[:-2]
    )
    T[0] = Tizq
    T[-1] = Tder

    if n in snapshot_steps:
        snapshots[n] = T.copy()

# Texact = T0 + A*np.sin(np.pi*x/L) * np.exp(
#     -alpha*(np.pi/L)**2*t_end
# )

# err_inf = np.max(np.abs(T - Texact))
print(f"dx = {dx:.5e} m")
print(f"dt = {dt:.5e} s")
print(f"r  = {r:.5f}")
# print(f"Error L_inf = {err_inf:.6e} °C")

fig, ax = plt.subplots(figsize=(8.5, 5.0))
colors = plt.cm.viridis(np.linspace(0.05, 0.90, len(snapshot_steps)))

for step, color in zip(snapshot_steps, colors):
    t = step * dt
    # exact = T0 + A*np.sin(np.pi*x/L) * np.exp(
    #     -alpha*(np.pi/L)**2*t
    # )
    # ax.plot(x, exact, color=color, lw=2.0,
            # label=f"t = {t/60:.0f} min")
    ax.plot(x, snapshots[step], color=color, lw=2.0,
            label=f"t = {t/60:.0f} min")

time_legend = ax.legend(title="Perfil", loc="upper right", ncol=2,
                        fontsize=9, title_fontsize=10)
ax.add_artist(time_legend)
ax.legend(handles=[
    Line2D([0], [0], color="black", lw=2, label="exacta"),
    Line2D([0], [0], marker="o", color="black", lw=0, mfc="white",
           label="numérica (cada 5 nodos)"),
], loc="upper center", fontsize=9)

plt.xlabel("x [m]")
plt.ylabel("T [°C]")
plt.title("Difusión 1D: decaimiento del perfil inicial")
plt.grid(True, alpha=0.35)
plt.tight_layout()
# plt.savefig(FIG / "01_1d_dirichlet.png", dpi=180)
plt.show()
plt.close(fig)
