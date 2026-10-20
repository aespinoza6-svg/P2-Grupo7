"""
1D: losa con condición Robin y ambiente oscilante en el lado izquierdo.

El sólido parte a 20 °C. El aire a la izquierda oscila entre 10 y 30 °C;
el aire del lado derecho se mantiene a 20 °C.

Abre un gráfico interactivo con historia temporal y perfiles T(x) para
experimentar con h, amplitud, período y tiempos de muestreo.
Figura: images/03_1d_robin.png
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "images"   # figuras que se insertan en esta guía
OUT = ROOT / "output"   # campos y datos guardados
FIG.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

# Losa 1D inicialmente caliente, expuesta por ambos extremos a convección.
L = 1.0             # m
k = 2.0             # W/(m K)
rho = 2400.0        # kg/m^3
cp = 900.0          # J/(kg K)
alpha = k/(rho*cp)

h = 12.0                # W/(m^2 K)
Tini = 20.0             # °C: temperatura inicial de la losa
Tinf_right = 20.0       # °C: aire del lado derecho
Tinf_left_mean = 20.0   # °C: valor medio del aire del lado izquierdo
Tinf_left_amp = 10.0    # °C: amplitud; T_inf,izq varía entre 10 y 30 °C
period = 12.0*3600.0    # s: período de la oscilación

nx = 81
x = np.linspace(0.0, L, nx)
dx = x[1]-x[0]

# Paso estable conservador.
dt_diff = 0.40*dx**2/alpha
dt_conv = 0.40*(rho*cp*dx)/(2*h)
dt = min(dt_diff, dt_conv)

t_end = 48*3600.0
nsteps = int(np.ceil(t_end/dt))
dt = t_end/nsteps
r = alpha*dt/dx**2

T = np.full(nx, Tini)

# Perfiles espaciales: 3 h y 9 h son extremos del ciclo de 12 h;
# 6 h y 12 h corresponden a los cruces por la temperatura media.
snapshot_hours = np.linspace(0., 10., 20)
snapshot_steps = np.rint(snapshot_hours*3600.0/dt).astype(int)
snapshots = {0: T.copy()}

times = []
Tcenter = []
Tsurface_left = []
Tsurface_right = []
Tinf_left_history = []

for n in range(nsteps+1):
    if n % max(1, nsteps//400) == 0:
        t = n*dt
        times.append(t/3600.0)
        Tcenter.append(T[nx//2])
        Tsurface_left.append(T[0])
        Tsurface_right.append(T[-1])
        Tinf_left_history.append(
            Tinf_left_mean + Tinf_left_amp*np.sin(2*np.pi*t/period)
        )

    if n == nsteps:
        break

    Told = T.copy()
    t = n*dt
    Tinf_left = Tinf_left_mean + Tinf_left_amp*np.sin(2*np.pi*t/period)

    T[1:-1] = Told[1:-1] + r*(
        Told[2:] - 2*Told[1:-1] + Told[:-2]
    )

    # Balance sobre medio volumen de control en cada borde.
    T[0] = Told[0] + dt * (
        2*k/(rho*cp*dx**2) * (Told[1]-Told[0])
        + 2*h/(rho*cp*dx) * (Tinf_left-Told[0])
    )

    T[-1] = Told[-1] + dt * (
        2*k/(rho*cp*dx**2) * (Told[-2]-Told[-1])
        + 2*h/(rho*cp*dx) * (Tinf_right-Told[-1])
    )

    # Tras el update, T representa el instante (n+1)·dt.
    if (n + 1) in snapshot_steps:
        snapshots[n + 1] = T.copy()

fig, (ax_history, ax_profile) = plt.subplots(1, 2, figsize=(14, 5.0))

# Historia temporal en puntos representativos.
ax_history.plot(times, Tcenter, lw=2.2, label="centro")
ax_history.plot(times, Tsurface_left, lw=2.0, label="superficie izquierda")
ax_history.plot(times, Tsurface_right, lw=2.0, label="superficie derecha")
ax_history.plot(times, Tinf_left_history, color="black", linestyle="--", lw=1.2,
                label="aire izquierdo: 20 ± 10 °C")
ax_history.axhline(Tinf_right, color="gray", linestyle=":", lw=1.3,
                   label="aire derecho: 20 °C")
ax_history.set_xlabel("tiempo [h]")
ax_history.set_ylabel("T [°C]")
ax_history.set_title("Historia temporal")
# ax_history.legend(fontsize=8)
ax_history.grid(True, alpha=0.35)

# Perfiles dentro de la losa para varios instantes.
colors = plt.cm.viridis(np.linspace(0.05, 0.90, len(snapshot_steps)))
for step, hour, color in zip(snapshot_steps, snapshot_hours, colors):
    ax_profile.plot(x, snapshots[step], color=color, lw=2.0,
                    marker="o", ms=2.8, mfc="white", mew=0.8,
                    markevery=4, label=f"t = {hour:g} h")
ax_profile.set_xlabel("x [m]")
ax_profile.set_ylabel("T [°C]")
ax_profile.set_title("Perfiles de temperatura en la losa")
# ax_profile.legend(title="Perfil", fontsize=8, title_fontsize=9, ncol=2)
ax_profile.grid(True, alpha=0.35)

fig.suptitle(f"1D: Robin con ambiente oscilante ($h={h:.0f}$ W/(m² K))")
fig.tight_layout()
# fig.savefig(FIG/"03_1d_robin.png", dpi=180)
plt.show()
plt.close(fig)
