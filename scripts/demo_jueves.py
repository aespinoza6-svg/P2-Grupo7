import numpy as np
from mcplan.solver1d import Material, Malla1D, Dirichlet, resolver

N = 40  # cantidad de celdas: este será el valor que cambiaremos

material = Material("Ejemplo", k=1.0, rho=1000.0, cp=10.0)
malla = Malla1D.uniforme(L=1.0, n=N, material=material)

T_inicial = 20 + 50 * np.sin(np.pi * malla.x)
dt = 0.25 * malla.dx**2 / material.alpha

resultado = resolver(
    malla, T_inicial, Dirichlet(20), Dirichlet(20),
    t_final=1500, dt=dt,
)

T_exacta = (
    20 + 50 * np.sin(np.pi * resultado.x)
    * np.exp(-material.alpha * np.pi**2 * resultado.t[-1])
)
error = np.max(np.abs(resultado.T[-1] - T_exacta))

print(f"Celdas: {N}")
print(f"Delta x: {malla.dx:.4f} m")
print(f"Delta t usado: {resultado.dt:.4f} s")
print(f"Pasos de tiempo: {resultado.n_pasos}")
print(f"Error máximo: {error:.6f} °C")