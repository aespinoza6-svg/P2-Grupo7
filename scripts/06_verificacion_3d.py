"""
Paso 6 — Verificación del solver 2D/3D (mcplan/solver3d.py).

    python scripts/06_verificacion_3d.py

1) fig07: caso seno 2D (bordes Dirichlet). Campo del solver, solución exacta
   y mapa del error. Muestra que el stencil 2D reproduce la solución exacta.

2) fig08: convergencia espacial en 2D y 3D con Fo fijo (dt ∝ dx²). Se espera
   orden 2, igual que en 1D: duplicar las celdas divide el error por 4.

3) fig09: bloque de hormigón de 3 x 3 x 2 m con una fuente de calor que decae
   en el tiempo (ILUSTRATIVA, no es todavía el modelo de hidratación de la
   semana 3), aire por los costados y arriba (Robin) y base aislada. Muestra
   lo que el 1D no puede: núcleo, cara y ARISTA se enfrían distinto.

Salida: resultados/figuras/fig07..fig09 y resultados/convergencia_2d_3d.csv
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mcplan import analiticas as exacta
from mcplan.solver3d import Dirichlet, Malla3D, Material, Robin, resolver

# --- Estilo y salida de las figuras (igual que scripts/04) -------------------
RAIZ = Path(__file__).resolve().parents[1]
CARPETA_FIGURAS = RAIZ / "resultados" / "figuras"
AZUL, NARANJO, VERDE_AGUA, GRIS = "#2a78d6", "#eb6834", "#1baf7a", "#52514e"
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2.0, "savefig.dpi": 160})


def guardar(fig, nombre):
    CARPETA_FIGURAS.mkdir(parents=True, exist_ok=True)
    ruta = CARPETA_FIGURAS / nombre
    fig.savefig(ruta, bbox_inches="tight")
    print(f"  figura -> {ruta.relative_to(RAIZ)}")


DIDACTICO = Material("didactico", k=1.0, rho=1000.0, cp=10.0)
HORMIGON = Material("hormigon", k=2.0, rho=2400.0, cp=900.0)  # config/caso_base.yaml
T_B, A, T_FINAL = 20.0, 40.0, 300.0


def seno(n, dims):
    """Corre el caso seno en un cuadrado/cubo de 1 m con n celdas por lado.
    Devuelve (resultado, solución exacta en el instante final)."""
    forma = (n, n, n) if dims == 3 else (n, n, 1)
    malla = Malla3D.uniforme((1.0, 1.0, 1.0), forma, DIDACTICO)
    X, Y, Z = np.meshgrid(*malla.centros, indexing="ij")
    coords, largos = ((X, Y, Z), (1, 1, 1)) if dims == 3 else ((X, Y), (1, 1))
    caras = ("x-", "x+", "y-", "y+", "z-", "z+")[: 2 * dims]
    bordes = {c: Dirichlet(T_B) for c in caras}
    T0 = exacta.seno_dirichlet_nd(coords, 0.0, largos, DIDACTICO.alpha, T_B, A)
    dt = 0.1 * malla.d[0] ** 2 / DIDACTICO.alpha  # Fo = 0.1 (esquina 3D Dirichlet: <= 1/9)
    res = resolver(malla, T0, bordes, T_FINAL, dt=dt, n_guardados=1)
    return res, exacta.seno_dirichlet_nd(coords, res.t[-1], largos, DIDACTICO.alpha, T_B, A)


def figura_seno_2d():
    res, T_ex = seno(40, dims=2)
    num, ex = res.T[-1][:, :, 0], T_ex[:, :, 0]
    fig, ejes = plt.subplots(1, 3, figsize=(15, 4.3), layout="constrained")
    extension = (0, 1, 0, 1)
    for eje, campo, titulo in ((ejes[0], num, "solver 2D"), (ejes[1], ex, "solución exacta")):
        im = eje.imshow(campo.T, origin="lower", extent=extension, cmap="inferno",
                        vmin=T_B, vmax=num.max())
        eje.set(title=f"{titulo}, t = {T_FINAL:.0f} s", xlabel="x [m]")
        eje.grid(False)
        fig.colorbar(im, ax=eje, label="T [°C]")
    im = ejes[2].imshow((num - ex).T, origin="lower", extent=extension, cmap="RdBu_r",
                        vmin=-np.abs(num - ex).max(), vmax=np.abs(num - ex).max())
    ejes[0].set_ylabel("y [m]")
    ejes[2].set(title=f"error (máx {np.abs(num - ex).max():.1e} °C)", xlabel="x [m]")
    ejes[2].grid(False)
    fig.colorbar(im, ax=ejes[2], label="solver - exacta [°C]")
    fig.suptitle("Verificación 2D: decaimiento del seno con bordes a 20 °C (40 x 40 celdas)")
    guardar(fig, "fig07_seno_2d.png")


def figura_convergencia():
    ns = {2: np.array([10, 20, 40, 80]), 3: np.array([8, 16, 32, 48])}
    filas = []
    fig, eje = plt.subplots(figsize=(6.5, 4.6))
    print("Convergencia espacial (Fo = 0.1)")
    for dims, color in ((2, AZUL), (3, NARANJO)):
        dxs = 1.0 / ns[dims]
        err = np.array([np.max(np.abs(r.T[-1] - ex)) for r, ex in (seno(n, dims) for n in ns[dims])])
        orden = np.log(err[:-1] / err[1:]) / np.log(dxs[:-1] / dxs[1:])
        for n, dx, e, p in zip(ns[dims], dxs, err, [np.nan, *orden]):
            print(f"  {dims}D  N = {n:3d}   dx = {dx:.4f} m   error = {e:.3e} °C   orden = {p:.3f}")
            filas.append(f"{dims}D,{n},{dx},{e}\n")
        eje.loglog(dxs, err, "o-", color=color, ms=6, label=f"{dims}D (orden ≈ {orden[-1]:.2f})")
        eje.loglog(dxs, err[-1] * (dxs / dxs[-1]) ** 2, "--", color=GRIS, lw=1.2)
    eje.plot([], [], "--", color=GRIS, lw=1.2, label="pendiente 2 (referencia)")
    eje.set(xlabel="dx [m]", ylabel="error máximo [°C]",
            title="Convergencia espacial del solver 2D/3D")
    eje.set_xticks([0.01, 0.02, 0.05, 0.1], labels=["0.01", "0.02", "0.05", "0.1"])
    eje.minorticks_off()
    eje.legend(fontsize=9)
    guardar(fig, "fig08_convergencia_2d_3d.png")

    csv = RAIZ / "resultados" / "convergencia_2d_3d.csv"
    with open(csv, "w", encoding="utf-8") as f:
        f.write("dimension,N,dx_m,error_C\n")
        f.writelines(filas)
    print(f"  tabla  -> {csv.relative_to(RAIZ)}")


def figura_bloque_hormigon():
    largos, n = (3.0, 3.0, 2.0), (30, 30, 20)          # celdas de 10 cm
    malla = Malla3D.uniforme(largos, n, HORMIGON)
    aire = Robin(h=10.0, T_inf=15.0)
    bordes = {"x-": aire, "x+": aire, "y-": aire, "y+": aire, "z+": aire}  # z- aislada
    Q0, tau = 1500.0, 1.0 * 86400.0                     # W/m³ y 1 día (ILUSTRATIVO)
    res = resolver(malla, 18.0, bordes, 7 * 86400.0, n_guardados=168,
                   fuente=lambda t, T: Q0 * np.exp(-t / tau))
    dias = res.t / 86400.0
    i_c, j_c = n[0] // 2, n[1] // 2
    nucleo = res.T[:, i_c, j_c, 0]                       # centro de la base aislada
    cara = res.T[:, 0, j_c, 0]                           # celda junto a la cara x-
    arista = res.T[:, 0, 0, -1]                          # esquina superior
    print(f"Bloque: dt = {res.dt:.0f} s, Fo = {res.fourier:.3f}, "
          f"núcleo máx {nucleo.max():.1f} °C, arista máx {arista.max():.1f} °C")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    a1.plot(dias, nucleo, color=NARANJO, label="núcleo (centro, abajo)")
    a1.plot(dias, cara, color=AZUL, label="junto a una cara lateral")
    a1.plot(dias, arista, color=VERDE_AGUA, label="esquina superior")
    a1.axhline(15.0, color=GRIS, ls="--", lw=1.2, label="aire 15 °C")
    a1.set(xlabel="tiempo [días]", ylabel="T [°C]", title="Tres puntos del bloque")
    a1.legend(fontsize=8)
    k_max = int(np.argmax(nucleo))
    x, _, z = malla.centros
    im = a2.imshow(res.T[k_max][:, j_c, :].T, origin="lower", cmap="inferno",
                   extent=(0, largos[0], 0, largos[2]), aspect="equal")
    a2.set(title=f"corte vertical y = {largos[1] / 2:.1f} m, día {dias[k_max]:.1f}",
           xlabel="x [m]", ylabel="z [m]")
    a2.grid(False)
    fig.colorbar(im, ax=a2, label="T [°C]")
    fig.suptitle("Bloque de hormigón 3 x 3 x 2 m con fuente de calor ilustrativa "
                 "(Robin h = 10 W/m²K, base aislada)")
    guardar(fig, "fig09_bloque_hormigon_3d.png")


def main():
    figura_seno_2d()
    figura_convergencia()
    figura_bloque_hormigon()


if __name__ == "__main__":
    main()
