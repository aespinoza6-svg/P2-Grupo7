"""
Solver térmico 3D (y 2D) explícito por balance de energía en celdas.

Es el MISMO método que mcplan/solver1d.py, repetido en las tres direcciones.
Ecuación que se resuelve:

    rho * cp * dT/dt = div( k * grad T ) + Q

    (en 3D: d/dx(k dT/dx) + d/dy(k dT/dy) + d/dz(k dT/dz) + Q)

La deducción está en docs/1_teoria.md (Parte D). Resumen:

1) MALLA CARTESIANA DE CELDAS (voxeles)
   La caja [0, Lx] x [0, Ly] x [0, Lz] se divide en nx·ny·nz celdas de tamaño
   dx·dy·dz. T[i, j, l] es la temperatura de la celda (i, j, l): i avanza en x,
   j en y, l en z (z hacia arriba). Cada celda tiene su propio material, así
   más adelante la voxelización del IFC solo tiene que decir "qué material
   hay en cada celda".

   2D = una sola capa de celdas en z (nz = 1) con las caras z aisladas.

2) BALANCE DE ENERGÍA EN LA CELDA (acumulación = entra - sale + genera)

        C · dT/dt = (Fx_izq - Fx_der) + (Fy_izq - Fy_der) + (Fz_abajo - Fz_arriba) + Q·V

        V = dx·dy·dz                 volumen de la celda          [m³]
        C = rho·cp·V                 capacidad térmica             [J/K]
        F = flujo por una cara, POSITIVO hacia +x, +y o +z          [W]

   A diferencia del 1D (que trabajaba "por m²"), aquí los flujos están en W
   porque las caras de x, y, z tienen áreas distintas:
        cara x: A = dy·dz     cara y: A = dx·dz     cara z: A = dx·dy

3) FLUJO ENTRE VECINAS = conductancia × diferencia de temperatura

        F = G · (T_celda - T_vecina)
        G = A / ( d/(2·k_celda) + d/(2·k_vecina) )      (dos medias celdas en serie)

   d es dx, dy o dz según la dirección. Si las dos celdas son del mismo
   material, G = k·A/d. Si no (hormigón-suelo, vertida nueva-vieja), la
   fórmula da automáticamente la media armónica de k (Parte D de la teoría).

4) BORDES: cada una de las 6 caras exteriores de la caja tiene su condición,
   con las mismas clases del 1D: Dirichlet, Neumann o Robin. Una cara que no
   se indica queda AISLADA (Neumann(0)).

        "x-" (x = 0)   "x+" (x = Lx)
        "y-" (y = 0)   "y+" (y = Ly)
        "z-" (z = 0, abajo)   "z+" (z = Lz, arriba)

5) AVANCE EN EL TIEMPO: Euler explícito, estable si dt <= C/S en todas las
   celdas (S = suma de conductancias de la celda). En malla uniforme y celda
   interior: Fo = alpha·dt/dx² <= 1/4 en 2D y <= 1/6 en 3D.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence, Union

import numpy as np

from mcplan.solver1d import (
    Borde,
    Dirichlet,
    Material,
    Neumann,
    Robin,
    _conductancia_borde,
    _flujo_entrante,
)

__all__ = ["Borde", "Dirichlet", "Material", "Neumann", "Robin",
           "Malla3D", "Resultado3D", "CARAS", "paso_tiempo_maximo", "resolver"]

# La fuente Q [W/m³] puede ser constante, un arreglo (un valor por celda) o una
# función Q(t, T) del tiempo y del campo de temperaturas actual. Que pueda
# depender de T es lo que necesitará el calor de hidratación (semana 3).
Fuente = Union[float, np.ndarray, Callable[[float, np.ndarray], np.ndarray]]

# Nombres de las 6 caras exteriores: (eje, lado). lado 0 = inicio, 1 = final.
CARAS = {"x-": (0, 0), "x+": (0, 1),
         "y-": (1, 0), "y+": (1, 1),
         "z-": (2, 0), "z+": (2, 1)}


# ---------------------------------------------------------------------------
# Malla
# ---------------------------------------------------------------------------

@dataclass(eq=False)
class Malla3D:
    """Malla cartesiana de nx·ny·nz celdas con propiedades por celda.

    Crear con Malla3D.uniforme(...) (un material) o Malla3D.por_materiales(...)
    (un número de material por celda, p. ej. hormigón sobre suelo).
    """

    d: tuple[float, float, float]  # (dx, dy, dz) [m]
    k: np.ndarray    # conductividad de cada celda, forma (nx, ny, nz)
    rho: np.ndarray  # densidad de cada celda
    cp: np.ndarray   # calor específico de cada celda

    @classmethod
    def uniforme(cls, largos: Sequence[float], n: Sequence[int],
                 material: Material) -> "Malla3D":
        """Caja de largos (Lx, Ly, Lz) [m] de un solo material, con
        n = (nx, ny, nz) celdas. Para 2D usar nz = 1."""
        return cls.por_materiales(largos, np.zeros(tuple(n), dtype=int), [material])

    @classmethod
    def por_materiales(cls, largos: Sequence[float], indice: np.ndarray,
                       materiales: Sequence[Material]) -> "Malla3D":
        """Caja de largos (Lx, Ly, Lz) [m] donde indice[i, j, l] dice qué
        material (posición en la lista materiales) ocupa cada celda."""
        indice = np.asarray(indice, dtype=int)
        if indice.ndim != 3:
            raise ValueError("indice debe tener forma (nx, ny, nz); para 2D use nz = 1.")
        d = tuple(float(L) / n for L, n in zip(largos, indice.shape))
        k = np.array([m.k for m in materiales], dtype=float)[indice]
        rho = np.array([m.rho for m in materiales], dtype=float)[indice]
        cp = np.array([m.cp for m in materiales], dtype=float)[indice]
        return cls(d=d, k=k, rho=rho, cp=cp)

    @property
    def forma(self) -> tuple[int, int, int]:
        return self.k.shape

    @property
    def n_celdas(self) -> int:
        return self.k.size

    @property
    def centros(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Coordenadas de los centros de celda en x, y, z [m] (tres vectores)."""
        return tuple((np.arange(n) + 0.5) * d for n, d in zip(self.forma, self.d))

    @property
    def volumen_celda(self) -> float:
        dx, dy, dz = self.d
        return dx * dy * dz

    @property
    def capacidad(self) -> np.ndarray:
        """C = rho·cp·V [J/K]: energía para subir 1 K cada celda."""
        return self.rho * self.cp * self.volumen_celda

    @property
    def alpha(self) -> np.ndarray:
        return self.k / (self.rho * self.cp)

    def area_cara(self, eje: int) -> float:
        """Área [m²] de una cara perpendicular al eje (0 = x, 1 = y, 2 = z)."""
        dx, dy, dz = self.d
        return (dy * dz, dx * dz, dx * dy)[eje]

    def conductancia_interior(self, eje: int) -> np.ndarray:
        """G [W/K] de las caras interiores perpendiculares a un eje.

        Forma: la de la malla con una celda menos en ese eje (hay n-1 caras
        interiores entre n celdas). Dos medias celdas en serie.
        """
        r_media = self.d[eje] / (2.0 * self.k)          # [m² K/W] por celda
        r_a = np.take(r_media, range(self.forma[eje] - 1), axis=eje)
        r_b = np.take(r_media, range(1, self.forma[eje]), axis=eje)
        return self.area_cara(eje) / (r_a + r_b)


def _celdas_de_cara(arreglo: np.ndarray, eje: int, lado: int) -> np.ndarray:
    """La capa de celdas que toca la cara exterior (eje, lado)."""
    return np.take(arreglo, 0 if lado == 0 else -1, axis=eje)


def _completar_bordes(bordes: Mapping[str, Borde] | None) -> dict[str, Borde]:
    """Las 6 caras con su condición; las que no se indican quedan aisladas."""
    bordes = dict(bordes or {})
    desconocidas = set(bordes) - set(CARAS)
    if desconocidas:
        raise KeyError(f"Caras desconocidas {sorted(desconocidas)}; use {list(CARAS)}.")
    return {cara: bordes.get(cara, Neumann(0.0)) for cara in CARAS}


def _conductancias_borde(malla: Malla3D, bordes: dict[str, Borde]) -> dict[str, np.ndarray]:
    """G [W/K] entre cada celda de una cara exterior y su referencia exterior
    (superficie en Dirichlet, aire en Robin; 0 en Neumann)."""
    G = {}
    for cara, (eje, lado) in CARAS.items():
        k_cara = _celdas_de_cara(malla.k, eje, lado)
        G_m2 = _conductancia_borde(bordes[cara], k_cara, malla.d[eje])  # [W/(m² K)]
        G[cara] = np.broadcast_to(G_m2 * malla.area_cara(eje), k_cara.shape)
    return G


# ---------------------------------------------------------------------------
# Estabilidad
# ---------------------------------------------------------------------------

def paso_tiempo_maximo(malla: Malla3D, bordes: Mapping[str, Borde] | None = None) -> float:
    """Mayor dt [s] estable: dt <= C/S en toda celda (mismo criterio del 1D).

    S = suma de las conductancias de la celda (vecinas en x, y, z + bordes).
    Malla uniforme, celda interior:
        2D: S = 4kA/d, C = rho·cp·V  ->  Fo <= 1/4
        3D: S = 6kA/d               ->  Fo <= 1/6
    """
    bordes = _completar_bordes(bordes)
    S = np.zeros(malla.forma)
    for eje in range(3):
        G = malla.conductancia_interior(eje)
        n = malla.forma[eje]
        # cada cara interior suma a la celda de cada lado
        S[(slice(None),) * eje + (slice(0, n - 1),)] += G
        S[(slice(None),) * eje + (slice(1, n),)] += G
    for cara, G in _conductancias_borde(malla, bordes).items():
        eje, lado = CARAS[cara]
        S[(slice(None),) * eje + (0 if lado == 0 else -1,)] += G
    with np.errstate(divide="ignore"):
        return float(np.min(malla.capacidad / S))  # S = 0 (todo aislado y 1 celda) -> inf


# ---------------------------------------------------------------------------
# Resultado y solver
# ---------------------------------------------------------------------------

@dataclass(eq=False)
class Resultado3D:
    """Salida de resolver(): campos T guardados en algunos instantes."""

    t: np.ndarray       # instantes guardados [s], forma (M,)
    T: np.ndarray       # temperaturas [°C], forma (M, nx, ny, nz)
    dt: float           # paso de tiempo usado [s]
    n_pasos: int
    dt_max: float       # paso máximo estable [s]
    malla: Malla3D

    @property
    def centros(self):
        return self.malla.centros

    @property
    def fourier(self) -> float:
        """Fo = alpha·dt/d² con el menor tamaño de celda de los ejes con más
        de una celda (el que limita la estabilidad)."""
        d_min = min(d for d, n in zip(self.malla.d, self.malla.forma) if n > 1)
        return float(np.max(self.malla.alpha) * self.dt / d_min**2)

    def energia(self) -> np.ndarray:
        """Energía interna (referida a 0 °C) en cada instante: E = sum C·T [J]."""
        return np.tensordot(self.T, self.malla.capacidad, axes=3)


def resolver(
    malla: Malla3D,
    T_inicial: float | np.ndarray,
    bordes: Mapping[str, Borde] | None,
    t_final: float,
    *,
    dt: float | None = None,
    fuente: Fuente = 0.0,
    n_guardados: int = 20,
    tiempos_guardado: Sequence[float] | None = None,
    permitir_inestable: bool = False,
) -> Resultado3D:
    """Integra la conducción 3D desde t = 0 hasta t_final [s].

    Parámetros
    ----------
    malla        : Malla3D con geometría y propiedades.
    T_inicial    : condición inicial [°C]: un número o un arreglo (nx, ny, nz).
    bordes       : dict {"x-": Robin(...), "z-": Dirichlet(...), ...}.
                   Las caras que no aparecen quedan aisladas (Neumann(0)).
    t_final      : duración de la simulación [s].
    dt           : paso [s]; None = 90 % del máximo estable.
    fuente       : Q [W/m³]: número, arreglo por celda o función Q(t, T).
    n_guardados / tiempos_guardado / permitir_inestable: igual que en 1D.
    """
    bordes = _completar_bordes(bordes)
    T = np.array(np.broadcast_to(np.asarray(T_inicial, dtype=float), malla.forma))

    # --- Paso de tiempo (igual que en 1D) -----------------------------------
    dt_max = paso_tiempo_maximo(malla, bordes)
    if dt is None:
        dt = 0.9 * dt_max
    n_pasos = int(np.ceil(t_final / dt - 1e-9))
    dt = t_final / n_pasos
    if dt > dt_max * (1 + 1e-12) and not permitir_inestable:
        raise ValueError(
            f"dt = {dt:.4g} s supera el paso máximo estable {dt_max:.4g} s. "
            "Reduzca dt o use una malla más gruesa."
        )

    if tiempos_guardado is None:
        pasos = np.linspace(0, n_pasos, n_guardados + 1)
    else:
        pasos = np.asarray(tiempos_guardado, dtype=float) / dt
    pasos_guardado = np.unique(np.clip(np.rint(pasos).astype(int), 0, n_pasos))

    # --- Cantidades que no cambian en el tiempo -----------------------------
    C = malla.capacidad
    G_int = [malla.conductancia_interior(eje) for eje in range(3)]
    G_borde = _conductancias_borde(malla, bordes)
    # Flujo en cada cara [W], positivo hacia + del eje. En el eje x hay nx+1
    # caras (nx-1 interiores + 2 de borde); igual en y, z.
    F = [np.zeros(tuple(n + (1 if e == eje else 0) for e, n in enumerate(malla.forma)))
         for eje in range(3)]

    t_guardados: list[float] = []
    T_guardados: list[np.ndarray] = []
    siguiente = 0

    # --- Bucle temporal ------------------------------------------------------
    for paso in range(n_pasos + 1):
        t = paso * dt
        if siguiente < len(pasos_guardado) and paso == pasos_guardado[siguiente]:
            t_guardados.append(t)
            T_guardados.append(T.copy())
            siguiente += 1
        if paso == n_pasos:
            break

        balance = np.zeros(malla.forma)  # calor neto que entra a cada celda [W]
        for eje in range(3):
            n = malla.forma[eje]
            sl = lambda a, b: (slice(None),) * eje + (slice(a, b),)
            # 1) Caras interiores: del más caliente al más frío.
            F[eje][sl(1, n)] = G_int[eje] * (T[sl(0, n - 1)] - T[sl(1, n)])
            # 2) Caras de borde: por el inicio lo que entra va hacia +,
            #    por el final lo que entra va hacia -  (igual que en 1D).
            for lado, signo, cara_idx in ((0, +1, 0), (1, -1, n)):
                cara = "xyz"[eje] + ("-" if lado == 0 else "+")
                T_cara = _celdas_de_cara(T, eje, lado)
                q = _flujo_entrante(bordes[cara], G_borde[cara], T_cara, t)
                if isinstance(bordes[cara], Neumann):
                    q = q * malla.area_cara(eje)  # Neumann viene en W/m²
                F[eje][(slice(None),) * eje + (cara_idx,)] = signo * q
            # 3) Lo que entra por la cara de atrás menos lo que sale por la de adelante.
            balance += F[eje][sl(0, n)] - F[eje][sl(1, n + 1)]

        # 4) Fuente volumétrica y avance explícito.
        Q = fuente(t, T) if callable(fuente) else fuente
        T = T + dt * (balance + Q * malla.volumen_celda) / C

    return Resultado3D(
        t=np.array(t_guardados),
        T=np.array(T_guardados),
        dt=dt,
        n_pasos=n_pasos,
        dt_max=dt_max,
        malla=malla,
    )
