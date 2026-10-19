"""
Solver térmico 1D explícito por balance de energía en celdas.

Ecuación que se resuelve (conducción de calor con fuente interna):

    rho * cp * dT/dt = d/dx( k * dT/dx ) + Q          0 < x < L,  t > 0

    T    temperatura                     [°C]
    k    conductividad térmica           [W/(m K)]
    rho  densidad                        [kg/m³]
    cp   calor específico                [J/(kg K)]
    Q    fuente volumétrica de calor     [W/m³]  (más adelante: hidratación)

La deducción completa está en docs/1_teoria.md (Parte C). Resumen con la misma
notación que usa el código:

1) MALLA DE CELDAS (volúmenes finitos centrados en celdas)
   [0, L] se divide en N celdas de ancho dx = L/N. La incógnita T[i] es la
   temperatura de la celda i, cuyo centro está en x[i] = (i + 1/2)·dx.
   Trabajamos por unidad de área transversal (A = 1 m²): capacidades y
   conductancias quedan "por m²".

        cara 0     cara 1     cara 2               cara N
          | celda 0  | celda 1  |   ...   | celda N-1 |
         x=0                                         x=L

2) BALANCE DE ENERGÍA EN LA CELDA i  (acumulación = entra - sale + genera)

        C[i] · dT[i]/dt = F[i] - F[i+1] + Q[i]·dx

        C[i] = rho[i]·cp[i]·dx    capacidad térmica de la celda   [J/(m² K)]
        F[j] = flujo de calor por la cara j, POSITIVO hacia +x     [W/m²]

3) FLUJO ENTRE VECINAS = conductancia × diferencia de temperatura

        F[i+1] = G[i+1] · (T[i] - T[i+1])        (el calor va de caliente a frío)
        G = 1/R,  R = dx/(2·k[i]) + dx/(2·k[i+1])   (dos medias celdas en serie)

   Si ambas celdas tienen el mismo k, G = k/dx: la diferencia finita clásica.

4) BORDES (caras 0 y N): Dirichlet, Neumann o Robin (clases más abajo).

5) AVANCE EN EL TIEMPO: Euler explícito (FTCS)

        T[i]^(n+1) = T[i]^n + dt/C[i] · ( F[i]^n - F[i+1]^n + Q[i]^n·dx )

   Es estable si dt <= C[i] / S[i] en TODAS las celdas, donde S[i] es la suma
   de las conductancias de la celda i (ver paso_tiempo_maximo). En una malla
   uniforme y lejos de los bordes esto equivale a Fo = alpha·dt/dx² <= 1/2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence, Union

import numpy as np

# Un valor de borde puede ser constante o una función del tiempo t [s].
# Ej.: T_inf = 15.0   o   T_inf = lambda t: 15 + 8*np.sin(2*np.pi*t/86400)
Valor = Union[float, Callable[[float], float]]

# La fuente Q [W/m³] puede ser constante o una función Q(t, x) que devuelve
# un arreglo con un valor por celda.
Fuente = Union[float, Callable[[float, np.ndarray], np.ndarray]]


def _evaluar(valor: Valor, t: float) -> float:
    """Valor en el instante t (si es una constante, se devuelve tal cual)."""
    return float(valor(t)) if callable(valor) else float(valor)


# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Material:
    """Propiedades térmicas de un material homogéneo e isótropo."""

    nombre: str
    k: float    # conductividad térmica [W/(m K)]: qué tan fácil conduce calor
    rho: float  # densidad [kg/m³]
    cp: float   # calor específico [J/(kg K)]: energía para subir 1 K un kg

    @property
    def alpha(self) -> float:
        """Difusividad térmica alpha = k/(rho·cp) [m²/s]: qué tan rápido se
        propaga un cambio de temperatura (conduce / almacena)."""
        return self.k / (self.rho * self.cp)


# ---------------------------------------------------------------------------
# Condiciones de borde
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Dirichlet:
    """Temperatura impuesta en la superficie: T(borde, t) = T(t).

    Física: contacto con algo de temperatura conocida y resistencia
    despreciable. La celda del borde "ve" la superficie a dx/2 de distancia,
    por eso su conductancia es 2k/dx.
    """

    T: Valor


@dataclass(frozen=True)
class Neumann:
    """Flujo de calor impuesto que ENTRA al dominio: q(t) [W/m²].

    q = 0 -> borde adiabático (aislación perfecta o plano de simetría).
    q > 0 -> entra calor (p. ej. radiación solar absorbida).
    El flujo no depende de la temperatura de la celda.
    """

    q: Valor = 0.0


@dataclass(frozen=True)
class Robin:
    """Convección con el ambiente (ley de enfriamiento de Newton):

        q_entrante = (T_inf - T_superficie) / (1/h + R_extra)

    h        coeficiente de convección [W/(m² K)] (crece con el viento)
    T_inf    temperatura del aire [°C], constante o función del tiempo
    R_extra  resistencia térmica adicional en serie [m² K/W], p. ej.
             moldaje o aislación: R = espesor/k_material. 0 = superficie desnuda.
    """

    h: float
    T_inf: Valor
    R_extra: float = 0.0


Borde = Union[Dirichlet, Neumann, Robin]


def _conductancia_borde(borde: Borde, k_celda: float, dx: float) -> float:
    """Conductancia [W/(m² K)] entre el centro de la celda del borde y la
    referencia exterior (superficie en Dirichlet, aire en Robin).

    Las resistencias en serie se suman:
        media celda (dx/2k) [+ moldaje/aislación R_extra + película de aire 1/h].
    Neumann devuelve 0 porque su flujo no depende de la temperatura.
    """
    r_media_celda = dx / (2.0 * k_celda)
    if isinstance(borde, Dirichlet):
        return 1.0 / r_media_celda
    if isinstance(borde, Robin):
        return 1.0 / (r_media_celda + borde.R_extra + 1.0 / borde.h)
    if isinstance(borde, Neumann):
        return 0.0
    raise TypeError(f"Condición de borde desconocida: {borde!r}")


def _flujo_entrante(borde: Borde, G_borde: float, T_celda: float, t: float) -> float:
    """Flujo de calor [W/m²] que ENTRA al dominio por ese borde en el instante t."""
    if isinstance(borde, Dirichlet):
        return G_borde * (_evaluar(borde.T, t) - T_celda)
    if isinstance(borde, Robin):
        return G_borde * (_evaluar(borde.T_inf, t) - T_celda)
    return _evaluar(borde.q, t)  # Neumann


# ---------------------------------------------------------------------------
# Malla
# ---------------------------------------------------------------------------

@dataclass(eq=False)
class Malla1D:
    """Malla 1D de N celdas de igual ancho dx, con propiedades por celda.

    Crear con Malla1D.uniforme(...) (un material) o Malla1D.por_capas(...)
    (varios materiales en serie, p. ej. hormigón sobre suelo).
    """

    dx: float
    k: np.ndarray    # conductividad de cada celda
    rho: np.ndarray  # densidad de cada celda
    cp: np.ndarray   # calor específico de cada celda

    @classmethod
    def uniforme(cls, L: float, n: int, material: Material) -> "Malla1D":
        """Barra de largo L [m] de un solo material, dividida en n celdas."""
        return cls.por_capas([(L, material)], n)

    @classmethod
    def por_capas(cls, capas: Sequence[tuple[float, Material]], n: int) -> "Malla1D":
        """Capas en serie [(espesor_m, Material), ...] ordenadas desde x = 0.

        Exigimos que cada interfaz entre materiales caiga justo en una cara
        de celda (número entero de celdas por capa). Así cada celda tiene un
        solo material y la interfaz se trata con resistencias en serie.
        """
        if n < 2:
            raise ValueError("Se necesitan al menos 2 celdas.")
        espesores = np.array([e for e, _ in capas], dtype=float)
        L = espesores.sum()
        dx = L / n
        celdas_por_capa = espesores / dx
        if not np.allclose(celdas_por_capa, np.round(celdas_por_capa), atol=1e-9):
            raise ValueError(
                f"Con n={n} (dx={dx:.4g} m) las capas no tienen un número entero "
                f"de celdas ({celdas_por_capa}). Elija otro n."
            )
        repeticiones = np.round(celdas_por_capa).astype(int)
        materiales = [m for _, m in capas]
        k = np.repeat([m.k for m in materiales], repeticiones)
        rho = np.repeat([m.rho for m in materiales], repeticiones)
        cp = np.repeat([m.cp for m in materiales], repeticiones)
        return cls(dx=dx, k=k, rho=rho, cp=cp)

    @property
    def n(self) -> int:
        return len(self.k)

    @property
    def L(self) -> float:
        return self.n * self.dx

    @property
    def x(self) -> np.ndarray:
        """Posición del centro de cada celda [m]."""
        return (np.arange(self.n) + 0.5) * self.dx

    @property
    def capacidad(self) -> np.ndarray:
        """C[i] = rho·cp·dx [J/(m² K)]: energía para subir 1 K la celda (por m²)."""
        return self.rho * self.cp * self.dx

    @property
    def conductancia_interior(self) -> np.ndarray:
        """G de las N-1 caras interiores [W/(m² K)]: dos medias celdas en serie."""
        r_media = self.dx / (2.0 * self.k)
        return 1.0 / (r_media[:-1] + r_media[1:])

    @property
    def alpha(self) -> np.ndarray:
        return self.k / (self.rho * self.cp)


# ---------------------------------------------------------------------------
# Estabilidad
# ---------------------------------------------------------------------------

def paso_tiempo_maximo(malla: Malla1D, borde_izq: Borde, borde_der: Borde) -> float:
    """Mayor dt [s] estable para el esquema explícito.

    Reordenando la actualización explícita de la celda i:

        T[i]^(n+1) = (1 - dt·S[i]/C[i])·T[i]^n + dt/C[i]·(sum_vecinos G·T_vecino + fuentes)

    con S[i] = suma de las conductancias de la celda i. Si el coeficiente
    (1 - dt·S/C) fuera negativo, la nueva temperatura dependería al revés de
    la anterior: aparecen oscilaciones que crecen (inestabilidad). Exigimos
    entonces dt <= C[i]/S[i] en toda celda (criterio de coeficientes positivos).

    Malla uniforme, celda interior: S = 2k/dx, C = rho·cp·dx  ->  Fo <= 1/2.
    Celda junto a un borde Dirichlet: S = k/dx + 2k/dx      ->  Fo <= 1/3.
    """
    G = malla.conductancia_interior
    S = np.zeros(malla.n)
    S[:-1] += G  # cara derecha de las celdas 0..N-2
    S[1:] += G   # cara izquierda de las celdas 1..N-1
    S[0] += _conductancia_borde(borde_izq, malla.k[0], malla.dx)
    S[-1] += _conductancia_borde(borde_der, malla.k[-1], malla.dx)
    return float(np.min(malla.capacidad / S))


# ---------------------------------------------------------------------------
# Resultado y solver
# ---------------------------------------------------------------------------

@dataclass(eq=False)
class Resultado1D:
    """Salida de resolver(): perfiles T(x) guardados en algunos instantes."""

    x: np.ndarray       # centros de celda [m], forma (N,)
    t: np.ndarray       # instantes guardados [s], forma (M,)
    T: np.ndarray       # temperaturas [°C], forma (M, N): fila = instante
    dt: float           # paso de tiempo efectivamente usado [s]
    n_pasos: int        # número de pasos de tiempo
    dt_max: float       # paso máximo estable de esta malla y bordes [s]
    malla: Malla1D

    @property
    def fourier(self) -> float:
        """Número de Fourier Fo = alpha·dt/dx² (el mayor de la malla)."""
        return float(np.max(self.malla.alpha) * self.dt / self.malla.dx**2)

    def energia(self) -> np.ndarray:
        """Energía interna por m² (referida a 0 °C) en cada instante guardado:
        E = sum_i C[i]·T[i]  [J/m²]."""
        return self.T @ self.malla.capacidad


def resolver(
    malla: Malla1D,
    T_inicial: float | np.ndarray,
    borde_izq: Borde,
    borde_der: Borde,
    t_final: float,
    *,
    dt: float | None = None,
    fuente: Fuente = 0.0,
    n_guardados: int = 50,
    tiempos_guardado: Sequence[float] | None = None,
    permitir_inestable: bool = False,
) -> Resultado1D:
    """Integra la conducción 1D desde t = 0 hasta t_final [s].

    Parámetros
    ----------
    malla        : Malla1D con geometría y propiedades.
    T_inicial    : condición inicial [°C]: un número (uniforme) o N valores.
    borde_izq/der: Dirichlet, Neumann o Robin en x = 0 y x = L.
    t_final      : duración de la simulación [s].
    dt           : paso de tiempo [s]. Si es None se usa el 90 % del máximo
                   estable. Se reduce levemente para que n_pasos·dt = t_final.
    fuente       : Q [W/m³], constante o función Q(t, x).
    n_guardados  : cuántos perfiles guardar (equiespaciados), si no se da
                   tiempos_guardado.
    tiempos_guardado  : instantes [s] en que se quieren los perfiles.
    permitir_inestable: True solo para DEMOSTRAR la inestabilidad (lanza
                   error por defecto si dt supera el máximo estable).
    """
    T = np.array(np.broadcast_to(np.asarray(T_inicial, dtype=float), (malla.n,)))

    # --- Paso de tiempo ----------------------------------------------------
    dt_max = paso_tiempo_maximo(malla, borde_izq, borde_der)
    if dt is None:
        dt = 0.9 * dt_max
    n_pasos = int(np.ceil(t_final / dt - 1e-9))
    dt = t_final / n_pasos  # ajuste para terminar exactamente en t_final
    if dt > dt_max * (1 + 1e-12) and not permitir_inestable:
        raise ValueError(
            f"dt = {dt:.4g} s supera el paso máximo estable {dt_max:.4g} s "
            f"(Fo = {np.max(malla.alpha) * dt / malla.dx**2:.3f}). "
            "Reduzca dt o use una malla más gruesa."
        )

    # --- Instantes a guardar (convertidos a número de paso) -----------------
    if tiempos_guardado is None:
        pasos = np.linspace(0, n_pasos, n_guardados + 1)
    else:
        pasos = np.asarray(tiempos_guardado, dtype=float) / dt
    pasos_guardado = np.unique(np.clip(np.rint(pasos).astype(int), 0, n_pasos))

    # --- Cantidades que no cambian en el tiempo -----------------------------
    C = malla.capacidad
    G = malla.conductancia_interior
    x = malla.x
    G_izq = _conductancia_borde(borde_izq, malla.k[0], malla.dx)
    G_der = _conductancia_borde(borde_der, malla.k[-1], malla.dx)
    F = np.empty(malla.n + 1)  # flujo en cada una de las N+1 caras [W/m²]

    t_guardados: list[float] = []
    T_guardados: list[np.ndarray] = []
    siguiente = 0  # índice del próximo instante a guardar

    # --- Bucle temporal ------------------------------------------------------
    for paso in range(n_pasos + 1):
        t = paso * dt
        if siguiente < len(pasos_guardado) and paso == pasos_guardado[siguiente]:
            t_guardados.append(t)
            T_guardados.append(T.copy())
            siguiente += 1
        if paso == n_pasos:
            break

        # 1) Flujos en caras interiores: del más caliente al más frío.
        F[1:-1] = G * (T[:-1] - T[1:])
        # 2) Flujos en los bordes (F es positivo hacia +x):
        #    por la izquierda lo que entra va hacia +x  -> F[0]  = +q_entrante
        #    por la derecha  lo que entra va hacia -x  -> F[-1] = -q_entrante
        F[0] = _flujo_entrante(borde_izq, G_izq, T[0], t)
        F[-1] = -_flujo_entrante(borde_der, G_der, T[-1], t)
        # 3) Fuente volumétrica en t^n.
        Q = fuente(t, x) if callable(fuente) else fuente
        # 4) Balance de cada celda y avance explícito.
        T = T + dt * (F[:-1] - F[1:] + Q * malla.dx) / C

    return Resultado1D(
        x=x,
        t=np.array(t_guardados),
        T=np.array(T_guardados),
        dt=dt,
        n_pasos=n_pasos,
        dt_max=dt_max,
        malla=malla,
    )
