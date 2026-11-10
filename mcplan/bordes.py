"""Condiciones de borde ambientales: como el clima real entra a la fisica.

Todas terminan en la misma forma por cara expuesta de una celda:

    q_entrante [W] = G(t) * (T_ref(t) - T_celda) + Q(t)

- Temperatura del aire: es T_ref del termino Robin (conveccion).
- Viento: cambia h. h_conv = 5,7 + 3,8 V [W/m2K], con V a 10 m en m/s
  (correlacion clasica de McAdams/Jurges; aproximada).
- Radiacion de onda larga: se linealiza como h_rad = 4 eps sigma T^3 y se SUMA a h,
  suponiendo el cielo a la temperatura del aire. Con V = 0 da h ~ 10,6 W/m2K.
- Sol: el flujo absorbido alpha*G entra como temperatura sol-aire,
  T_sa = T_aire + alpha*G/h. Asi el sol y la conveccion no se cuentan dos veces:
  h solo multiplica diferencias de temperatura y el sol es un aporte aparte.
  Funciona tambien bajo una aislacion: el sol calienta la cara externa de la manta.
- Moldaje o aislacion: resistencia R [m2K/W] en serie con la pelicula de aire:
  G = A / (dx/2k + R + 1/h). Al desmoldar, R pasa a 0 en la hora de retiro de
  CADA vertida: la superficie queda de golpe expuesta al aire.

El conductor G se evalua al inicio de cada paso (Euler explicito); el paso estable
se calcula con la mayor conductancia posible del horizonte (R = 0 y viento maximo).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from mcplan.solver3d import CARAS

SIGMA = 5.670374419e-8  # W/(m2 K4), Stefan-Boltzmann


def h_conveccion(viento_m_s):
    """Coeficiente convectivo exterior [W/m2K] para viento V [m/s] (V >= 0)."""
    return 5.7 + 3.8 * np.maximum(viento_m_s, 0.0)


def h_radiacion(T_C, emisividad=0.9):
    """Radiacion de onda larga linealizada [W/m2K] alrededor de T [degC]."""
    return 4.0 * emisividad * SIGMA * (np.asarray(T_C, dtype=float) + 273.15) ** 3


@dataclass
class BordeAmbiental:
    """Aire + viento + sol + moldaje/aislacion, por celda y en el tiempo.

    clima: objeto con .en(t) -> (T_aire, viento, ghi) y series .T_aire_C, .viento_m_s.
    sol: True en caras que miran al cielo (z+). Las verticales no reciben sol aqui.
    resistencias: {vertida_id: (R [m2K/W], retiro_h o None)}; retiro_h se mide desde
        el inicio de ESA vertida. None = no se retira en el horizonte.
    h_fijo: si se da, reemplaza h_conv(V) + h_rad (sirve para reproducir E07).
    excluir: vertida_id que no intercambian por esta cara (p. ej. costados del suelo).
    """
    clima: object
    sol: bool = False
    absortividad: float = 0.65
    emisividad: float = 0.9
    resistencias: dict = field(default_factory=dict)
    h_fijo: float | None = None
    excluir: tuple = ()

    def __post_init__(self):
        if not 0 <= self.absortividad <= 1 or not 0 <= self.emisividad <= 1:
            raise ValueError('Absortividad y emisividad deben estar entre 0 y 1.')
        for vid, (R, retiro) in self.resistencias.items():
            if not np.isfinite(R) or R < 0 or (retiro is not None and retiro < 0):
                raise ValueError(f'Vertida {vid}: R >= 0 y retiro >= 0.')
        if self.h_fijo is not None and not self.h_fijo > 0:
            raise ValueError('h_fijo debe ser positivo.')

    def h(self, T_C, viento):
        if self.h_fijo is not None:
            return np.full(np.shape(T_C), float(self.h_fijo)) if np.ndim(T_C) else float(self.h_fijo)
        return h_conveccion(viento) + h_radiacion(T_C, self.emisividad)

    def R_en(self, t_s, ids, sistema):
        R = np.zeros(len(ids))
        for vid, (R0, retiro) in self.resistencias.items():
            if vid not in sistema.programa.vertidas:
                continue
            fin = np.inf if retiro is None else sistema.programa.vertidas[vid]['inicio_s'] + 3600.0 * retiro
            if t_s < fin:
                R[ids == vid] = R0
        return R

    def coeficientes(self, t_s, cara, mascara, sistema):
        e = CARAS[cara][0]
        ids = sistema.programa.grilla.vertida_id[mascara]
        T_aire, viento, ghi = self.clima.en(t_s)
        h = self.h(T_aire, viento)
        G = sistema.areas[e] / (sistema.d[e] / (2 * sistema.k[mascara]) + self.R_en(t_s, ids, sistema) + 1.0 / h)
        if self.excluir:
            G[np.isin(ids, self.excluir)] = 0.0
        T_ref = T_aire + (self.absortividad * ghi / h if self.sol else 0.0)
        return G, np.full(len(G), T_ref), np.zeros(len(G))

    def G_maximo(self, cara, mascara, sistema):
        e = CARAS[cara][0]
        h_max = float(np.max(self.h(np.asarray(self.clima.T_aire_C), np.asarray(self.clima.viento_m_s))))
        G = sistema.areas[e] / (sistema.d[e] / (2 * sistema.k[mascara]) + 1.0 / h_max)
        if self.excluir:
            G[np.isin(sistema.programa.grilla.vertida_id[mascara], self.excluir)] = 0.0
        return G


@dataclass
class BordeTemperatura:
    """Temperatura impuesta (Dirichlet) solo en las celdas de ciertas vertidas.

    Sirve para el fondo del bloque de suelo: a esa profundidad la temperatura casi
    no cambia en 14 dias. En las demas celdas de la cara el borde es adiabatico.
    """
    T_C: float
    solo: tuple

    def coeficientes(self, t_s, cara, mascara, sistema):
        G = self.G_maximo(cara, mascara, sistema)
        return G, np.full(len(G), float(self.T_C)), np.zeros(len(G))

    def G_maximo(self, cara, mascara, sistema):
        e = CARAS[cara][0]
        G = sistema.areas[e] / (sistema.d[e] / (2 * sistema.k[mascara]))
        G[~np.isin(sistema.programa.grilla.vertida_id[mascara], self.solo)] = 0.0
        return G


def caras_ambientales(clima, *, moldaje=None, aislacion=None, suelo=None, absortividad=0.65,
                      emisividad=0.9, h_fijo=None):
    """Arma las seis caras del caso con clima real.

    moldaje / aislacion: {vertida_id: (R, retiro_h)} en lados / arriba.
    suelo: dict(id=..., T_fondo_C=...) si hay bloque de suelo; sin suelo la base es adiabatica.
    """
    from mcplan.solver1d import Neumann
    excluir = (suelo['id'],) if suelo else ()
    lados = BordeAmbiental(clima, sol=False, absortividad=absortividad, emisividad=emisividad,
                           resistencias=dict(moldaje or {}), h_fijo=h_fijo, excluir=excluir)
    arriba = BordeAmbiental(clima, sol=True, absortividad=absortividad, emisividad=emisividad,
                            resistencias=dict(aislacion or {}), h_fijo=h_fijo)
    base = BordeTemperatura(suelo['T_fondo_C'], (suelo['id'],)) if suelo else Neumann(0.0)
    return {'x-': lados, 'x+': lados, 'y-': lados, 'y+': lados, 'z-': base, 'z+': arriba}
