"""
Soluciones exactas de la conducción 1D, usadas para VERIFICAR el solver.

Verificar = comprobar que el programa resuelve bien la ecuación que dice
resolver (comparando con respuestas conocidas). Es distinto de VALIDAR,
que es comprobar que la ecuación representa bien la realidad (eso se hará
comparando con las mediciones de Ju y Lei, 2019).

Cada función indica el problema exacto que resuelve, para poder reproducirlo
con el solver numérico. La deducción de cada una está en docs/03_solver_1d.md.
"""

from __future__ import annotations

import numpy as np


def seno_dirichlet(x, t, L, alpha, T_borde, amplitud):
    """Decaimiento de un perfil sinusoidal con ambos bordes a T_borde.

    Problema:  dT/dt = alpha·d²T/dx²,  T(0,t) = T(L,t) = T_borde,
               T(x,0) = T_borde + amplitud·sin(pi·x/L)

    Solución:  T(x,t) = T_borde + amplitud·sin(pi·x/L)·exp(-alpha·(pi/L)²·t)

    (Por separación de variables: el modo sin(pi·x/L) cumple los bordes y
    decae con su propia tasa alpha·(pi/L)².)
    """
    x = np.asarray(x, dtype=float)
    return T_borde + amplitud * np.sin(np.pi * x / L) * np.exp(-alpha * (np.pi / L) ** 2 * t)


def estacionario_dirichlet_robin(x, L, k, T_izq, h, T_inf, R_extra=0.0):
    """Estado estacionario: x = 0 a T_izq, x = L con convección (h, T_inf).

    Sin fuente y en régimen permanente el flujo q es el mismo en todo x
    (perfil lineal). Es un circuito de resistencias en serie:

        q = (T_izq - T_inf) / (L/k + R_extra + 1/h)       [W/m²]
        T(x) = T_izq - q·x/k
    """
    x = np.asarray(x, dtype=float)
    q = (T_izq - T_inf) / (L / k + R_extra + 1.0 / h)
    return T_izq - q * x / k


def estacionario_dos_capas(x, L1, k1, L2, k2, T_izq, T_der):
    """Estado estacionario en dos capas en serie con bordes Dirichlet.

        q = (T_izq - T_der) / (L1/k1 + L2/k2)
        capa 1 (x <= L1):  T = T_izq - q·x/k1
        capa 2 (x >  L1):  T = T_interfaz - q·(x - L1)/k2
    """
    x = np.asarray(x, dtype=float)
    q = (T_izq - T_der) / (L1 / k1 + L2 / k2)
    T_interfaz = T_izq - q * L1 / k1
    return np.where(x <= L1, T_izq - q * x / k1, T_interfaz - q * (x - L1) / k2)


def adiabatico_fuente_uniforme(t, T0, Q, rho, cp):
    """Bloque aislado (sin pérdidas) con fuente uniforme Q [W/m³].

    Todo el calor generado se acumula:  rho·cp·dT/dt = Q
        T(t) = T0 + Q·t/(rho·cp)        (aumento lineal)

    Es la versión más simple del "aumento adiabático" del hormigón.
    """
    return T0 + Q * np.asarray(t, dtype=float) / (rho * cp)
