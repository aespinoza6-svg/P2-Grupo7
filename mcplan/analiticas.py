"""
Soluciones exactas de la conducción (1D, 2D y 3D), usadas para VERIFICAR el solver.

Verificar = comprobar que el programa resuelve bien la ecuación que dice
resolver (comparando con respuestas conocidas). Es distinto de VALIDAR,
que es comprobar que la ecuación representa bien la realidad (eso se hará
comparando con las mediciones de Ju y Lei, 2019).

Cada función indica el problema exacto que resuelve, para poder reproducirlo
con el solver numérico. La explicación está en docs/1_teoria.md (Partes C y D).

Notación del enunciado E05 (la misma que usamos aquí):
    T0  temperatura de las caras (y de referencia)        [°C]
    A   amplitud inicial del seno                         [°C]
    λ   lambda = (pi/Lx)² + (pi/Ly)² [+ (pi/Lz)²]          [1/m²]
    q0  fuente volumétrica uniforme (en el solver: fuente) [W/m³]
"""

from __future__ import annotations

import numpy as np


def lambda_seno(largos):
    """λ = (pi/Lx)² [+ (pi/Ly)² + (pi/Lz)²]  [1/m²], una suma por dirección.

    Dice qué tan rápido se apaga el seno: la amplitud decae como exp(-alpha·λ·t).
    Mientras más chico el dominio o más direcciones, más grande λ y más rápido
    se enfría.
    """
    return float(sum((np.pi / L) ** 2 for L in largos))


def seno_dirichlet(x, t, L, alpha, T0, A):
    """Decaimiento de un perfil sinusoidal con ambos bordes a T0 (caso 1D).

    Problema:  dT/dt = alpha·d²T/dx²,  T(0,t) = T(L,t) = T0,
               T(x,0) = T0 + A·sin(pi·x/L)

    Solución:  T(x,t) = T0 + A·sin(pi·x/L)·exp(-alpha·λ·t),   λ = (pi/L)²

    (Por separación de variables: el modo sin(pi·x/L) cumple los bordes y
    decae con su propia tasa alpha·λ.)
    """
    x = np.asarray(x, dtype=float)
    return T0 + A * np.sin(np.pi * x / L) * np.exp(-alpha * lambda_seno([L]) * t)


def seno_dirichlet_nd(coordenadas, t, largos, alpha, T0, A):
    """Casos 2D y 3D del enunciado E05 (todas las caras a T0).

    2D:  T(x,y,t)   = T0 + A·sin(pi·x/Lx)·sin(pi·y/Ly)·exp(-alpha·λ·t)
    3D:  T(x,y,z,t) = T0 + A·sin(pi·x/Lx)·sin(pi·y/Ly)·sin(pi·z/Lz)·exp(-alpha·λ·t)

    con λ = (pi/Lx)² + (pi/Ly)² [+ (pi/Lz)²]  (ver lambda_seno).

    coordenadas: (X, Y) o (X, Y, Z) en forma de malla (np.meshgrid con
    indexing="ij"); largos: (Lx, Ly) o (Lx, Ly, Lz). Cada dirección suma su
    término a λ, por eso en 3D el seno se apaga más rápido.
    """
    modo = np.ones(np.shape(coordenadas[0]))
    for X, L in zip(coordenadas, largos):
        modo = modo * np.sin(np.pi * np.asarray(X, dtype=float) / L)
    return T0 + A * modo * np.exp(-alpha * lambda_seno(largos) * t)


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


def adiabatico_fuente_uniforme(t, T0, q0, rho, cp):
    """Bloque aislado (sin pérdidas) con fuente uniforme q0 [W/m³].

    Todo el calor generado se acumula:  rho·cp·dT/dt = q0
        dT/dt = q0/(rho·cp)   ->   T(t) = T0 + q0·t/(rho·cp)   (aumento lineal)

    q0 es VOLUMÉTRICA (W/m³): en el solver es el argumento fuente. No confundir
    con el flujo superficial q de Neumann (W/m²).
    Es la versión más simple del "aumento adiabático" del hormigón.
    """
    return T0 + q0 * np.asarray(t, dtype=float) / (rho * cp)
