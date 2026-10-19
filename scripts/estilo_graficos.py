"""Estilo común de las figuras (lo importan los scripts 03 y 04)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
CARPETA_FIGURAS = RAIZ / "resultados" / "figuras"

# Colores categóricos (curvas que representan COSAS distintas).
AZUL = "#2a78d6"
NARANJO = "#eb6834"
VERDE_AGUA = "#1baf7a"
GRIS = "#52514e"   # límites, referencias y texto secundario

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 160,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 2.0,
    "font.size": 10,
})


def colores_tiempo(n):
    """Rampa de un solo tono (claro = temprano, oscuro = tarde) para perfiles
    en distintos instantes: el tiempo es una magnitud ordenada."""
    return plt.cm.Blues(np.linspace(0.35, 0.95, n))


def guardar(fig, nombre):
    CARPETA_FIGURAS.mkdir(parents=True, exist_ok=True)
    ruta = CARPETA_FIGURAS / nombre
    fig.savefig(ruta, bbox_inches="tight")
    print(f"  figura -> {ruta.relative_to(RAIZ)}")
    return ruta
