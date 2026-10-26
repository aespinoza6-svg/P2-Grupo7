"""mcplan — Planificación de vertidas de hormigón masivo mediante simulación.

Proyecto 2, Métodos Computacionales 2026, Grupo 7 (Renaico / P.E. Vergara).

Módulos disponibles (semanas 1 y 2):

    mcplan.solver1d      solver térmico 1D explícito (balance por celdas)
    mcplan.solver3d      el mismo solver en 2D/3D (malla cartesiana de voxeles)
    mcplan.analiticas    soluciones exactas usadas para verificar el solver
    mcplan.bim           creación y lectura del modelo IFC (IfcOpenShell)

Regla de arquitectura: ningún módulo de mcplan importa Blender ni Bonsai.
El único puente entre el modelo BIM y el cálculo es el archivo IFC.
"""

__version__ = "0.1.0"
