"""Diseno inicial E07: curva adiabatica -> fuente de hidratacion por EDAD.

Modelo pedagogico de la guia del curso, secciones 23-25:
theta(a) = theta_inf (1-exp(-a/tau)), a >= 0.
q_h(a) = rho cp theta_inf/tau exp(-a/tau) [W/m3].
La edad a se mide desde la colocacion de CADA vertida, no desde t=0 global.
No hay calibracion experimental ni dependencia de temperatura implementada.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class HidratacionExponencial:
    elevacion_final_C: float
    tau_s: float
    rho_kg_m3: float
    cp_J_kgK: float

    def __post_init__(self):
        valores=[self.elevacion_final_C,self.tau_s,self.rho_kg_m3,self.cp_J_kgK]
        if not np.isfinite(valores).all() or self.elevacion_final_C<0 or min(valores[1:])<=0:
            raise ValueError('Parametros de hidratacion invalidos: elevacion >=0; tau, rho y cp >0; todos finitos.')

    @property
    def capacidad_volumetrica(self):
        return self.rho_kg_m3*self.cp_J_kgK # J/(m3 K)

    def elevacion_adiabatica(self,edad_s):
        edad=np.asarray(edad_s,dtype=float)
        if not np.isfinite(edad).all():
            raise ValueError('La edad debe ser finita.')
        # Antes de colocacion, elevacion cero. expm1 evita perdida de precision.
        a=np.maximum(edad,0)
        return self.elevacion_final_C*(-np.expm1(-a/self.tau_s))

    def fuente_W_m3(self,edad_s):
        edad=np.asarray(edad_s,dtype=float)
        if not np.isfinite(edad).all():
            raise ValueError('La edad debe ser finita.')
        a=np.maximum(edad,0)
        q=self.capacidad_volumetrica*self.elevacion_final_C/self.tau_s*np.exp(-a/self.tau_s)
        return np.where(edad>=0,q,0.)

    def energia_entre_J_m3(self,edad_inicial_s,edad_final_s):
        a,b=np.broadcast_arrays(np.asarray(edad_inicial_s,dtype=float),np.asarray(edad_final_s,dtype=float))
        if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(b<a):
            raise ValueError('El intervalo de edad debe ser finito y no decreciente.')
        a=np.maximum(a,0);b=np.maximum(b,0)
        # Forma estable de rho cp [theta(b)-theta(a)], tambien antes de colocar.
        return self.capacidad_volumetrica*self.elevacion_final_C*np.exp(-a/self.tau_s)*(-np.expm1(-(b-a)/self.tau_s))


def modelos_desde_configuracion(caso:dict,configuracion:dict):
    """Parametros de curva por mezcla; rho y cp del caso base, comunes en E07."""
    if configuracion.get('referencia_edad')!='desde_colocacion':
        raise ValueError('La referencia de edad debe ser desde_colocacion.')
    h=caso['hormigon'];modelos={}
    for nombre,param in configuracion['mezclas'].items():
        if param.get('modelo')!='adiabatico_exponencial':
            raise ValueError(f'Modelo no implementado para mezcla {nombre}.')
        modelos[nombre]=HidratacionExponencial(
            elevacion_final_C=float(param['elevacion_final_C']),tau_s=3600*float(param['tau_h']),
            rho_kg_m3=float(h['densidad_kg_m3']),cp_J_kgK=float(h['calor_especifico_J_kgK']))
    return modelos
