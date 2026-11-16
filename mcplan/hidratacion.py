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


class HidratacionTabulada:
    """E13: fuente medida Q(a) [W/m3 de hormigon] por edad, tal como la entrega el ensayo.

    Entre puntos, Q es lineal; la energia de cada paso es la integral EXACTA de esa recta
    (trapecio sin error), asi el balance del solver sigue cerrando a precision de maquina.
    Antes de la colocacion (a < 0) y despues del ultimo dato, Q = 0: no se inventa cola.
    rho y cp son los del ensayo: la curva y sus propiedades viajan juntas.
    """
    def __init__(self,edad_h,fuente_W_m3,rho_kg_m3,cp_J_kgK,nombre=''):
        a=3600.*np.asarray(edad_h,dtype=float);q=np.asarray(fuente_W_m3,dtype=float)
        ok=np.isfinite(a)&np.isfinite(q);a,q=a[ok],q[ok]
        if a.size<2 or np.any(np.diff(a)<=0) or a[0]<0 or np.any(q<0):
            raise ValueError(f'Curva {nombre}: edades crecientes desde 0 y Q >= 0.')
        if not (rho_kg_m3>0 and cp_J_kgK>0):raise ValueError('rho y cp deben ser positivos.')
        self.nombre=nombre;self.a_s=a;self.q=q;self.rho_kg_m3=float(rho_kg_m3);self.cp_J_kgK=float(cp_J_kgK)
        # Energia acumulada en cada punto de la tabla [J/m3]
        self.E=np.r_[0.,np.cumsum(.5*(q[1:]+q[:-1])*np.diff(a))]

    @property
    def capacidad_volumetrica(self):
        return self.rho_kg_m3*self.cp_J_kgK

    @property
    def energia_total_J_m3(self):
        return float(self.E[-1])

    def fuente_W_m3(self,edad_s):
        a=np.asarray(edad_s,dtype=float)
        return np.interp(a,self.a_s,self.q,left=0.,right=0.)*((a>=self.a_s[0])&(a<=self.a_s[-1]))

    def _energia_hasta(self,edad_s):
        a=np.clip(np.asarray(edad_s,dtype=float),self.a_s[0],self.a_s[-1])
        i=np.clip(np.searchsorted(self.a_s,a,side='right')-1,0,len(self.a_s)-2)
        s=a-self.a_s[i];h=self.a_s[i+1]-self.a_s[i]
        pendiente=(self.q[i+1]-self.q[i])/h
        return self.E[i]+self.q[i]*s+.5*pendiente*s*s

    def energia_entre_J_m3(self,edad_inicial_s,edad_final_s):
        a,b=np.broadcast_arrays(np.asarray(edad_inicial_s,dtype=float),np.asarray(edad_final_s,dtype=float))
        if not np.isfinite(a).all() or not np.isfinite(b).all() or np.any(b<a):
            raise ValueError('El intervalo de edad debe ser finito y no decreciente.')
        return self._energia_hasta(b)-self._energia_hasta(a)

    def elevacion_adiabatica(self,edad_s):
        return self._energia_hasta(edad_s)/self.capacidad_volumetrica


def leer_curvas_tabuladas(ruta):
    """Lee el archivo del profesor (source_terms.txt): columnas time_h y Q [W/m3] por curva.

    rho y cp salen del encabezado del mismo archivo ('rho = ... kg/m^3, cp = ... J/(kg K)').
    Devuelve {nombre: HidratacionTabulada}. NaN = fin del registro de esa curva.
    """
    import re
    texto=open(ruta,encoding='utf-8').read()
    m=re.search(r'rho\s*=\s*([\d.]+)\s*kg/m\^3,\s*cp\s*=\s*([\d.]+)',texto)
    if not m:raise ValueError('No encuentro rho y cp en el encabezado.')
    rho,cp=float(m.group(1)),float(m.group(2))
    cabeza=[l for l in texto.splitlines() if l.startswith('# time_h')]
    if not cabeza:raise ValueError('No encuentro la fila de nombres "# time_h ...".')
    nombres=cabeza[0].lstrip('#').split()[1:]
    datos=np.genfromtxt(ruta,comments='#')
    if datos.shape[1]!=len(nombres)+1:raise ValueError('Numero de columnas distinto a los nombres.')
    return {n:HidratacionTabulada(datos[:,0],datos[:,j+1],rho,cp,n) for j,n in enumerate(nombres)}


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
