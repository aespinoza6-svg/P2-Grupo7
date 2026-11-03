"""E07: estructura temporal de vertidas sobre la grilla del IFC.

Separa identidad geometrica (vertida_id), existencia temporal (activa) y
temperatura. Aun NO integra conduccion ni condiciones de borde irregulares.
Las celdas sin colocar y el exterior tienen temperatura NaN: no son bloques
frios ficticios. La temperatura inicial solo se asigna al activar una vez.
"""
import numpy as np
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import HidratacionExponencial


class ProgramaVertidas:
    def __init__(self,grilla:GrillaVoxeles,modelos:dict[str,HidratacionExponencial]):
        self.grilla=grilla
        self.modelos=modelos
        self.vertidas={};self.mascaras={}
        for vid,props in grilla.vertidas.items():
            try:
                inicio_s=3600*float(props['Inicio_h'])
                temperatura=float(props['TemperaturaColocacion_C'])
                mezcla=props['TipoHormigon']
            except (KeyError,TypeError,ValueError) as e:
                raise ValueError(f'Vertida {vid}: faltan Inicio_h, TemperaturaColocacion_C o TipoHormigon validos.') from e
            if not np.isfinite([inicio_s,temperatura]).all() or inicio_s<0:
                raise ValueError(f'Vertida {vid}: horario o temperatura de colocacion invalida.')
            if mezcla not in modelos:
                raise ValueError(f'No hay modelo de hidratacion para mezcla {mezcla}.')
            self.vertidas[vid]=dict(inicio_s=inicio_s,temperatura_colocacion_C=temperatura,mezcla=mezcla)
            self.mascaras[vid]=grilla.vertida_id==vid

    @property
    def eventos_s(self):
        return np.unique([v['inicio_s'] for v in self.vertidas.values()])

    @staticmethod
    def _tiempo(t):
        if not np.isfinite(t) or t<0:
            raise ValueError('Tiempo global: se esperan segundos finitos no negativos.')
        return float(t)

    def mascara_activa(self,t_s):
        t=self._tiempo(t_s)
        activa=np.zeros(self.grilla.vertida_id.shape,dtype=bool)
        for vid,v in self.vertidas.items():
            if t>=v['inicio_s']:activa|=self.mascaras[vid]
        return activa

    def edad_s(self,t_s):
        t=self._tiempo(t_s)
        edad=np.full(self.grilla.vertida_id.shape,np.nan)
        for vid,v in self.vertidas.items():
            if t>=v['inicio_s']:edad[self.mascaras[vid]]=t-v['inicio_s']
        return edad

    def fuente_W_m3(self,t_s):
        t=self._tiempo(t_s)
        q=np.zeros(self.grilla.vertida_id.shape)
        for vid,v in self.vertidas.items():
            if t>=v['inicio_s']:
                q[self.mascaras[vid]]=self.modelos[v['mezcla']].fuente_W_m3(t-v['inicio_s'])
        return q

    def energia_generada_entre_J_m3(self,t0_s,t1_s):
        t0=self._tiempo(t0_s);t1=self._tiempo(t1_s)
        if t1<t0:raise ValueError('El tiempo final no puede ser menor que el inicial.')
        energia=np.zeros(self.grilla.vertida_id.shape)
        for vid,v in self.vertidas.items():
            energia[self.mascaras[vid]]=self.modelos[v['mezcla']].energia_entre_J_m3(t0-v['inicio_s'],t1-v['inicio_s'])
        return energia

    def partir_intervalo(self,t0_s,t1_s):
        """Extremos para que ningun paso de conduccion cruce una colocacion."""
        t0=self._tiempo(t0_s);t1=self._tiempo(t1_s)
        if t1<=t0:raise ValueError('El intervalo debe tener duracion positiva.')
        return np.r_[t0,self.eventos_s[(self.eventos_s>t0)&(self.eventos_s<t1)],t1]


class EstadoConstruccion:
    def __init__(self,programa:ProgramaVertidas):
        self.programa=programa
        self.activa=np.zeros(programa.grilla.vertida_id.shape,dtype=bool)
        self.temperatura_C=np.full(programa.grilla.vertida_id.shape,np.nan)
        self.t_actual_s=None

    def activar_hasta(self,t_s):
        """Incorpora vertidas pendientes sin reiniciar temperaturas existentes.
        No simula la evolucion termica entre ambos tiempos.
        """
        t=self.programa._tiempo(t_s)
        if self.t_actual_s is not None and t<self.t_actual_s:
            raise ValueError('El estado de construccion no puede retroceder en el tiempo.')
        objetivo=self.programa.mascara_activa(t)
        nuevas=objetivo&~self.activa
        for vid,v in self.programa.vertidas.items():
            self.temperatura_C[nuevas&self.programa.mascaras[vid]]=v['temperatura_colocacion_C']
        self.activa=objetivo
        self.t_actual_s=t
        return nuevas
