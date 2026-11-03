"""Difusion explicita sobre el IFC voxelizado, con colocacion por eventos.

Solo las celdas activas tienen capacidad y participan en intercambios.
Una cara activa-inactiva es EXTERIOR; al colocar su vecina pasa a interfaz.
Los bordes se asignan por orientacion de la cara voxelizada, no por la caja.
Contacto entre hormigones perfecto, sin resistencia de junta adicional.
Hidratacion prescrita por edad; energia integrada exactamente por intervalo.
"""
from dataclasses import dataclass
from typing import Callable
import numpy as np
from mcplan.construccion import ProgramaVertidas,EstadoConstruccion
from mcplan.solver1d import Material,Dirichlet,Neumann,Robin,_conductancia_borde,_evaluar
from mcplan.solver3d import CARAS


def _lados(eje,n):
    izq=[slice(None)]*3;der=[slice(None)]*3
    izq[eje]=slice(0,n-1);der[eje]=slice(1,n)
    return tuple(izq),tuple(der)


class SistemaActivo:
    def __init__(self,programa:ProgramaVertidas,materiales:dict[str,Material],bordes=None):
        self.programa=programa;g=programa.grilla;self.forma=g.vertida_id.shape
        self.d=g.paso;self.V=g.volumen_celda
        # Valores auxiliares finitos fuera del hormigon. Nunca se activan alli.
        self.k=np.ones(self.forma);self.C=np.ones(self.forma)
        for vid,v in programa.vertidas.items():
            if v['mezcla'] not in materiales:raise ValueError(f'Falta material {v["mezcla"]}.')
            m=materiales[v['mezcla']]
            if not np.isfinite([m.k,m.rho,m.cp]).all() or min(m.k,m.rho,m.cp)<=0:
                raise ValueError('k, rho y cp deben ser finitos y positivos.')
            h=programa.modelos[v['mezcla']]
            if not np.isclose(m.rho,h.rho_kg_m3) or not np.isclose(m.cp,h.cp_J_kgK):
                raise ValueError('La capacidad termica del material y del modelo de hidratacion debe coincidir.')
            mask=programa.mascaras[vid]
            self.k[mask]=m.k;self.C[mask]=m.rho*m.cp*self.V
        bordes=dict(bordes or {})
        if set(bordes)-set(CARAS):raise ValueError('Nombre de cara desconocido.')
        self.bordes={c:bordes.get(c,Neumann(0)) for c in CARAS}
        self.areas=[self.V/d for d in self.d]
        self.Gbase=[];self.Gb={}
        for e,n in enumerate(self.forma):
            a,b=_lados(e,n)
            self.Gbase.append(self.areas[e]/(self.d[e]/(2*self.k[a])+self.d[e]/(2*self.k[b])))
        for cara,(e,_) in CARAS.items():
            borde=self.bordes[cara]
            if isinstance(borde,Robin) and (not np.isfinite([borde.h,borde.R_extra]).all() or borde.h<=0 or borde.R_extra<0):
                raise ValueError('Robin requiere h>0 y resistencia adicional >=0.')
            self.Gb[cara]=np.broadcast_to(_conductancia_borde(borde,self.k,self.d[e])*self.areas[e],self.forma)

    def actualizar(self,activa):
        self.activa=activa.copy();self.G=[];self.expuestas={};S=np.zeros(self.forma)
        for e,n in enumerate(self.forma):
            a,b=_lados(e,n)
            G=np.where(activa[a]&activa[b],self.Gbase[e],0.)
            self.G.append(G);S[a]+=G;S[b]+=G
        for cara,(e,lado) in CARAS.items():
            vecina=np.zeros(self.forma,dtype=bool)
            a,b=_lados(e,self.forma[e])
            if lado==0:vecina[b]=activa[a]
            else:vecina[a]=activa[b]
            mask=activa&~vecina
            self.expuestas[cara]=mask
            S[mask]+=self.Gb[cara][mask]
        seleccion=activa&(S>0)
        self.dt_critico=float(np.min(self.C[seleccion]/S[seleccion])) if seleccion.any() else np.inf
        return self.dt_critico

    def potencias_W(self,T,t):
        """Potencia neta por celda y total exterior entrante [W]."""
        segura=np.nan_to_num(T,nan=0.)
        balance=np.zeros(self.forma);exterior=0.
        for e,n in enumerate(self.forma):
            a,b=_lados(e,n)
            F=self.G[e]*(segura[b]-segura[a])
            balance[a]+=F;balance[b]-=F
        for cara,(e,_) in CARAS.items():
            mask=self.expuestas[cara];borde=self.bordes[cara]
            if isinstance(borde,Neumann):
                q=_evaluar(borde.q,t)*self.areas[e]
                if q:
                    balance[mask]+=q;exterior+=q*np.count_nonzero(mask)
            else:
                referencia=_evaluar(borde.T if isinstance(borde,Dirichlet) else borde.T_inf,t)
                q=self.Gb[cara][mask]*(referencia-T[mask])
                balance[mask]+=q;exterior+=float(q.sum())
        return balance,exterior


@dataclass
class ResultadoVertidas:
    t:np.ndarray
    T_final:np.ndarray
    activa_final:np.ndarray
    cortes_xz:np.ndarray
    y_corte_m:float
    balance:list[dict]
    por_vertida:list[dict]
    maximo_C:float
    tiempo_maximo_s:float
    n_pasos:int
    T:np.ndarray|None=None


def resolver_vertidas(programa,materiales,bordes,t_final_s,*,dt_s=None,max_dt_s=600.,
                      tiempos_guardado_s=None,guardar_campos=False,progreso:Callable|None=None):
    """Euler explicito con eventos, fuentes por edad y balance de energia.

    dt_s explicito se verifica; None usa min(max_dt_s,0.8*C/S).
    Siempre se divide el paso en colocaciones, tiempos de salida y final.
    Las salidas en un evento incluyen el nuevo hormigon ya incorporado.
    Ledger: E = E_colocacion + E_hidratacion + E_bordes, referido a 0 C.
    No se calcula aun un diferencial nucleo-superficie normativo.
    """
    if not np.isfinite(t_final_s) or t_final_s<=0:raise ValueError('Tiempo final positivo y finito.')
    if not np.isfinite(max_dt_s) or max_dt_s<=0:raise ValueError('max_dt_s debe ser positivo y finito.')
    if dt_s is not None and (not np.isfinite(dt_s) or dt_s<=0):raise ValueError('dt_s debe ser positivo y finito.')
    if tiempos_guardado_s is None:tiempos_guardado_s=np.linspace(0,t_final_s,97)
    salidas=np.asarray(tiempos_guardado_s,dtype=float)
    if not np.isfinite(salidas).all() or np.any(salidas<0) or np.any(salidas>t_final_s):
        raise ValueError('Tiempos de salida fuera del horizonte.')
    eventos=programa.eventos_s[programa.eventos_s<=t_final_s]
    paradas=np.unique(np.r_[0.,salidas,eventos,t_final_s])
    if guardar_campos and len(paradas)*programa.grilla.vertida_id.size*8>512*1024**2:
        raise ValueError('Historial 3D mayor a 512 MB. Use cortes o menos salidas.')
    sistema=SistemaActivo(programa,materiales,bordes);estado=EstadoConstruccion(programa)
    g=programa.grilla;centro_y=g.origen[1]+g.paso[1]*g.vertida_id.shape[1]/2
    iy=int(np.argmin(abs(g.centros[1]-centro_y)))
    Ecol=0.;Ehid=0.;Ebor=0.;maximo=-np.inf;tmax=0.;n_pasos=0
    historial=[];historias=[];cortes=[];campos=[];ts=[]
    def incorporar(t):
        nonlocal Ecol
        nuevas=estado.activar_hasta(t)
        if nuevas.any():
            Ecol+=float(np.sum(sistema.C[nuevas]*estado.temperatura_C[nuevas]))
        if nuevas.any() or not hasattr(sistema,'dt_critico'):
            sistema.actualizar(estado.activa)
        if progreso and nuevas.any():progreso(f'Colocacion en t={t/3600:g} h: {nuevas.sum()} celdas nuevas')
    def medir_maximo(t):
        nonlocal maximo,tmax
        if estado.activa.any():
            valor=float(estado.temperatura_C[estado.activa].max())
            if valor>maximo:maximo=valor;tmax=t
    def guardar(t):
        activa=estado.activa;T=estado.temperatura_C
        energia=float(np.sum(sistema.C[activa]*T[activa]))
        residuo=energia-Ecol-Ehid-Ebor
        historial.append(dict(tiempo_s=float(t),tiempo_h=float(t/3600),
                              celdas_activas=int(activa.sum()),volumen_activo_m3=float(activa.sum()*sistema.V),
                              T_min_C=float(T[activa].min()) if activa.any() else None,
                              T_max_C=float(T[activa].max()) if activa.any() else None,
                              energia_J=energia,energia_colocacion_J=Ecol,energia_hidratacion_J=Ehid,
                              energia_bordes_J=Ebor,residuo_balance_J=residuo,
                              dt_critico_s=sistema.dt_critico if np.isfinite(sistema.dt_critico) else None))
        for vid,v in programa.vertidas.items():
            mask=programa.mascaras[vid]&activa
            historias.append(dict(tiempo_h=float(t/3600),vertida_id=vid,activa=bool(mask.any()),
                                  T_min_C=float(T[mask].min()) if mask.any() else None,
                                  T_max_C=float(T[mask].max()) if mask.any() else None,
                                  T_media_C=float(T[mask].mean()) if mask.any() else None))
        ts.append(t);cortes.append(T[:,iy,:].T.astype(np.float32))
        if guardar_campos:campos.append(T.copy())
    incorporar(0.);medir_maximo(0.);guardar(0.)
    t=0.
    for parada in paradas[1:]:
        while t<parada-1e-8:
            preferido=dt_s if dt_s is not None else min(max_dt_s,.8*sistema.dt_critico)
            dt=min(preferido,parada-t)
            if dt>sistema.dt_critico*(1+1e-12):
                raise ValueError(f'dt={dt:g} s supera el criterio estable {sistema.dt_critico:g} s para el hormigon activo.')
            if dt<=0 or not np.isfinite(dt):raise ValueError('Paso temporal invalido.')
            potencia,pot_ext=sistema.potencias_W(estado.temperatura_C,t)
            energia_h=programa.energia_generada_entre_J_m3(t,t+dt)*sistema.V
            activa=estado.activa
            estado.temperatura_C[activa]+=(dt*potencia[activa]+energia_h[activa])/sistema.C[activa]
            Ehid+=float(energia_h[activa].sum());Ebor+=dt*pot_ext
            t+=dt;n_pasos+=1;medir_maximo(t)
        # Evita que error de redondeo retrase una activacion exacta.
        t=float(parada);incorporar(t);medir_maximo(t);guardar(t)
        if progreso:progreso(f't={t/3600:.1f} h / {t_final_s/3600:.1f} h; Tmax acumulada={maximo:.2f} C')
    return ResultadoVertidas(np.array(ts),estado.temperatura_C.copy(),estado.activa.copy(),
                            np.array(cortes),float(g.centros[1][iy]),historial,historias,
                            float(maximo) if np.isfinite(maximo) else np.nan,float(tmax),n_pasos,
                            np.array(campos) if guardar_campos else None)
