"""E07: activacion del IFC voxelizado y diseno pedagogico de hidratacion.

NO simula la difusion termica de la fundacion. Las curvas de temperatura
son referencias ADIABATICAS por vertida, sin intercambio entre ellas.
Ejecutar primero scripts/07_voxelizar_ifc.py y despues este archivo.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from mcplan import bim
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import modelos_desde_configuracion
from mcplan.construccion import ProgramaVertidas,EstadoConstruccion

RAIZ=Path(__file__).resolve().parents[1]
COLORES=['#287bb5','#d97922','#258b62','#985bad','#cd4c44']


def escribir_csv(ruta,filas):
    with ruta.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(filas[0]));w.writeheader();w.writerows(filas)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--grilla',type=Path,default=RAIZ/'resultados'/'e07_voxelizacion'/'grilla_vertidas.npz')
    p.add_argument('--config',type=Path,default=RAIZ/'config'/'hidratacion_e07.yaml')
    p.add_argument('--salida',type=Path,default=RAIZ/'resultados'/'e07_construccion')
    p.add_argument('--horizonte-h',type=float,default=72)
    args=p.parse_args()
    if not args.grilla.is_file():
        p.error('Primero ejecute scripts/07_voxelizar_ifc.py para generar la grilla.')
    if not np.isfinite(args.horizonte_h) or args.horizonte_h<=0:
        p.error('El horizonte debe ser positivo y finito.')
    resumen_grilla=args.grilla.with_name('resumen.json')
    if resumen_grilla.is_file():
        trazabilidad=json.loads(resumen_grilla.read_text(encoding='utf-8'))
        ifc=Path(trazabilidad['ifc'])
        if ifc.is_file() and hashlib.sha256(ifc.read_bytes()).hexdigest()!=trazabilidad['ifc_sha256']:
            p.error('El IFC cambio desde la voxelizacion. Repita scripts/07_voxelizar_ifc.py.')
    g=GrillaVoxeles.cargar(args.grilla)
    caso=bim.cargar_yaml(RAIZ/'config'/'caso_base.yaml')
    config=bim.cargar_yaml(args.config)
    modelos=modelos_desde_configuracion(caso,config)
    programa=ProgramaVertidas(g,modelos)
    out=args.salida;out.mkdir(parents=True,exist_ok=True)
    print('E07: DISENO PEDAGOGICO NO CALIBRADO. No es una prediccion termica de la fundacion.',flush=True)
    print(f'Grilla leida: {args.grilla}',flush=True)
    filas=[]
    for vid,v in sorted(programa.vertidas.items()):
        m=modelos[v['mezcla']]
        filas.append(dict(vertida_id=vid,inicio_h=v['inicio_s']/3600,
                          temperatura_colocacion_C=v['temperatura_colocacion_C'],mezcla=v['mezcla'],
                          volumen_voxel_m3=float(np.count_nonzero(programa.mascaras[vid])*g.volumen_celda),
                          elevacion_adiabatica_final_C=m.elevacion_final_C,tau_h=m.tau_s/3600,
                          fuente_inicial_W_m3=float(m.fuente_W_m3(0))))
        print(f"Vertida {vid}: inicio {v['inicio_s']/3600:g} h; colocacion {v['temperatura_colocacion_C']:g} C; q inicial {m.fuente_W_m3(0):.1f} W/m3",flush=True)
    escribir_csv(out/'programa_vertidas.csv',filas)

    # Activacion de celdas: no evoluciona temperaturas entre eventos.
    estado=EstadoConstruccion(programa);instantaneas=[];eventos=[]
    for t in programa.eventos_s:
        nuevas=estado.activar_hasta(float(t))
        activos=np.unique(g.vertida_id[estado.activa]).tolist()
        eventos.append(dict(t_h=float(t/3600),vertidas_activas=activos,
                            celdas_nuevas=int(nuevas.sum()),celdas_activas=int(estado.activa.sum()),
                            volumen_activo_m3=float(estado.activa.sum()*g.volumen_celda)))
        print(f'Evento {t/3600:g} h: activas {activos}; volumen colocado {eventos[-1]["volumen_activo_m3"]:.3f} m3',flush=True)
        instantaneas.append((float(t),estado.activa.copy()))
    np.savez_compressed(out/'estructura_estado_ultima_colocacion.npz',vertida_id=g.vertida_id,
                        activa=estado.activa,temperatura_colocacion_asignada_C=estado.temperatura_C,
                        edad_s=programa.edad_s(estado.t_actual_s),
                        fuente_W_m3=programa.fuente_W_m3(estado.t_actual_s),
                        t_global_s=np.array(estado.t_actual_s),
                        alcance=np.array('Estado de incorporacion: no incluye evolucion termica por conduccion'))

    # Historias analiticas ADIABATICAS: cada vertida lleva su propio reloj.
    tiempos=np.unique(np.r_[np.linspace(0,args.horizonte_h*3600,289),programa.eventos_s])
    tiempos=tiempos[tiempos<=args.horizonte_h*3600]
    historias=[];curvas={};total_energia=np.zeros_like(tiempos)
    for fila in filas:
        vid=fila['vertida_id'];v=programa.vertidas[vid];m=modelos[v['mezcla']]
        edad=tiempos-v['inicio_s'];activa=edad>=0
        q=m.fuente_W_m3(edad)
        elevacion=m.elevacion_adiabatica(edad)
        temperatura=np.where(activa,v['temperatura_colocacion_C']+elevacion,np.nan)
        energia=m.capacidad_volumetrica*elevacion*fila['volumen_voxel_m3']
        total_energia+=energia;curvas[vid]=(q,temperatura)
        for t,a,existe,fuente,temp,e in zip(tiempos,edad,activa,q,temperatura,energia):
            historias.append(dict(t_global_h=float(t/3600),vertida_id=vid,activa=bool(existe),
                                  edad_h=float(a/3600) if existe else '',fuente_W_m3=float(fuente),
                                  T_adiabatica_referencia_C=float(temp) if existe else '',
                                  energia_hidratacion_acumulada_J=float(e)))
    escribir_csv(out/'historias_adiabaticas_por_vertida.csv',historias)
    fig,axs=plt.subplots(3,1,figsize=(11,10),layout='constrained')
    for i,(vid,(q,temp)) in enumerate(curvas.items()):
        c=COLORES[i%len(COLORES)]
        axs[0].plot(tiempos/3600,q,label=f'Vertida {vid}',color=c)
        axs[1].plot(tiempos/3600,temp,label=f'Vertida {vid}',color=c)
    volumenes=np.zeros_like(tiempos)
    for fila in filas:volumenes+=np.where(tiempos>=fila['inicio_h']*3600,fila['volumen_voxel_m3'],0)
    axs[2].step(tiempos/3600,volumenes,where='post',color='#287bb5')
    axs[0].set(title='Fuente de hidratación desplazada al inicio de cada vertida',ylabel='qₕ [W/m³]')
    axs[1].set(title='Referencia adiabática por vertida: sin pérdidas NI intercambio entre vertidas',ylabel='T colocación + θad [°C]')
    axs[2].set(title='Volumen voxelizado incorporado a la construcción',ylabel='Volumen activo [m³]')
    for a in axs:a.set_xlabel('Tiempo desde la primera vertida [h]');a.grid(alpha=.2)
    axs[0].legend(ncol=5,fontsize=9)
    fig.suptitle('E07 · modelo pedagógico de hidratación · parámetros pendientes de calibrar')
    fig.savefig(out/'hidratacion_y_activacion.png',dpi=150);plt.close(fig)

    # Cortes de la mascara temporal, conservando etiquetas geometricas.
    ids=sorted(programa.vertidas);cmap=ListedColormap(['white','#e0e2e5']+[COLORES[i%len(COLORES)] for i in range(len(ids))])
    norm=BoundaryNorm(np.arange(len(ids)+3)-.5,len(ids)+2)
    iy=int(np.argmin(abs(g.centros[1])));extension=(g.origen[0],g.origen[0]+g.paso[0]*g.vertida_id.shape[0],g.origen[2],g.origen[2]+g.paso[2]*g.vertida_id.shape[2])
    columnas=min(3,len(instantaneas));filas_plot=int(np.ceil(len(instantaneas)/columnas))
    fig,axs=plt.subplots(filas_plot,columnas,figsize=(12,2.3*filas_plot),layout='constrained',squeeze=False)
    for a,(t,activa) in zip(axs.flat,instantaneas):
        muestra=np.where(g.vertida_id[:,iy,:]>0,1,0)
        for i,vid in enumerate(ids):muestra[(g.vertida_id[:,iy,:]==vid)&activa[:,iy,:]]=i+2
        a.imshow(muestra.T,origin='lower',extent=extension,cmap=cmap,norm=norm,interpolation='nearest')
        a.set(title=f't = {t/3600:g} h',xlabel='x [m]',ylabel='z [m]')
    for a in list(axs.flat)[len(instantaneas):]:a.set_visible(False)
    fig.legend(handles=[Patch(color='#e0e2e5',label='Pendiente')]+[Patch(color=COLORES[i%len(COLORES)],label=f'Vertida {vid}') for i,vid in enumerate(ids)],loc='outside lower center',ncol=min(6,len(ids)+1),fontsize=9)
    fig.suptitle('Activación temporal del IFC voxelizado · color = colocado; gris = pendiente\nNo es un campo de temperatura',fontsize=12)
    fig.savefig(out/'activacion_vertidas_cortes.png',dpi=150);plt.close(fig)

    resumen=dict(alcance='Diseno de hidratacion y estructura de activacion. Sin simulacion de conduccion de la fundacion.',
                 estado_parametros=config['estado'],fuente_parametros=config['fuente'],
                 modelo='theta_ad(a)=theta_inf*(1-exp(-a/tau)); q_h=rho*cp*dtheta_ad/da',
                 edad='a=t_global-t_colocacion, independiente para cada vertida',
                 dependencia_de_temperatura='No implementada; futura extension de edad equivalente/Arrhenius a calibrar',
                 grilla_sha256=hashlib.sha256(args.grilla.read_bytes()).hexdigest(),
                 configuracion_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                 programa_vertidas=filas,eventos=eventos,
                 energia_total_generada_al_final_J=float(total_energia[-1]))
    (out/'resumen.json').write_text(json.dumps(resumen,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'Resultados: {out}',flush=True)


if __name__=='__main__':main()
