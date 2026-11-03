"""Primera simulacion acoplada: IFC -> voxeles -> vertidas -> calor -> difusion.

DEMOSTRACION: hidratacion pedagogica no calibrada, ambiente constante,
base adiabatica y contacto perfecto. No es el Plan de Control final.
Requiere la grilla de scripts/07_voxelizar_ifc.py.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mcplan import bim
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import modelos_desde_configuracion
from mcplan.construccion import ProgramaVertidas
from mcplan.solver1d import Material,Robin,Neumann
from mcplan.solver3d import CARAS
from mcplan.solver_vertidas import resolver_vertidas

RAIZ=Path(__file__).resolve().parents[1]


def csv_salida(ruta,filas):
    with ruta.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(filas[0]));w.writeheader();w.writerows(filas)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--grilla',type=Path,default=RAIZ/'resultados'/'e07_voxelizacion'/'grilla_vertidas.npz')
    p.add_argument('--termica',type=Path,default=RAIZ/'config'/'termica_e07.yaml')
    p.add_argument('--hidratacion',type=Path,default=RAIZ/'config'/'hidratacion_e07.yaml')
    p.add_argument('--horizonte-h',type=float,default=None)
    p.add_argument('--max-dt-s',type=float,default=None,help='Limita el paso temporal sin editar el YAML; sirve para comparar precision')
    p.add_argument('--salida',type=Path,default=RAIZ/'resultados'/'e07_simulacion_acoplada')
    args=p.parse_args()
    if not args.grilla.is_file():p.error('Ejecute primero scripts/07_voxelizar_ifc.py.')
    rastreo=args.grilla.with_name('resumen.json')
    if rastreo.is_file():
        datos=json.loads(rastreo.read_text(encoding='utf-8'));ruta=Path(datos['ifc'])
        if ruta.is_file() and hashlib.sha256(ruta.read_bytes()).hexdigest()!=datos['ifc_sha256']:
            p.error('El IFC cambio. Repita la voxelizacion antes de simular.')
    g=GrillaVoxeles.cargar(args.grilla)
    caso=bim.cargar_yaml(RAIZ/'config'/'caso_base.yaml')
    hc=bim.cargar_yaml(args.hidratacion);tc=bim.cargar_yaml(args.termica)
    modelos=modelos_desde_configuracion(caso,hc)
    programa=ProgramaVertidas(g,modelos)
    h=caso['hormigon']
    materiales={mezcla:Material(mezcla,h['conductividad_W_mK'],h['densidad_kg_m3'],h['calor_especifico_J_kgK']) for mezcla in modelos}
    if tc['base']['condicion']!='adiabatica':p.error('Esta demostracion solo configura base adiabatica; debe modelarse el suelo posteriormente.')
    aire=Robin(float(tc['ambiente']['h_W_m2K']),float(tc['ambiente']['temperatura_C']))
    bordes={c:aire for c in CARAS};bordes['z-']=Neumann(0.)
    horizonte=float(args.horizonte_h if args.horizonte_h is not None else tc['horizonte_h'])
    max_dt=float(args.max_dt_s if args.max_dt_s is not None else tc['max_dt_s'])
    intervalo=float(tc['intervalo_salida_h'])
    if not np.isfinite([horizonte,intervalo,max_dt]).all() or min(horizonte,intervalo,max_dt)<=0:p.error('Horizonte, intervalo de salida y paso temporal positivos y finitos.')
    tc={**tc,'horizonte_h':horizonte,'max_dt_s':max_dt}
    tiempos=np.unique(np.r_[np.arange(0,horizonte,intervalo),horizonte])*3600
    print('SIMULACION DEMOSTRATIVA: hidratacion no calibrada; aire constante; base aislada; sin suelo explicito.',flush=True)
    print(f'Grilla: {g.vertida_id.shape}; {np.count_nonzero(g.vertida_id)} celdas de hormigon; horizonte {horizonte:g} h.',flush=True)
    reloj=time.monotonic();ultimo=[reloj]
    def progreso(msg):
        ahora=time.monotonic()
        if msg.startswith('Colocacion') or ahora-ultimo[0]>12:
            print(msg,flush=True);ultimo[0]=ahora
    r=resolver_vertidas(programa,materiales,bordes,horizonte*3600,
                       max_dt_s=max_dt,tiempos_guardado_s=tiempos,progreso=progreso)
    escala=np.array([max(1.,abs(f['energia_colocacion_J'])+abs(f['energia_hidratacion_J'])+abs(f['energia_bordes_J'])) for f in r.balance])
    residuos=np.array([f['residuo_balance_J'] for f in r.balance])
    relativo=float(np.max(abs(residuos)/escala))
    assert relativo<1e-8,'No cierra el balance energetico global.'
    out=args.salida;out.mkdir(parents=True,exist_ok=True)
    csv_salida(out/'balance_energia.csv',r.balance);csv_salida(out/'historias_por_vertida.csv',r.por_vertida)
    np.savez_compressed(out/'resultado_termico.npz',tiempo_s=r.t,T_final_C=r.T_final,
                        activa_final=r.activa_final,vertida_id=g.vertida_id,
                        cortes_xz_C=r.cortes_xz,y_corte_m=np.array(r.y_corte_m),
                        origen_m=g.origen,paso_m=np.array(g.paso))
    resumen=dict(alcance='Simulacion numerica acoplada demostrativa; no es prediccion validada ni plan constructivo recomendado.',
                 hidratacion=hc,termica=tc,horizonte_h=horizonte,
                 grilla_sha256=hashlib.sha256(args.grilla.read_bytes()).hexdigest(),
                 maximo_global_C=r.maximo_C,tiempo_maximo_h=r.tiempo_maximo_s/3600,
                 pasos=r.n_pasos,residuo_balance_max_J=float(np.max(abs(residuos))),
                 residuo_balance_relativo_max=relativo,
                 energia_final_J=r.balance[-1]['energia_J'],
                 energia_colocacion_J=r.balance[-1]['energia_colocacion_J'],
                 energia_hidratacion_J=r.balance[-1]['energia_hidratacion_J'],
                 energia_bordes_J=r.balance[-1]['energia_bordes_J'],
                 diferencial_normativo='No evaluado: faltan definicion y muestreo de nucleo-superficie a 25 mm',
                 limitaciones=['Parametros termicos e hidratacion pendientes de calibracion','Aire constante ilustrativo, no clima de Renaico',
                               'Base aislada provisional, sin suelo','Superficie escalonada de voxeles',
                               'Contacto perfecto y sin retiro de moldajes','No se incorpora dependencia de hidratacion con temperatura'])
    (out/'resumen.json').write_text(json.dumps(resumen,ensure_ascii=False,indent=2),encoding='utf-8')
    horas=r.t/3600
    fig,axs=plt.subplots(2,1,figsize=(10,8),layout='constrained')
    colores=['#287bb5','#d97922','#258b62','#985bad','#cd4c44']
    for i,vid in enumerate(sorted(programa.vertidas)):
        filas=[f for f in r.por_vertida if f['vertida_id']==vid]
        maxima=[np.nan if f['T_max_C'] is None else f['T_max_C'] for f in filas]
        media=[np.nan if f['T_media_C'] is None else f['T_media_C'] for f in filas]
        axs[0].plot(horas,maxima,color=colores[i%5],label=f'Vertida {vid}: máxima')
        axs[0].plot(horas,media,color=colores[i%5],ls='--',alpha=.7)
    axs[0].set(title='Temperaturas calculadas: máxima (continua) y media (discontinua)',ylabel='T [°C]')
    axs[0].legend(fontsize=8,ncol=3)
    axs[1].step(horas,[f['volumen_activo_m3'] for f in r.balance],where='post',color='#287bb5')
    axs[1].set(title='Volumen incorporado según el IFC',ylabel='Volumen activo [m³]')
    for a in axs:a.set_xlabel('Tiempo desde primera vertida [h]');a.grid(alpha=.2)
    fig.suptitle('Fundación voxelizada con difusión y vertidas · supuestos pedagógicos')
    fig.savefig(out/'temperaturas_por_vertida.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(10,7),layout='constrained')
    energia=np.array([f['energia_J'] for f in r.balance])
    for clave,etiqueta in [('energia_colocacion_J','Incorporada con las vertidas'),('energia_hidratacion_J','Generada por hidratación'),('energia_bordes_J','Intercambio neto por bordes')]:
        axs[0].plot(horas,np.array([f[clave] for f in r.balance])/1e9,label=etiqueta)
    axs[0].plot(horas,energia/1e9,color='black',lw=2,label='Energía almacenada')
    axs[0].set(ylabel='Energía [GJ]',title='Balance global, energía interna referida a 0 °C');axs[0].legend(fontsize=8)
    axs[1].plot(horas,residuos,color='#985bad');axs[1].set(ylabel='Residuo del balance [J]',title=f'Error relativo máximo: {relativo:.2e}')
    for a in axs:a.set_xlabel('Tiempo [h]');a.grid(alpha=.2)
    fig.savefig(out/'balance_energia.png',dpi=150);plt.close(fig)
    solicitar=[0,6,12,24,48,horizonte]
    indices=list(dict.fromkeys(int(np.argmin(abs(horas-t))) for t in solicitar if t<=horizonte))
    columnas=3;filas=int(np.ceil(len(indices)/columnas))
    fig,axs=plt.subplots(filas,columnas,figsize=(12,2.4*filas),layout='constrained',squeeze=False)
    extension=(g.origen[0],g.origen[0]+g.paso[0]*g.vertida_id.shape[0],g.origen[2],g.origen[2]+g.paso[2]*g.vertida_id.shape[2])
    cmap=plt.get_cmap('inferno').copy();cmap.set_bad('#e5e7e9')
    for a,j in zip(axs.flat,indices):
        im=a.imshow(np.ma.masked_invalid(r.cortes_xz[j]),origin='lower',extent=extension,cmap=cmap,
                    vmin=float(tc['ambiente']['temperatura_C']),vmax=max(20.,r.maximo_C),interpolation='nearest')
        a.set(title=f't = {horas[j]:g} h',xlabel='x [m]',ylabel='z [m]')
    for a in list(axs.flat)[len(indices):]:a.set_visible(False)
    fig.colorbar(im,ax=list(axs.flat),label='Temperatura [°C]',shrink=.8)
    fig.suptitle(f'Cortes térmicos y = {r.y_corte_m:.3f} m · gris = celda no activa o exterior',fontsize=12)
    fig.savefig(out/'cortes_termicos.png',dpi=150);plt.close(fig)
    print(f'Tmax global={r.maximo_C:.3f} C, en t={r.tiempo_maximo_s/3600:.2f} h; {r.n_pasos} pasos.',flush=True)
    print(f'Residuo energetico relativo maximo: {relativo:.3e}. Tiempo de calculo y graficos: {time.monotonic()-reloj:.1f} s.',flush=True)
    print(f'Resultados: {out}',flush=True)


if __name__=='__main__':main()
