"""Ejercicio independiente del jueves: bloque 3 x 3 x 0,60 m.

Ejecutar este archivo con el Python .venv del repositorio P2-Grupo7.
Usa el solver propio mcplan.solver3d y no modifica el IFC.
La fuente dura 2 dias, con maximo de 1300 W/m3 en el dia 1.

BETA es el coeficiente de dT/dn = BETA*(T-20), normal EXTERIOR [1/m].
Para conveccion pasiva BETA < 0 y h=-k*BETA [W/(m2 K)].
Si BETA=None se ejecutan los primeros tres casos y se indica el pendiente.
"""
from pathlib import Path
from datetime import datetime
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from mcplan.solver3d import Material, Malla3D, Dirichlet, Neumann, Robin, CARAS, resolver, paso_tiempo_maximo

# ---------------- DATOS QUE SE PUEDEN MODIFICAR EN VIVO ----------------
LARGOS = (3.0, 3.0, 0.60)  # x: ancho; y: profundidad; z: altura [m]
CELDAS = (30, 30, 6)       # dx=dy=dz=0,10 m; 5400 voxeles
K, RHO, CP = 2.0, 2400.0, 900.0
T_INICIAL = 20.0
QMAX = 1300.0             # W/m3, fuente uniforme en todo el volumen
DIA = 86400.0
DURACION_PULSO = 2*DIA
T_FINAL = 2*DIA
DT = 120.0                # 2 minutos; menor al criterio de los cuatro casos
BETA = -5.0                # Primer ensayo: coincide con adiabatico. Normal exterior [1/m]
# --------------------------------------------------------------------


def fraccion_fuente(t):
    t = np.asarray(t)
    return np.where((t >= 0) & (t < DURACION_PULSO),
                    .5*(1-np.cos(2*np.pi*t/DURACION_PULSO)), 0.)


def integral_fuente(t):
    """Energia generada por m3 entre 0 y t [J/m3]."""
    u = np.clip(np.asarray(t,dtype=float),0,DURACION_PULSO)
    return QMAX/2*(u-DURACION_PULSO/(2*np.pi)*np.sin(2*np.pi*u/DURACION_PULSO))


def temperatura_centro(resultado):
    """Interpolacion al centro geometrico (1,5; 1,5; 0,30) m.
    Promedio de las dos celdas centrales por eje cuando el numero es par.
    """
    a=resultado.T
    for eje,n in reversed(list(enumerate(resultado.malla.forma))):
        indices=[n//2] if n%2 else [n//2-1,n//2]
        a=np.take(a,indices,axis=eje+1).mean(axis=eje+1)
    return a


def configurar_casos(beta):
    casos=[('dirichlet_20','Dirichlet: todas las caras a 20 °C',Dirichlet(20)),
           ('dirichlet_30','Dirichlet: todas las caras a 30 °C',Dirichlet(30)),
           ('adiabatico','Adiabático: dT/dn = 0',Neumann(0))]
    if beta is not None:
        if beta > 0:
            raise ValueError('Con normal exterior, beta positivo representa aporte de calor cuando T>20, no enfriamiento. Confirmar el signo antes de usar Robin pasivo.')
        borde=Neumann(0) if beta==0 else Robin(h=-K*beta,T_inf=20)
        casos.append(('beta',f'Robin: beta = {beta:g} 1/m'+(' (adiabático)' if beta==0 else ''),borde))
    return casos


def generar_movies(out, registros):
    """Una animacion conjunta y una por caso; misma escala 20..Tmax en todas.
    Solo se redibuja una figura por instante y luego se recortan los paneles.
    """
    ncasos=len(registros)
    vmax=max(70.,max(float(r['cortes'].max()) for r in registros))
    fig,axes=plt.subplots(ncasos,1,figsize=(9,2*ncasos),layout='constrained',squeeze=False)
    imagenes=[]
    for a,r in zip(axes[:,0],registros):
        im=a.imshow(r['cortes'][0],origin='lower',extent=(0,LARGOS[0],0,LARGOS[2]),
                    aspect='equal',cmap='inferno',vmin=20,vmax=vmax,interpolation='nearest')
        a.set(xlabel='x [m]',ylabel='z [m]',title=r['nombre'])
        imagenes.append(im)
    barra=fig.colorbar(im,ax=list(axes[:,0]),label='Temperatura [°C]',shrink=.85)
    fig.suptitle('Cortes verticales en y = 1,5 m · misma escala de color',fontsize=12)
    fig.set_dpi(85)
    indices=np.unique(np.rint(np.linspace(0,len(registros[0]['t'])-1,49)).astype(int))
    todos=[];individuales=[[] for _ in registros]
    for num,j in enumerate(indices):
        for a,im,r in zip(axes[:,0],imagenes,registros):
            im.set_data(r['cortes'][j])
            a.set_title(f"{r['nombre']} | día {r['t'][j]/DIA:.2f} | centro {r['centro'][j]:.2f} °C",fontsize=10)
        fig.canvas.draw()
        cuadro=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
        todos.append(cuadro)
        renderer=fig.canvas.get_renderer()
        alto=cuadro.height
        for i,a in enumerate(axes[:,0]):
            b=a.get_tightbbox(renderer)
            # Mantener ancho y altura constantes; incluir los ejes y el titulo.
            caja=(max(0,int(b.x0)-6),max(0,int(alto-b.y1)-8),min(cuadro.width,int(b.x1)+8),min(alto,int(alto-b.y0)+8))
            panel=cuadro.crop(caja)
            individuales[i].append(panel)
        if num%16==0:
            print(f'Animaciones: {num+1}/{len(indices)} cuadros preparados',flush=True)
    plt.close(fig)
    todos[0].save(out/'movie_cortes_comparacion.gif',save_all=True,append_images=todos[1:],duration=160,loop=0)
    for r,frames in zip(registros,individuales):
        # Adjuntar una leyenda explicita de la escala comun para cada pelicula.
        from PIL import ImageDraw
        normalizados=[]
        w=max(f.width for f in frames);h=max(f.height for f in frames)
        for frame in frames:
            p=Image.new('RGB',(w,h+28),'white');p.paste(frame,(0,0))
            ImageDraw.Draw(p).text((8,h+5),f'Escala fija: 20 a {vmax:.1f} C. Plano y=1.5 m.',fill='black')
            normalizados.append(p)
        normalizados[0].save(out/f"movie_{r['clave']}.gif",save_all=True,append_images=normalizados[1:],duration=160,loop=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--beta',type=float,default=BETA,help='Coeficiente con normal exterior [1/m], negativo para conveccion')
    parser.add_argument('--refinar',action='store_true',help='Repite con celdas de 5 cm para comprobar sensibilidad')
    parser.add_argument('--sin-movies',action='store_true',help='Solo calcula y genera graficos estaticos, para ejecucion rapida')
    parser.add_argument('--salida',type=Path,default=None,help='Carpeta explicita; por defecto se crea una sesion nueva para no sobrescribir archivos abiertos')
    args=parser.parse_args()
    casos=configurar_casos(args.beta)
    out=args.salida or (Path(__file__).resolve().parent/'resultados_ejercicio_jueves_3D'
                       /datetime.now().strftime('ejecucion_%Y%m%d_%H%M%S_%f'))
    out.mkdir(parents=True,exist_ok=True)
    material=Material('Hormigon',k=K,rho=RHO,cp=CP)
    n=tuple(2*v for v in CELDAS) if args.refinar else CELDAS
    malla=Malla3D.uniforme(LARGOS,n,material)
    bordes_casos=[{c:b for c in CARAS} for _,_,b in casos]
    limites=[paso_tiempo_maximo(malla,b) for b in bordes_casos]
    dt_obj=min(DT,.8*min(limites))
    dt=T_FINAL/int(np.ceil(T_FINAL/dt_obj))
    # Fuente promediada exactamente en cada paso: evita error al integrar el pulso.
    # La conduccion se sigue resolviendo con Euler explicito.
    def fuente(t,T):
        return float((integral_fuente(t+dt)-integral_fuente(t))/dt)
    resumen={'dimensiones_m':LARGOS,'centro_m':tuple(L/2 for L in LARGOS),
             'celdas':n,'pasos_espaciales_m':malla.d,'alpha_m2_s':material.alpha,
             'Qmax_W_m3':QMAX,'duracion_pulso_dias':DURACION_PULSO/DIA,
             'duracion_simulacion_dias':T_FINAL/DIA,'dt_s':dt,'beta_1_m':args.beta,
             'nota':'El Qmax de otro dominio NO impone aqui una temperatura central de 70 C.',
             'casos':[]}
    print(f'Dominio {LARGOS} m; celdas {n}; alpha={material.alpha:.7g} m2/s; dt={dt:g} s',flush=True)
    fig,ax=plt.subplots(2,1,figsize=(10,8),layout='constrained')
    colores=['#287bb5','#d97922','#258b62','#9e4eab']
    cortes=[];nombres=[];registros=[]
    for (clave,nombre,_),bordes,limite,color in zip(casos,bordes_casos,limites,colores):
        res=resolver(malla,T_INICIAL,bordes,T_FINAL,dt=dt,fuente=fuente,n_guardados=288)
        centro=temperatura_centro(res)
        E=res.energia();energiag=integral_fuente(res.t)*np.prod(LARGOS)
        intercambio=E-E[0]-energiag
        tmin=res.T.min(axis=(1,2,3));tmax=res.T.max(axis=(1,2,3))
        tmedia=res.T.mean(axis=(1,2,3))
        fila={'caso':clave,'condicion':nombre,'dt_critico_coeficientes_positivos_s':limite,
              'Tmax_centro_C':float(centro.max()),'dia_max_centro':float(res.t[centro.argmax()]/DIA),
              'T_centro_final_C':float(centro[-1]),'Tmax_global_guardada_C':float(res.T.max()),
              'energia_generada_final_J':float(energiag[-1]),
              'intercambio_neto_por_bordes_J':float(intercambio[-1]),
              'T_min_final_C':float(tmin[-1]),'T_max_final_C':float(tmax[-1]),
              'T_media_final_C':float(tmedia[-1]),'energia_inicial_J':float(E[0]),'energia_final_J':float(E[-1])}
        if clave=='adiabatico':
            exacta=T_INICIAL+integral_fuente(res.t)/(RHO*CP)
            error=float(np.max(np.abs(res.T-exacta[:,None,None,None])))
            fila['error_adiabatico_C']=error
            assert error < 1e-8
            assert np.max(np.abs(intercambio)) < 1e-5*max(1.,energiag[-1])
        resumen['casos'].append(fila)
        print(f'{nombre}: limite dt={limite:.1f} s; max centro={centro.max():.3f} C (dia {res.t[centro.argmax()]/DIA:.3f}); final={centro[-1]:.3f} C',flush=True)
        with (out/f'{clave}_historia.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.writer(f);w.writerow(['t_s','t_dias','Q_W_m3','T_centro_C','T_min_C','T_max_C','T_media_C','energia_J','energia_generada_J','intercambio_bordes_J'])
            w.writerows(zip(res.t,res.t/DIA,QMAX*fraccion_fuente(res.t),centro,tmin,tmax,tmedia,E,energiag,intercambio))
        ax[0].plot(res.t/DIA,centro,label=nombre,color=color,lw=2)
        ax[1].plot(res.t/DIA,intercambio/1e6,label=nombre,color=color,lw=2)
        # Corte vertical por el plano medio y=1,5 m; promedio de capas centrales.
        js=[n[1]//2-1,n[1]//2] if n[1]%2==0 else [n[1]//2]
        cortes_t=np.take(res.T,js,axis=2).mean(axis=2).transpose(0,2,1)
        cortes.append(cortes_t[-1])
        nombres.append(nombre)
        registros.append(dict(clave=clave,nombre=nombre,t=res.t,centro=centro,
                              minimo=tmin,maximo=tmax,media=tmedia,energia=E,
                              generada=energiag,intercambio=intercambio,cortes=cortes_t))
        # Campo 3D completo para volver a analizar cualquier punto o plano.
        np.savez_compressed(out/f'{clave}_campos_3D.npz',t_s=res.t,T_C=res.T,
                            x_m=malla.centros[0],y_m=malla.centros[1],z_m=malla.centros[2])
    ax[0].axhline(70,color='gray',ls='--',label='Referencia de 70 °C (no impuesta)')
    ax[0].set(title='Historia del centro (1,5; 1,5; 0,30) m | Qmax = 1300 W/m³',xlabel='Tiempo [días]',ylabel='Temperatura [°C]')
    ax[1].set(title='Energía intercambiada por las caras: positiva = entrada neta',xlabel='Tiempo [días]',ylabel='Intercambio acumulado [MJ]')
    for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
    fig.savefig(out/'comparacion_condiciones.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(len(cortes),1,figsize=(10,2.2*len(cortes)),layout='constrained',squeeze=False)
    vmax=max(np.max(c) for c in cortes)
    for a,c,nombre in zip(axes[:,0],cortes,nombres):
        im=a.imshow(c,origin='lower',extent=(0,3,0,.6),aspect='equal',cmap='inferno',vmin=20,vmax=vmax)
        a.set(title=nombre+' | día 2',xlabel='x [m]',ylabel='z [m]')
    fig.colorbar(im,ax=list(axes[:,0]),label='Temperatura [°C]',shrink=.8)
    fig.savefig(out/'cortes_verticales_dia2.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(10,10),layout='constrained')
    t=registros[0]['t']
    axs[0].plot(t/DIA,QMAX*fraccion_fuente(t),color='#c36c19',lw=2)
    axs[0].set(title='Fuente volumétrica común a los cuatro casos',ylabel='Q(t) [W/m³]')
    for r,color in zip(registros,colores):
        axs[1].plot(r['t']/DIA,r['centro'],color=color,label=r['nombre'])
        axs[1].fill_between(r['t']/DIA,r['minimo'],r['maximo'],color=color,alpha=.10)
        axs[2].plot(r['t']/DIA,r['energia']/1e6,color=color,label=r['nombre'])
    axs[1].axhline(70,color='gray',ls='--')
    axs[1].set(title='Centro y rango mínimo–máximo de las celdas (sombreado)',ylabel='Temperatura [°C]')
    axs[2].set(title='Energía interna referida a 0 °C',ylabel='Energía [MJ]')
    for a in axs:a.set_xlabel('Tiempo [días]');a.grid(alpha=.2)
    axs[1].legend(fontsize=8)
    fig.savefig(out/'fuente_temperaturas_energia.png',dpi=150);plt.close(fig)
    campos_csv=sorted(set().union(*(fila.keys() for fila in resumen['casos'])))
    with (out/'resumen_casos.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=campos_csv);w.writeheader();w.writerows(resumen['casos'])
    if args.beta is None:
        resumen['pendiente']='Falta valor y signo de beta: cuarto caso no ejecutado.'
        print('PENDIENTE: indicar BETA o --beta para ejecutar el cuarto caso.',flush=True)
    (out/'resumen.json').write_text(json.dumps(resumen,indent=2,ensure_ascii=False),encoding='utf-8')
    if not args.sin_movies:
        generar_movies(out,registros)
    (out/'LEEME.txt').write_text(
        'RESULTADOS DEL EJERCICIO 3D\n'
        'movie_cortes_comparacion.gif: cuatro cortes verticales sincronizados.\n'
        'movie_dirichlet_20.gif, movie_dirichlet_30.gif, movie_adiabatico.gif, movie_beta.gif: peliculas individuales.\n'
        'Las peliculas cubren 0 a 2 dias, en el plano y=1,5 m; escala de temperatura comun y fija.\n'
        'beta=0 coincide con adiabatico: no hay intercambio por las caras.\n'
        'comparacion_condiciones.png: centro e intercambio acumulado por bordes.\n'
        'fuente_temperaturas_energia.png: fuente, temperaturas y energia interna.\n'
        'cortes_verticales_dia2.png: cortes finales.\n'
        '*_historia.csv: temperaturas, fuente y balance de energia para cada caso.\n'
        '*_campos_3D.npz: todas las temperaturas 3D guardadas, tiempos y coordenadas.\n'
        'resumen_casos.csv y resumen.json: parametros y resultados resumidos.\n'
        'Los campos son temperaturas en centros de celda, no en la superficie exacta.\n'
        'La condicion Dirichlet se impone en la cara: su temperatura no tiene que coincidir con la de la celda adyacente.\n',encoding='utf-8')
    print(f'Resultados guardados en {out}',flush=True)


if __name__=='__main__':main()
