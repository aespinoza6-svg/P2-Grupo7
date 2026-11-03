"""E07: leer IFC -> grilla con vertida_id -> volumenes y evidencia visual.

Desde la raiz y con el entorno instalado:
python scripts/07_voxelizar_ifc.py
python scripts/07_voxelizar_ifc.py --dx .1 --dy .1 --dz .05
python scripts/07_voxelizar_ifc.py --comparar

Por defecto se usa una celda de 20 x 20 x 5 cm: dz=5 cm permite representar
las juntas a z=0.55,1.10,2.20,3.30,4.30 sin que centros caigan sobre ellas.
No modifica el IFC ni los parametros del solver. No activa vertidas todavia.
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
from mcplan.voxelizacion import extraer_mallas,voxelizar_mallas

RAIZ=Path(__file__).resolve().parents[1]


def exportar_vtk(grilla,ruta):
    """VTK legacy: x cambia mas rapido, etiquetas almacenadas como CELL_DATA.
    Abrir en ParaView/PyVista; 0 es exterior, filtrar IDs >=1 para hormigon.
    """
    nx,ny,nz=grilla.vertida_id.shape
    with ruta.open('w',encoding='ascii') as f:
        f.write('# vtk DataFile Version 3.0\nIFC voxelizado - VertidaID\nASCII\nDATASET STRUCTURED_POINTS\n')
        f.write(f'DIMENSIONS {nx+1} {ny+1} {nz+1}\n')
        f.write('ORIGIN '+' '.join(map(str,grilla.origen))+'\n')
        f.write('SPACING '+' '.join(map(str,grilla.paso))+'\n')
        f.write(f'CELL_DATA {grilla.vertida_id.size}\nSCALARS vertida_id int 1\nLOOKUP_TABLE default\n')
        np.savetxt(f,grilla.vertida_id.ravel(order='F'),fmt='%d')


def visualizar(g,out):
    ids=[0]+sorted(g.vertidas)
    paleta=['#f2f3f5','#287bb5','#d97922','#258b62','#985bad','#cd4c44','#b3942a']
    colores=[paleta[i%len(paleta)] for i in range(len(ids))]
    cmap=ListedColormap(colores);norm=BoundaryNorm(np.arange(len(ids)+1)-.5,len(ids))
    valores=np.searchsorted(ids,g.vertida_id)
    x,y,z=g.centros
    ix=int(np.argmin(abs(x)));iy=int(np.argmin(abs(y)));iz=int(np.argmin(abs(z-1.65)))
    limites=[(g.origen[e],g.origen[e]+g.paso[e]*g.vertida_id.shape[e]) for e in range(3)]
    fig,axs=plt.subplots(1,3,figsize=(14,5),layout='constrained')
    vistas=[(valores[:,iy,:].T,(*limites[0],*limites[2]),'x [m]','z [m]',f'Corte vertical y = {y[iy]:.3f} m'),
            (valores[ix,:,:].T,(*limites[1],*limites[2]),'y [m]','z [m]',f'Corte vertical x = {x[ix]:.3f} m'),
            (valores[:,:,iz].T,(*limites[0],*limites[1]),'x [m]','y [m]',f'Planta z = {z[iz]:.3f} m')]
    for a,(campo,extension,xlabel,ylabel,titulo) in zip(axs,vistas):
        a.imshow(campo,origin='lower',extent=extension,cmap=cmap,norm=norm,interpolation='nearest')
        a.set(xlabel=xlabel,ylabel=ylabel,title=titulo)
    fig.legend(handles=[Patch(color=colores[i],label='Exterior' if vid==0 else f'Vertida {vid}') for i,vid in enumerate(ids)],loc='outside lower center',ncol=min(6,len(ids)),fontsize=9)
    fig.suptitle('IFC → vóxeles: cada celda conserva la identidad de su vertida')
    fig.savefig(out/'voxelizacion_cortes.png',dpi=150);plt.close(fig)
    ocupadas=g.vertida_id>0
    interior=ocupadas.copy()
    for eje in range(3):
        interior&=np.roll(ocupadas,1,axis=eje)&np.roll(ocupadas,-1,axis=eje)
        inicio=[slice(None)]*3;fin=[slice(None)]*3;inicio[eje]=0;fin[eje]=-1
        interior[tuple(inicio)]=False;interior[tuple(fin)]=False
    puntos=np.column_stack(np.nonzero(ocupadas&~interior))
    if len(puntos)>14000:
        puntos=puntos[np.random.default_rng(7).choice(len(puntos),14000,replace=False)]
    fig=plt.figure(figsize=(10,7),layout='constrained');a=fig.add_subplot(projection='3d')
    xyz=np.column_stack([g.centros[e][puntos[:,e]] for e in range(3)])
    for indice,vid in enumerate(ids[1:],1):
        elegir=g.vertida_id[tuple(puntos.T)]==vid
        if elegir.any():
            a.scatter(*xyz[elegir].T,s=5,marker='s',color=colores[indice],alpha=.8,label=f'Vertida {vid}',rasterized=True)
    a.set(xlabel='x [m]',ylabel='y [m]',zlabel='z [m]',title='Centros de vóxeles en la superficie exterior · color por vertida')
    a.set_box_aspect([b-a for a,b in limites]);a.view_init(elev=24,azim=-48);a.legend(loc='upper left',fontsize=9)
    fig.savefig(out/'voxelizacion_3d.png',dpi=150);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ifc',type=Path,default=RAIZ/'fundaciones'/'fundacion.ifc')
    p.add_argument('--dx',type=float,default=.2);p.add_argument('--dy',type=float,default=.2)
    p.add_argument('--dz',type=float,default=.05)
    p.add_argument('--salida',type=Path,default=RAIZ/'resultados'/'e07_voxelizacion')
    p.add_argument('--comparar',action='store_true',help='Compara pasos en x/y doble, actual y mitad, manteniendo dz')
    p.add_argument('--sin-figuras',action='store_true')
    args=p.parse_args();out=args.salida
    mallas=extraer_mallas(args.ifc)
    g=voxelizar_mallas(mallas,(args.dx,args.dy,args.dz),max_celdas=4_000_000)
    filas=g.resumen_volumenes();exacto=sum(f['volumen_ifc_m3'] for f in filas)
    aproximado=sum(f['volumen_voxel_m3'] for f in filas)
    out.mkdir(parents=True,exist_ok=True)
    g.guardar(out/'grilla_vertidas.npz');exportar_vtk(g,out/'grilla_vertidas.vtk')
    with (out/'volumenes_por_vertida.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(filas[0]));w.writeheader();w.writerows(filas)
    print(f'IFC: {args.ifc}\nGrilla (x,y,z): {g.vertida_id.shape}; paso: {g.paso} m',flush=True)
    print(' ID   celdas     IFC [m3]   voxeles [m3]    error [%]',flush=True)
    for fila in filas:
        print(f"{fila['vertida_id']:3d} {fila['celdas']:9d} {fila['volumen_ifc_m3']:12.3f} {fila['volumen_voxel_m3']:13.3f} {fila['error_relativo_porcentaje']:12.3f}",flush=True)
    print(f'TOTAL: IFC={exacto:.3f} m3; voxeles={aproximado:.3f} m3; error={100*(aproximado-exacto)/exacto:+.3f} %',flush=True)
    print(f'Centros sobre superficies compartidas: {g.centros_compartidos}. Exterior = 0.',flush=True)
    info=dict(ifc=str(args.ifc.resolve()),ifc_sha256=hashlib.sha256(args.ifc.read_bytes()).hexdigest(),
              forma=list(g.vertida_id.shape),origen_m=g.origen.tolist(),paso_m=g.paso,
              celdas_hormigon=int(np.count_nonzero(g.vertida_id)),volumen_ifc_m3=exacto,
              volumen_voxel_m3=aproximado,error_relativo_porcentaje=100*(aproximado-exacto)/exacto,
              centros_compartidos=g.centros_compartidos,volumenes_por_vertida=filas,
              metodo='Centro de celda y angulo solido de triangulos IFC cerrados y orientados',
              exterior='ID 0 significa fuera del hormigon; no define todavia aire, suelo ni condicion termica',
              activacion='Solo geometria; Inicio_h y temperatura de colocacion se conservan en metadata, sin activar celdas aun')
    if args.comparar:
        comparacion=[]
        for factor in [2,1,.5]:
            prueba=g if factor==1 else voxelizar_mallas(mallas,(args.dx*factor,args.dy*factor,args.dz),max_celdas=4_000_000)
            vp=sum(r['volumen_voxel_m3'] for r in prueba.resumen_volumenes())
            comparacion.append(dict(dx_m=prueba.paso[0],dy_m=prueba.paso[1],dz_m=prueba.paso[2],
                                   celdas=int(prueba.vertida_id.size),volumen_voxel_m3=vp,
                                   error_relativo_porcentaje=100*(vp-exacto)/exacto))
            print(f'Comparacion {prueba.paso}: volumen={vp:.3f} m3; error={100*(vp-exacto)/exacto:+.3f} %',flush=True)
        info['comparacion_mallas']=comparacion
        with (out/'comparacion_mallas.csv').open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=list(comparacion[0]));w.writeheader();w.writerows(comparacion)
    (out/'resumen.json').write_text(json.dumps(info,indent=2,ensure_ascii=False,default=lambda x:int(x) if isinstance(x,np.integer) else float(x)),encoding='utf-8')
    if not args.sin_figuras:visualizar(g,out)
    print(f'Resultados: {out}',flush=True)


if __name__=='__main__':main()
