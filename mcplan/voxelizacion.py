"""Primer voxelizador IFC: triangulos cerrados -> grilla cartesiana etiquetada.

Metodo: centro de celda. Se usa el angulo solido orientado de TODOS los
triangulos para decidir pertenencia (numero de enrollamiento). No supone
que la fundacion sea un prisma ni que los objetos sean convexos.

0 = fuera del hormigon modelado (no implica automaticamente aire o suelo).
ID positivo = VertidaID leido del IFC. Ejes del arreglo: (x,y,z), como solver3d.
No activa vertidas ni resuelve temperatura: solo prepara geometria y datos.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import json
import numpy as np
import ifcopenshell.geom
import ifcopenshell.util.shape
from mcplan import bim


@dataclass
class MallaVertida:
    vertida_id: int
    nombre: str
    global_id: str
    vertices: np.ndarray
    triangulos: np.ndarray
    volumen_ifc_m3: float
    propiedades: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.vertida_id <= 0:
            raise ValueError('VertidaID debe ser positivo; 0 se reserva para exterior.')
        self.vertices=np.asarray(self.vertices,dtype=float)
        self.triangulos=np.asarray(self.triangulos,dtype=np.int64)
        if self.vertices.ndim!=2 or self.vertices.shape[1]!=3 or not np.isfinite(self.vertices).all():
            raise ValueError('Vertices: se espera un arreglo finito (n,3).')
        if self.triangulos.ndim!=2 or self.triangulos.shape[1]!=3 or not len(self.triangulos):
            raise ValueError('Triangulos: se espera un arreglo no vacio (n,3).')
        if self.triangulos.min()<0 or self.triangulos.max()>=len(self.vertices):
            raise ValueError('Los indices de triangulos no corresponden a los vertices.')
        # Soldar vertices coincidentes a escala de 1 nanometro; mantiene la
        # posicion original y permite validar mallas con vertices duplicados.
        _,primeros,inversa=np.unique(np.round(self.vertices,9),axis=0,return_index=True,return_inverse=True)
        self.vertices=self.vertices[primeros]
        self.triangulos=inversa[self.triangulos]
        tri=self.vertices[self.triangulos]
        if np.any(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)<1e-14):
            raise ValueError(f'{self.nombre}: triangulo degenerado.')
        aristas=np.concatenate([self.triangulos[:,[0,1]],self.triangulos[:,[1,2]],self.triangulos[:,[2,0]]])
        _,indices,cuentas=np.unique(np.sort(aristas,axis=1),axis=0,return_inverse=True,return_counts=True)
        signos=np.where(aristas[:,0]<aristas[:,1],1,-1)
        orientacion=np.bincount(indices,weights=signos)
        if np.any(cuentas!=2) or np.any(orientacion!=0):
            raise ValueError(f'{self.nombre}: la malla debe ser cerrada y tener triangulos orientados consistentemente.')
        if not np.isfinite(self.volumen_ifc_m3) or self.volumen_ifc_m3<=0:
            raise ValueError(f'{self.nombre}: el volumen IFC debe ser positivo.')

    @property
    def limites(self):
        return self.vertices.min(axis=0),self.vertices.max(axis=0)


def extraer_mallas(modelo_o_ruta) -> list[MallaVertida]:
    """Lee el IFC REAL: triangulacion y posiciones mundiales en metros.
    Conserva GlobalId, ID y propiedades. No reconstruye formas desde YAML.
    """
    modelo=bim.abrir(modelo_o_ruta)
    propiedades={v['GlobalId']:v for v in bim.leer_vertidas(modelo)}
    ajustes=ifcopenshell.geom.settings()
    ajustes.set('use-world-coords',True)
    mallas=[]
    for elemento in bim.elementos_vertida(modelo):
        props=propiedades[elemento.GlobalId]
        forma=ifcopenshell.geom.create_shape(ajustes,elemento)
        mallas.append(MallaVertida(
            vertida_id=bim._vertida_id_numero(props['VertidaID']),
            nombre=elemento.Name or '',global_id=elemento.GlobalId,
            vertices=np.asarray(forma.geometry.verts,dtype=float).reshape(-1,3),
            triangulos=np.asarray(forma.geometry.faces,dtype=np.int64).reshape(-1,3),
            volumen_ifc_m3=abs(float(ifcopenshell.util.shape.get_volume(forma.geometry))),
            propiedades=props,
        ))
    if not mallas:
        raise ValueError('No hay vertidas con el property set esperado en el IFC.')
    ids=[m.vertida_id for m in mallas]
    if len(set(ids))!=len(ids):
        raise ValueError('VertidaID duplicado: cada solido debe tener un ID unico.')
    return sorted(mallas,key=lambda m:m.vertida_id)


def clasificar_puntos(malla:MallaVertida,puntos:np.ndarray,tolerancia_m=1e-8):
    """Devuelve (interior_estricto, superficie) para un lote de puntos.
    Un punto interior acumula +/-4 pi de angulo solido; exterior, cero.
    Superficie se comprueba aparte para resolver centros en juntas sin azar.
    Las cavidades requieren cascaras interiores correctamente orientadas.
    """
    p=np.asarray(puntos,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():
        raise ValueError('Se esperan puntos finitos de forma (n,3).')
    if tolerancia_m<=0:
        raise ValueError('La tolerancia geometrica debe ser positiva.')
    total=np.zeros(len(p));superficie=np.zeros(len(p),dtype=bool)
    for tri in malla.vertices[malla.triangulos]:
        a,b,c=tri[0]-p,tri[1]-p,tri[2]-p
        la,lb,lc=(np.linalg.norm(v,axis=1) for v in (a,b,c))
        numerador=np.einsum('ij,ij->i',a,np.cross(b,c))
        denominador=(la*lb*lc+np.einsum('ij,ij->i',a,b)*lc
                     +np.einsum('ij,ij->i',b,c)*la+np.einsum('ij,ij->i',c,a)*lb)
        total+=2*np.arctan2(numerador,denominador)
        u,v=tri[1]-tri[0],tri[2]-tri[0]
        normal=np.cross(u,v);longitud=np.linalg.norm(normal)
        cercanos=np.flatnonzero(np.abs((p-tri[0])@normal)<=tolerancia_m*longitud)
        if len(cercanos):
            w=p[cercanos]-tri[0]
            uu,uv,vv=np.dot(u,u),np.dot(u,v),np.dot(v,v)
            wu,wv=w@u,w@v;det=uu*vv-uv**2
            s=(vv*wu-uv*wv)/det;t=(uu*wv-uv*wu)/det
            eps=tolerancia_m/max(np.linalg.norm(u),np.linalg.norm(v))
            en_cara=(s>=-eps)&(t>=-eps)&(s+t<=1+eps)
            superficie[cercanos[en_cara]]=True
    interior=(np.abs(total)>2*np.pi)&~superficie
    return interior,superficie


@dataclass
class GrillaVoxeles:
    origen: np.ndarray
    paso: tuple[float,float,float]
    vertida_id: np.ndarray
    vertidas: dict[int,dict]
    centros_compartidos: int = 0

    @property
    def centros(self):
        return tuple(self.origen[e]+(np.arange(n)+.5)*self.paso[e]
                     for e,n in enumerate(self.vertida_id.shape))

    @property
    def volumen_celda(self):
        return float(np.prod(self.paso))

    def resumen_volumenes(self):
        filas=[]
        for vid,props in sorted(self.vertidas.items()):
            celdas=int(np.count_nonzero(self.vertida_id==vid))
            voxel=celdas*self.volumen_celda;exacto=props['volumen_ifc_m3']
            filas.append(dict(vertida_id=vid,nombre=props['Nombre'],celdas=celdas,
                              volumen_ifc_m3=exacto,volumen_voxel_m3=voxel,
                              error_relativo_porcentaje=100*(voxel-exacto)/exacto))
        return filas

    def guardar(self,ruta):
        """NPZ sin pickle; incluye coordenadas, etiquetas y metadata JSON."""
        metadata=json.dumps({'vertidas':self.vertidas,'centros_compartidos':self.centros_compartidos},ensure_ascii=False)
        np.savez_compressed(ruta,origen=self.origen,paso=np.array(self.paso),
                            vertida_id=self.vertida_id,metadata_json=np.array(metadata))

    @classmethod
    def cargar(cls,ruta):
        """Recupera la grilla y sus propiedades sin permitir pickle."""
        with np.load(ruta,allow_pickle=False) as f:
            origen=np.array(f['origen'],dtype=float)
            paso=np.array(f['paso'],dtype=float)
            ids=np.array(f['vertida_id'])
            metadata=json.loads(str(f['metadata_json']))
        if origen.shape!=(3,) or paso.shape!=(3,) or not np.isfinite(origen).all() or not np.isfinite(paso).all() or np.any(paso<=0):
            raise ValueError('Origen y pasos invalidos en el archivo de grilla.')
        if ids.ndim!=3 or not np.issubdtype(ids.dtype,np.integer) or np.any(ids<0):
            raise ValueError('vertida_id debe ser un arreglo entero 3D no negativo.')
        vertidas={int(k):v for k,v in metadata['vertidas'].items()}
        if not (set(np.unique(ids))-{0}).issubset(vertidas):
            raise ValueError('Hay celdas sin propiedades de vertida en la grilla.')
        return cls(origen,tuple(float(x) for x in paso),ids,vertidas,
                   int(metadata.get('centros_compartidos',0)))


def voxelizar_mallas(mallas:list[MallaVertida],paso=.20,*,max_celdas=2_000_000,lote=4096):
    """Centros interiores o sobre superficie -> etiqueta de vertida.

    Grilla anclada al minimo global IFC; paso fijo, caja ampliada con ceil.
    Centros estrictamente interiores a dos vertidas => ERROR de traslape.
    En juntas exactas, interior estricto prevalece sobre superficie; si solo
    hay superficies compartidas, gana el menor ID (regla documentada).
    Refinar no garantiza error monotono por el escalonado de la superficie.
    """
    if not mallas:
        raise ValueError('No hay mallas para voxelizar.')
    d=np.broadcast_to(np.asarray(paso,dtype=float),(3,)).copy()
    if not np.isfinite(d).all() or np.any(d<=0) or lote<=0:
        raise ValueError('El paso debe ser finito y positivo, y el lote mayor que cero.')
    if len({m.vertida_id for m in mallas})!=len(mallas):
        raise ValueError('VertidaID duplicado.')
    origen=np.min([m.limites[0] for m in mallas],axis=0)
    superior=np.max([m.limites[1] for m in mallas],axis=0)
    forma=np.maximum(1,np.ceil((superior-origen)/d-1e-10).astype(int))
    if np.prod(forma,dtype=np.int64)>max_celdas:
        raise ValueError(f'La grilla requiere {np.prod(forma)} celdas (limite {max_celdas}). Use un paso mayor o aumente max_celdas conscientemente.')
    etiquetas=np.zeros(tuple(forma),dtype=np.int32)
    estricto=np.zeros(tuple(forma),dtype=bool)
    centros=tuple(origen[e]+(np.arange(n)+.5)*d[e] for e,n in enumerate(forma))
    compartidos=0;metadatos={}
    for m in sorted(mallas,key=lambda m:m.vertida_id):
        lo,hi=m.limites
        indices=[np.flatnonzero((x>=lo[e]-1e-8)&(x<=hi[e]+1e-8)) for e,x in enumerate(centros)]
        f=tuple(len(i) for i in indices)
        total=int(np.prod(f))
        for inicio in range(0,total,lote):
            locales=np.unravel_index(np.arange(inicio,min(inicio+lote,total)),f)
            glob=tuple(indices[e][locales[e]] for e in range(3))
            p=np.column_stack([centros[e][glob[e]] for e in range(3)])
            dentro,borde=clasificar_puntos(m,p)
            anterior=etiquetas[glob];era_estricto=estricto[glob]
            if np.any(dentro&era_estricto):
                otros=np.unique(anterior[dentro&era_estricto]).tolist()
                raise ValueError(f'Traslape real: centros interiores a vertida {m.vertida_id} y a {otros}.')
            compartidos+=int(np.count_nonzero(borde&(anterior!=0)))
            tomar=dentro|(borde&(anterior==0))
            anterior[tomar]=m.vertida_id
            era_estricto[tomar]=dentro[tomar]
            etiquetas[glob]=anterior;estricto[glob]=era_estricto
        metadatos[m.vertida_id]={**m.propiedades,'Nombre':m.nombre,'GlobalId':m.global_id,
                                'volumen_ifc_m3':m.volumen_ifc_m3}
    return GrillaVoxeles(origen,tuple(float(x) for x in d),etiquetas,metadatos,compartidos)


def voxelizar_ifc(ruta,paso=.20,**opciones):
    return voxelizar_mallas(extraer_mallas(ruta),paso,**opciones)
