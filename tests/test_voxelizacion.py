"""Verificacion geometrica: no depende de temperaturas ni de Blender."""
from pathlib import Path
import numpy as np
import pytest
from mcplan import bim
from mcplan.voxelizacion import MallaVertida,clasificar_puntos,extraer_mallas,voxelizar_mallas,voxelizar_ifc

RAIZ=Path(__file__).resolve().parents[1]


def caja(vid=1,origen=(0,0,0),lados=(1,1,1)):
    v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                [0,0,1],[1,0,1],[1,1,1],[0,1,1]],dtype=float)
    v=v*np.array(lados)+np.array(origen)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],
                [0,1,5],[0,5,4],[3,7,6],[3,6,2],
                [0,4,7],[0,7,3],[1,2,6],[1,6,5]])
    return MallaVertida(vid,f'Caja {vid}',f'guid-{vid}',v,f,float(np.prod(lados)),{'Inicio_h':vid*6})


def test_interior_exterior_y_superficie():
    dentro,borde=clasificar_puntos(caja(),np.array([[.5,.5,.5],[1.1,.5,.5],[0,.5,.5],[0,0,0]]))
    assert dentro.tolist()==[True,False,False,False]
    assert borde.tolist()==[False,False,True,True]


def test_cavidad_no_se_rellena_como_si_fuera_caja():
    exterior=caja(lados=(3,3,3))
    interior=caja(origen=(1,1,1))
    hueca=MallaVertida(1,'Hueca','hueca',np.vstack([exterior.vertices,interior.vertices]),
                      np.vstack([exterior.triangulos,interior.triangulos[:,::-1]+8]),26)
    dentro,borde=clasificar_puntos(hueca,np.array([[.5,1.5,1.5],[1.5,1.5,1.5],[3.1,1,1]]))
    assert dentro.tolist()==[True,False,False]
    assert not borde.any()
    g=voxelizar_mallas([hueca],.5)
    assert g.resumen_volumenes()[0]['volumen_voxel_m3']==pytest.approx(26)


def test_caja_trasladada_con_volumen_exacto_y_ejes_xyz():
    g=voxelizar_mallas([caja(7,(-3,4,2),(1,2,.5))],(.25,.5,.125))
    assert g.vertida_id.shape==(4,4,4)
    assert np.all(g.vertida_id==7)
    assert g.origen==pytest.approx([-3,4,2])
    assert g.centros[2][0]==pytest.approx(2.0625)
    assert g.resumen_volumenes()[0]['volumen_voxel_m3']==pytest.approx(1.)
    assert g.vertidas[7]['Inicio_h']==42


def test_rotacion_se_lee_desde_triangulos_no_desde_bounding_box():
    m=caja()
    ang=np.pi/4
    R=np.array([[np.cos(ang),-np.sin(ang),0],[np.sin(ang),np.cos(ang),0],[0,0,1]])
    m=MallaVertida(1,'Rotada','rot',m.vertices@R.T,m.triangulos,1)
    d,b=clasificar_puntos(m,np.array([[0,.7,.5],[.65,.1,.5]]))
    assert d.tolist()==[True,False]


def test_junta_exacta_tiene_regla_determinista_sin_duplicar_celdas():
    a=caja(1,lados=(1,1,.5))
    b=caja(2,origen=(0,0,.5),lados=(1,1,.5))
    g=voxelizar_mallas([b,a],(1,1,1)) # unico centro exactamente en la junta
    assert g.vertida_id.item()==1
    assert g.centros_compartidos==1
    assert sum(f['celdas'] for f in g.resumen_volumenes())==1


def test_traslape_real_no_se_oculta_con_el_orden_de_las_vertidas():
    with pytest.raises(ValueError,match='Traslape real'):
        voxelizar_mallas([caja(1),caja(2,origen=(.25,0,0))],.25)


@pytest.mark.parametrize('paso',[0,-.1,float('nan')])
def test_paso_invalido_se_rechaza(paso):
    with pytest.raises(ValueError,match='paso'):
        voxelizar_mallas([caja()],paso)


def test_malla_abierta_se_rechaza():
    m=caja()
    with pytest.raises(ValueError,match='cerrada'):
        MallaVertida(1,'Abierta','bad',m.vertices,m.triangulos[:-1],1)


def test_limite_de_memoria_se_comprueba_antes_de_crear_grilla():
    with pytest.raises(ValueError,match='limite'):
        voxelizar_mallas([caja()],.001,max_celdas=1000)


def test_npz_guarda_etiquetas_y_metadata_sin_pickle(tmp_path):
    g=voxelizar_mallas([caja(8)],.25)
    ruta=tmp_path/'g.npz';g.guardar(ruta)
    with np.load(ruta,allow_pickle=False) as f:
        assert np.array_equal(f['vertida_id'],g.vertida_id)
        assert 'guid-8' in str(f['metadata_json'])


def test_ifc_real_conserva_cinco_ids_y_volumenes_por_vertida():
    ruta=RAIZ/'fundaciones'/'fundacion.ifc'
    mallas=extraer_mallas(ruta)
    assert [m.vertida_id for m in mallas]==[1,2,3,4,5]
    g=voxelizar_ifc(ruta,(.4,.4,.05))
    assert set(np.unique(g.vertida_id))=={0,1,2,3,4,5}
    filas=g.resumen_volumenes()
    exacto=sum(f['volumen_ifc_m3'] for f in filas)
    voxel=sum(f['volumen_voxel_m3'] for f in filas)
    assert exacto==pytest.approx(636.77522695,rel=1e-7)
    assert abs(voxel-exacto)/exacto<.03
    assert all(f['celdas']>0 for f in filas)
    assert g.vertidas[3]['Inicio_h']==12


def test_extraccion_respeta_traslacion_guardada_en_ifc(tmp_path):
    import ifcopenshell.api.geometry
    modelo=bim.abrir(RAIZ/'fundaciones'/'fundacion.ifc')
    elemento=bim.elementos_vertida(modelo)[-1]
    antes=extraer_mallas(modelo)[-1].limites[0].copy()
    import ifcopenshell.util.placement
    matrix=ifcopenshell.util.placement.get_local_placement(elemento.ObjectPlacement).copy()
    matrix[:3,3]+=np.array([2.,3.,4.])
    ifcopenshell.api.geometry.edit_object_placement(modelo,product=elemento,matrix=matrix)
    ruta=tmp_path/'movido.ifc';modelo.write(str(ruta))
    despues=extraer_mallas(ruta)[-1].limites[0]
    assert despues-antes==pytest.approx([2,3,4])
