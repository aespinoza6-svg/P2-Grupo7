"""Verificacion del acoplamiento: caras expuestas, activacion y energia."""
import numpy as np
import pytest
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import HidratacionExponencial
from mcplan.construccion import ProgramaVertidas
from mcplan.solver1d import Material,Dirichlet,Neumann,Robin
from mcplan.solver3d import Malla3D,resolver,CARAS
from mcplan.solver_vertidas import SistemaActivo,resolver_vertidas


def programa(ids,inicios=None,temperaturas=None,elevacion=0.,tau=1000.,d=(.1,.1,.1)):
    ids=np.asarray(ids,dtype=np.int32)
    if ids.ndim==1:ids=ids[:,None,None]
    vids=set(np.unique(ids))-{0};inicios=inicios or {};temperaturas=temperaturas or {}
    props={int(v):dict(Inicio_h=inicios.get(v,0.)/3600,TemperaturaColocacion_C=temperaturas.get(v,20.),TipoHormigon='M1') for v in vids}
    g=GrillaVoxeles(np.zeros(3),d,ids,props)
    m=Material('M1',1.,1000.,10.)
    h=HidratacionExponencial(elevacion,tau,m.rho,m.cp)
    return ProgramaVertidas(g,{'M1':h}),{'M1':m}


def comprobar_balance(r):
    for fila in r.balance:
        escala=max(1,abs(fila['energia_colocacion_J'])+abs(fila['energia_hidratacion_J'])+abs(fila['energia_bordes_J']))
        assert abs(fila['residuo_balance_J'])/escala<1e-11


def test_caja_totalmente_activa_reproduce_solver3d():
    p,mat=programa(np.ones((8,3,2)),temperaturas={1:40.})
    bordes={'x-':Dirichlet(60.),'x+':Robin(15.,5.,.05)}
    r=resolver_vertidas(p,mat,bordes,100.,dt_s=2.,tiempos_guardado_s=[0,50,100],guardar_campos=True)
    m=Malla3D.uniforme((.8,.3,.2),(8,3,2),mat['M1'])
    ex=resolver(m,40.,bordes,100.,dt=2.,tiempos_guardado=[0,50,100])
    assert np.allclose(r.T,ex.T,rtol=0,atol=1e-10)
    comprobar_balance(r)


def test_exterior_vacio_no_absorbe_calor_y_se_conserva_energia():
    p,mat=programa([1,0,1],temperaturas={1:60.})
    r=resolver_vertidas(p,mat,{},100.,max_dt_s=10,tiempos_guardado_s=[0,100])
    assert r.T_final[[0,2],0,0]==pytest.approx([60.,60.])
    assert np.isnan(r.T_final[1,0,0])
    comprobar_balance(r)


def test_caras_de_un_hueco_reciben_conveccion():
    p,mat=programa([1,0,1],temperaturas={1:60.})
    s=SistemaActivo(p,mat,{'x+':Robin(10,20)})
    s.actualizar(p.mascara_activa(0))
    assert s.expuestas['x+'].ravel().tolist()==[True,False,True]
    r=resolver_vertidas(p,mat,{'x+':Robin(10,20)},10.,max_dt_s=1)
    assert r.T_final[0,0,0]<60
    assert r.T_final[2,0,0]<60
    comprobar_balance(r)


def test_caras_futuras_pasan_a_interfaz_en_colocacion():
    p,mat=programa([1,2],inicios={2:10.},temperaturas={1:60.,2:20.})
    s=SistemaActivo(p,mat,{'x+':Robin(10,20)})
    s.actualizar(p.mascara_activa(0))
    assert s.expuestas['x+'][0,0,0]
    assert s.G[0][0,0,0]==0
    s.actualizar(p.mascara_activa(10))
    assert not s.expuestas['x+'][0,0,0]
    assert s.G[0][0,0,0]>0


def test_no_hay_transferencia_antes_de_colocar_y_no_reinicia_lo_anterior():
    p,mat=programa([1,2],inicios={2:10.},temperaturas={1:60.,2:20.})
    r=resolver_vertidas(p,mat,{},20.,dt_s=1.,tiempos_guardado_s=[0,5,10,20],guardar_campos=True)
    assert np.isnan(r.T[1,1,0,0])
    j=int(np.where(r.t==10)[0][0])
    assert r.T[j,:,0,0]==pytest.approx([60.,20.])
    assert r.T_final[0,0,0]<60
    assert r.T_final[1,0,0]>20
    assert r.balance[j]['energia_colocacion_J']>r.balance[j-1]['energia_colocacion_J']
    comprobar_balance(r)


def test_antes_de_primera_vertida_no_hay_hormigon_ni_energia():
    p,mat=programa([1],inicios={1:10.},elevacion=45)
    r=resolver_vertidas(p,mat,{c:Robin(10,15) for c in CARAS},20.,dt_s=1,tiempos_guardado_s=[0,5,10,20],guardar_campos=True)
    assert np.isnan(r.T[0]).all() and np.isnan(r.T[1]).all()
    assert r.balance[1]['energia_J']==0
    assert r.balance[1]['energia_hidratacion_J']==0
    assert r.T[2].item()==20.
    comprobar_balance(r)


def test_hidratacion_adiabatica_uniforme_coincide_con_formula():
    p,mat=programa(np.ones((3,2,2)),temperaturas={1:18.},elevacion=45,tau=300)
    r=resolver_vertidas(p,mat,{},100.,dt_s=3.,tiempos_guardado_s=[0,50,100],guardar_campos=True)
    exacta=18+p.modelos['M1'].elevacion_adiabatica(r.t)
    assert np.allclose(r.T,exacta[:,None,None,None],atol=1e-10)
    comprobar_balance(r)


def test_fuente_en_vertidas_escalonadas_conserva_energia_total():
    p,mat=programa([1,2],inicios={2:10.},temperaturas={1:18.,2:16.},elevacion=45,tau=300)
    r=resolver_vertidas(p,mat,{},100.,dt_s=3.,tiempos_guardado_s=[0,10,50,100])
    e1=p.modelos['M1'].energia_entre_J_m3(0,100)
    e2=p.modelos['M1'].energia_entre_J_m3(0,90)
    assert r.balance[-1]['energia_hidratacion_J']==pytest.approx((e1+e2)*p.grilla.volumen_celda)
    comprobar_balance(r)


def test_neumann_positivo_aporta_q_por_area_por_tiempo():
    p,mat=programa([1,0,1])
    r=resolver_vertidas(p,mat,{'z+':Neumann(50)},10.,dt_s=1)
    assert r.balance[-1]['energia_bordes_J']==pytest.approx(50*.01*2*10)
    comprobar_balance(r)


def test_dt_excesivo_se_rechaza():
    p,mat=programa(np.ones((4,4,4)))
    with pytest.raises(ValueError,match='estable'):
        resolver_vertidas(p,mat,{c:Dirichlet(20) for c in CARAS},100.,dt_s=100,tiempos_guardado_s=[0,100])


def test_dos_materiales_conducen_con_resistencias_en_serie():
    p,mat=programa([1,2],temperaturas={1:60,2:20})
    p.vertidas[2]['mezcla']='M2'
    mat['M2']=Material('M2',.5,1000,10)
    p.modelos['M2']=HidratacionExponencial(0,1000,1000,10)
    s=SistemaActivo(p,mat,{});s.actualizar(p.mascara_activa(0))
    assert s.G[0][0,0,0]==pytest.approx(.01/(.05/1+.05/.5))
    r=resolver_vertidas(p,mat,{},10.,dt_s=1)
    comprobar_balance(r)
