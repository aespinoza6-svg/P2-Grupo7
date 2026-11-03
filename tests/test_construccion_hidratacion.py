"""E07: edades por vertida, activacion y energia de hidratacion."""
from pathlib import Path
import numpy as np
import pytest
from mcplan import bim
from mcplan.hidratacion import HidratacionExponencial,modelos_desde_configuracion
from mcplan.construccion import ProgramaVertidas,EstadoConstruccion
from mcplan.voxelizacion import GrillaVoxeles

RAIZ=Path(__file__).resolve().parents[1]
HORA=3600.


@pytest.fixture
def modelo():
    return HidratacionExponencial(45.,30*HORA,2400.,900.)


@pytest.fixture
def programa(modelo):
    ids=np.array([1,2,0],dtype=np.int32).reshape(3,1,1)
    props={1:{'Inicio_h':0,'TemperaturaColocacion_C':18,'TipoHormigon':'M1'},
           2:{'Inicio_h':6,'TemperaturaColocacion_C':16,'TipoHormigon':'M1'}}
    g=GrillaVoxeles(np.zeros(3),(1.,1.,1.),ids,props)
    return ProgramaVertidas(g,{'M1':modelo})


def test_curva_y_fuente_con_unidades_y_parametros_de_la_guia(modelo):
    assert modelo.elevacion_adiabatica(0)==0
    assert modelo.fuente_W_m3(0)==pytest.approx(900.)
    assert modelo.elevacion_adiabatica(modelo.tau_s)==pytest.approx(45*(1-np.exp(-1)))
    assert modelo.elevacion_adiabatica(100*modelo.tau_s)==pytest.approx(45)
    assert modelo.energia_entre_J_m3(0,100*modelo.tau_s)==pytest.approx(2400*900*45)


def test_no_hay_hidratacion_antes_de_colocar(modelo):
    assert modelo.fuente_W_m3(-100)==0
    assert modelo.elevacion_adiabatica(-100)==0
    assert modelo.energia_entre_J_m3(-100,-10)==0
    assert modelo.energia_entre_J_m3(-100,100)==pytest.approx(modelo.energia_entre_J_m3(0,100))


def test_fuente_es_derivada_energetica_de_la_curva(modelo):
    a=2*HORA;eps=.01
    derivada=(modelo.elevacion_adiabatica(a+eps)-modelo.elevacion_adiabatica(a-eps))/(2*eps)
    assert modelo.fuente_W_m3(a)==pytest.approx(modelo.capacidad_volumetrica*derivada,rel=1e-7)


@pytest.mark.parametrize('parametros',[(45,0,2400,900),(-1,100,2400,900),(45,100,-1,900),(45,100,2400,np.nan)])
def test_parametros_no_fisicos_se_rechazan(parametros):
    with pytest.raises(ValueError,match='Parametros'):
        HidratacionExponencial(*parametros)


def test_modelo_no_admite_intervalos_invertidos(modelo):
    with pytest.raises(ValueError,match='intervalo'):
        modelo.energia_entre_J_m3(10,5)


def test_activacion_justo_en_horario_y_exterior_nunca_activo(programa):
    assert programa.mascara_activa(0).ravel().tolist()==[True,False,False]
    assert programa.mascara_activa(6*HORA-1).ravel().tolist()==[True,False,False]
    assert programa.mascara_activa(6*HORA).ravel().tolist()==[True,True,False]
    assert np.isnan(programa.edad_s(0).ravel()[1:]).all()
    assert programa.edad_s(6*HORA).ravel()[:2]==pytest.approx([6*HORA,0])


def test_cada_vertida_tiene_su_propio_reloj(programa,modelo):
    q=programa.fuente_W_m3(6*HORA).ravel()
    assert q[0]==pytest.approx(modelo.fuente_W_m3(6*HORA))
    assert q[1]==pytest.approx(900.)
    assert q[1]>q[0]
    assert q[2]==0
    assert programa.fuente_W_m3(5*HORA).ravel()[1]==0


def test_celdas_nuevas_reciben_temperatura_sin_reiniciar_anteriores(programa):
    e=EstadoConstruccion(programa)
    assert np.isnan(e.temperatura_C).all()
    e.activar_hasta(0)
    assert e.temperatura_C.ravel()[0]==18
    assert np.isnan(e.temperatura_C.ravel()[1:]).all()
    e.temperatura_C[0,0,0]=27. # simula una temperatura que el futuro solver ya cambio
    nuevas=e.activar_hasta(6*HORA)
    assert nuevas.ravel().tolist()==[False,True,False]
    assert e.temperatura_C.ravel()[:2]==pytest.approx([27.,16.])
    assert not e.activar_hasta(6*HORA).any()
    assert e.temperatura_C[0,0,0]==27
    assert np.isnan(e.temperatura_C[2,0,0])


def test_estado_no_puede_retroceder(programa):
    e=EstadoConstruccion(programa);e.activar_hasta(10)
    with pytest.raises(ValueError,match='retroceder'):e.activar_hasta(9)


def test_intervalo_se_parte_en_eventos_exactos(programa):
    assert programa.partir_intervalo(5*HORA,7*HORA)==pytest.approx([5*HORA,6*HORA,7*HORA])
    assert programa.partir_intervalo(6*HORA,7*HORA)==pytest.approx([6*HORA,7*HORA])


def test_energia_no_incluye_tiempo_previo_a_colocacion(programa,modelo):
    e=programa.energia_generada_entre_J_m3(5*HORA,7*HORA).ravel()
    assert e[0]==pytest.approx(modelo.energia_entre_J_m3(5*HORA,7*HORA))
    assert e[1]==pytest.approx(modelo.energia_entre_J_m3(0,HORA))
    assert e[2]==0


def test_integracion_adiabatica_independiente_del_paso_y_balance(programa,modelo):
    finales=[]
    for paso in [HORA,7*HORA,24*HORA]:
        acumulado=np.zeros(programa.grilla.vertida_id.shape)
        for t0 in np.arange(0,24*HORA,paso):
            acumulado+=programa.energia_generada_entre_J_m3(t0,min(t0+paso,24*HORA))
        esperada=np.array([modelo.energia_entre_J_m3(0,24*HORA),modelo.energia_entre_J_m3(0,18*HORA),0]).reshape(3,1,1)
        assert acumulado==pytest.approx(esperada,rel=1e-12,abs=1e-8)
        finales.append(acumulado.copy())
    assert finales[0]==pytest.approx(finales[-1],rel=1e-12)


def test_configuracion_se_declara_pedagogica():
    caso=bim.cargar_yaml(RAIZ/'config'/'caso_base.yaml')
    config=bim.cargar_yaml(RAIZ/'config'/'hidratacion_e07.yaml')
    assert config['estado']=='pedagogico_no_calibrado'
    modelos=modelos_desde_configuracion(caso,config)
    assert modelos['M1'].fuente_W_m3(0)==pytest.approx(900.)


def test_npz_cargado_conserva_horarios_para_activacion(programa,tmp_path):
    ruta=tmp_path/'g.npz';programa.grilla.guardar(ruta)
    g=GrillaVoxeles.cargar(ruta)
    p=ProgramaVertidas(g,programa.modelos)
    assert p.eventos_s==pytest.approx([0,6*HORA])
    assert np.array_equal(g.vertida_id,programa.grilla.vertida_id)


def test_mezcla_sin_modelo_se_rechaza(programa):
    with pytest.raises(ValueError,match='mezcla'):
        ProgramaVertidas(programa.grilla,{})
