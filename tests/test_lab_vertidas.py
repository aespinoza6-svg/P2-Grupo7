"""Lab: la secuencia sale del IFC y el recorrido completo funciona en una malla chica."""
import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pytest
from mcplan.voxelizacion import GrillaVoxeles

RAIZ = Path(__file__).resolve().parents[1]
RUTA = RAIZ / 'scripts' / '11_lab_vertidas.py'
IFC = RAIZ / 'fundaciones' / 'fundacion.ifc'
espec = importlib.util.spec_from_file_location('lab_vertidas', RUTA)
lab = importlib.util.module_from_spec(espec)
espec.loader.exec_module(lab)


def correr(*argumentos):
    return subprocess.run([sys.executable, '-X', 'utf8', str(RUTA), *map(str, argumentos)],
                          capture_output=True, text=True, encoding='utf-8', cwd=RAIZ)


def grilla_dos_pisos(inicio_abajo, inicio_arriba):
    ids = np.array([[[1, 2]]], dtype=np.int32)
    props = {1: dict(Inicio_h=inicio_abajo, TemperaturaColocacion_C=18.), 2: dict(Inicio_h=inicio_arriba, TemperaturaColocacion_C=18.)}
    return GrillaVoxeles(np.zeros(3), (.1, .1, .1), ids, props)


def test_lee_cambios_validos():
    assert lab.leer_cambio('V02=24') == (2, 24.)
    assert lab.leer_cambio('v2=-3.5') == (2, -3.5)


@pytest.mark.parametrize('texto', ['V02', 'V02=hola', 'hola=3', 'V02=nan', 'V02=inf'])
def test_rechaza_cambios_invalidos(texto):
    with pytest.raises(argparse.ArgumentTypeError):
        lab.leer_cambio(texto)


def test_avisa_si_una_vertida_entra_antes_que_la_de_abajo():
    assert lab.avisos_programa(grilla_dos_pisos(0., 6.)) == []
    assert lab.avisos_programa(grilla_dos_pisos(0., 0.)) == []
    avisos = lab.avisos_programa(grilla_dos_pisos(8., 6.))
    assert len(avisos) == 1 and 'V02' in avisos[0] and 'sin apoyo' in avisos[0]


def test_avisa_si_una_vertida_queda_fuera_del_horizonte():
    avisos = lab.avisos_programa(grilla_dos_pisos(0., 100.), horizonte=72.)
    assert len(avisos) == 1 and 'horizonte' in avisos[0]


@pytest.mark.parametrize('argumentos', [['--inicio', 'V09=3'], ['--ocultar', 'V09'], ['--inicio', 'V02=-1'],
                                        ['--solo-voxelizar', '--temp', 'V02=10'], ['--dx', '0']])
def test_pedidos_imposibles_terminan_con_error_claro(argumentos, tmp_path):
    r = correr(*argumentos, '--salida', tmp_path)
    assert r.returncode == 2 and 'error:' in r.stderr
    assert not any(tmp_path.iterdir())


def test_recorrido_completo_cambia_solo_por_el_ifc(tmp_path):
    antes = IFC.read_bytes()
    r = correr('--dx', .8, '--dy', .8, '--horizonte-h', 6, '--inicio', 'V02=3', '--temp', 'V02=10',
               '--ocultar', 'V05', '--salida', tmp_path)
    assert r.returncode == 0, r.stderr
    assert IFC.read_bytes() == antes

    malla = tmp_path / 'dx0.8_dy0.8_dz0.05'
    escenario = malla / 'V02_inicio3__V02_temp10'
    assert (malla / 'vista_3d_sin_V05.png').is_file()
    assert (escenario / 'comparacion_6h.png').is_file()

    # El solver recibio el dato que se cambio en el IFC, sin tocar las demas vertidas.
    base = GrillaVoxeles.cargar(malla / 'base' / 'voxeles' / 'grilla_vertidas.npz')
    cambiada = GrillaVoxeles.cargar(escenario / 'voxeles' / 'grilla_vertidas.npz')
    assert base.vertidas[2]['Inicio_h'] == 6.0
    assert cambiada.vertidas[2]['Inicio_h'] == 3.0
    assert cambiada.vertidas[2]['TemperaturaColocacion_C'] == 10.0
    assert cambiada.vertidas[3]['Inicio_h'] == base.vertidas[3]['Inicio_h']
    np.testing.assert_array_equal(base.vertida_id, cambiada.vertida_id)

    # V02 existe en el escenario modificado a las 3 h y parte en 10 °C.
    historias = lab.leer_historias(escenario / 'simulacion_6h')
    fila = int(np.argmin(abs(historias[2]['t'] - 3)))
    assert historias[2]['media'][fila] == pytest.approx(10.0)
    assert np.isnan(historias[2]['media'][fila - 1])
    assert np.isnan(lab.leer_historias(malla / 'base' / 'simulacion_6h')[2]['media'][fila])
    resumen = json.loads((escenario / 'simulacion_6h' / 'resumen.json').read_text(encoding='utf-8'))
    assert resumen['residuo_balance_relativo_max'] < 1e-8

    # La segunda corrida reutiliza la base, pero no un resultado de otras entradas.
    r2 = correr('--dx', .8, '--dy', .8, '--horizonte-h', 6, '--temp', 'V02=12', '--salida', tmp_path)
    assert r2.returncode == 0, r2.stderr
    assert 'Simulacion ya calculada' in r2.stdout
    assert (malla / 'V02_temp12' / 'comparacion_6h.png').is_file()
