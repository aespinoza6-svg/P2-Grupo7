"""E11: exportacion VTK de campos 3D, superficie y lectura con PyVista."""
import importlib.util
import subprocess
import sys
from pathlib import Path
import numpy as np
import pytest
from mcplan.exportar import celdas_superficie, escribir_vtk, leer_vtk, campo_vtk
from mcplan.voxelizacion import GrillaVoxeles

RAIZ = Path(__file__).resolve().parents[1]


def grilla_chica():
    ids = np.zeros((5, 4, 3), dtype=np.int32)
    ids[:, :, 0] = 90            # suelo
    ids[1:4, 1:3, 1:] = 1        # hormigon
    return GrillaVoxeles(np.array([-1.0, -0.8, -0.2]), (0.4, 0.4, 0.2), ids, {})


def test_la_superficie_es_hormigon_con_cara_al_aire_y_no_hacia_abajo():
    g = grilla_chica()
    activa = g.vertida_id > 0
    sup = celdas_superficie(activa, g.vertida_id)
    assert not sup[g.vertida_id == 90].any()            # el suelo no es superficie del hormigon
    assert sup[g.vertida_id == 1].all()                 # todo el bloque chico toca el aire
    activa2 = activa.copy(); activa2[2, 1, 2] = False   # una celda sin colocar deja expuesta a su vecina
    assert celdas_superficie(activa2, g.vertida_id)[2, 1, 1]


def test_el_vtk_guarda_temperatura_y_etiquetas_en_orden_vtk(tmp_path):
    g = grilla_chica()
    activa = g.vertida_id > 0
    T = np.arange(g.vertida_id.size, dtype=float).reshape(g.vertida_id.shape)
    ruta = campo_vtk(tmp_path / 'c.vtk', g, T, activa, 12.0, 'prueba')
    texto = ruta.read_bytes()[:400].decode('ascii', 'replace')
    assert 'DIMENSIONS 6 5 4' in texto and 'CELL_DATA 60' in texto and 'BINARY' in texto
    d = leer_vtk(ruta)
    assert d['forma'] == [5, 4, 3] and d['paso'] == pytest.approx([0.4, 0.4, 0.2])
    np.testing.assert_array_equal(d['campos']['vertida_id'], g.vertida_id)
    np.testing.assert_allclose(d['campos']['temperatura_C'][activa], T[activa])
    assert np.isnan(d['campos']['temperatura_C'][~activa]).all()
    assert set(d['campos']) == {'temperatura_C', 'vertida_id', 'activa', 'superficie'}


def test_campos_de_formas_distintas_se_rechazan(tmp_path):
    with pytest.raises(ValueError):
        escribir_vtk(tmp_path / 'x.vtk', (0, 0, 0), (1, 1, 1), dict(a=np.zeros((2, 2, 2)), b=np.zeros((2, 2, 3))))


def test_pyvista_lee_el_mismo_campo(tmp_path):
    pv = pytest.importorskip('pyvista')
    g = grilla_chica()
    activa = g.vertida_id > 0
    T = np.random.default_rng(1).uniform(10, 50, g.vertida_id.shape)
    ruta = campo_vtk(tmp_path / 'c.vtk', g, T, activa, 24.0)
    m = pv.read(str(ruta))
    assert m.dimensions == (6, 5, 4) and m.n_cells == 60
    np.testing.assert_allclose(m.origin, g.origen)
    leido = np.asarray(m.cell_data['temperatura_C']).reshape(g.vertida_id.shape, order='F')
    np.testing.assert_allclose(leido[activa], T[activa], rtol=1e-6)


def test_el_visor_pyvista_encuentra_nucleo_y_superficie(tmp_path):
    pytest.importorskip('pyvista')
    g = grilla_chica()
    activa = g.vertida_id > 0
    T = np.full(g.vertida_id.shape, 12.0)
    T[g.vertida_id == 1] = 30.0
    T[2, 1, 1] = 45.0      # nucleo
    T[1, 1, 2] = 20.0      # superficie mas fria
    ruta = campo_vtk(tmp_path / 'c.vtk', g, T, activa, 24.0)
    png = tmp_path / 'v.png'
    r = subprocess.run([sys.executable, '-X', 'utf8', str(RAIZ / 'scripts' / '15_ver_vtk.py'), str(ruta), '--plano', 'y',
                        '--pos', '-0.2', '--png', str(png)], capture_output=True, text=True, encoding='utf-8', cwd=RAIZ)
    assert r.returncode == 0, r.stderr
    assert 'Nucleo 45.00 C' in r.stdout and 'superficie 20.00 C' in r.stdout and png.is_file()
