"""E13: curvas tabuladas, S25 a 25 mm, cronograma, trayectorias y enumeracion ponderada."""
from pathlib import Path
import numpy as np
import pytest
from mcplan.hidratacion import HidratacionTabulada, leer_curvas_tabuladas
from mcplan.riesgo import cronograma, temperatura_25mm, probabilidades, elegir_trayectorias
from mcplan.clima import serie_constante, snapshot_de_miembro, miembros_de
from mcplan.bordes import caras_ambientales
from mcplan.construccion import ProgramaVertidas
from mcplan.solver1d import Material
from mcplan.solver_vertidas import SistemaActivo
from mcplan.voxelizacion import GrillaVoxeles

RAIZ = Path(__file__).resolve().parents[1]
CURVAS = RAIZ / 'data' / 'hidratacion' / 'source_terms.txt'


def test_la_energia_de_una_curva_tabulada_es_exacta_y_no_inventa_cola():
    c = HidratacionTabulada([0, 1, 2], [0, 3600, 0], 2400, 900)     # triangulo: 0,5 * 7200 s * 3600 W/m3
    assert c.energia_total_J_m3 == pytest.approx(12.96e6)
    assert c.energia_entre_J_m3(1800, 5400) == pytest.approx(12.96e6 - 2 * 1.62e6)
    pasos = np.linspace(0, 7200, 37)
    assert np.sum(c.energia_entre_J_m3(pasos[:-1], pasos[1:])) == pytest.approx(12.96e6, rel=1e-12)
    assert c.energia_entre_J_m3(-3600, 0) == 0 and c.energia_entre_J_m3(7200, 1e6) == 0   # antes y despues: nada
    assert c.fuente_W_m3(3600) == 3600 and c.fuente_W_m3(9000) == 0
    assert c.elevacion_adiabatica(7200) == pytest.approx(12.96e6 / (2400 * 900))


def test_las_ocho_curvas_del_profesor_reproducen_el_calor_total_del_encabezado():
    curvas = leer_curvas_tabuladas(CURVAS)
    assert list(curvas) == ['M1_SS', 'M2_SS', 'M3_SS', 'M4_SS', 'M1_22A', 'M2_22A', 'M3_22A', 'M4_22A']
    encabezado = dict(M1_SS=(88.1, 34.8), M2_SS=(95.7, 37.8), M3_SS=(94.1, 37.2), M4_SS=(100.2, 39.5),
                      M1_22A=(61.2, 24.2), M2_22A=(82.4, 32.5), M3_22A=(73.3, 28.9), M4_22A=(73.6, 29.1))
    for n, (MJ, K) in encabezado.items():
        c = curvas[n]
        assert (c.rho_kg_m3, c.cp_J_kgK) == (2476.0, 1023.0)
        assert c.energia_total_J_m3 / 1e6 == pytest.approx(MJ, abs=0.06)
        assert c.elevacion_adiabatica(1e9) == pytest.approx(K, abs=0.06)


def test_cronograma_cuenta_la_duracion_de_cada_vertida_y_las_esperas():
    vol = {1: 80.0, 2: 40.0, 3: 120.0}
    filas, D = cronograma(vol, 40.0, espera_h=6)
    assert [filas[v]['duracion_h'] for v in vol] == [2, 1, 3]
    assert [filas[v]['inicio_h'] for v in vol] == [0, 8, 15] and D == 18
    assert filas[3]['activacion_h'] == 16.5                    # mitad de su colocacion
    unica, Du = cronograma(vol, 40.0, unica=True)
    assert Du == 6 and all(f['activacion_h'] == 3 for f in unica.values())


def test_S25_reconstruye_la_temperatura_a_25_mm_con_el_flujo_de_la_cara():
    """Bloque a 30 C, aire a 10 C, h = 10. Cara lateral (dx = 0,2 m, k = 2):
    q = 20 / (0,05 + 0,1) = 133,3 W/m2 y, con perfil lineal, T25 = 30 - q (0,1 - 0,025) / 2 = 25 C.
    En la cara superior dz/2 = 25 mm: T25 es el centro de la celda (30 C)."""
    ids = np.ones((3, 3, 2), dtype=np.int32)
    g = GrillaVoxeles(np.zeros(3), (0.2, 0.2, 0.05), ids, {1: dict(Inicio_h=0, TemperaturaColocacion_C=30, TipoHormigon='M1')})
    m = Material('M1', 2.0, 2400., 900.)
    p = ProgramaVertidas(g, {'M1': HidratacionTabulada([0, 1], [0, 0], 2400., 900.)})
    s = SistemaActivo(p, {'M1': m}, caras_ambientales(serie_constante(10.0, 2), h_fijo=10.0, absortividad=0.0))
    s.actualizar(ids > 0)
    T25, i, cara = temperatura_25mm(np.full(ids.shape, 30.0), s, 0.0)
    assert T25 == pytest.approx(25.0) and cara in ('x-', 'x+', 'y-', 'y+')


def test_la_probabilidad_conjunta_es_la_union_no_la_suma_de_marginales():
    filas = [dict(curva='a', clima='x', MT=75, MD=25), dict(curva='a', clima='y', MT=60, MD=25),
             dict(curva='b', clima='x', MT=60, MD=10), dict(curva='b', clima='y', MT=75, MD=10)]
    P = probabilidades(filas, dict(a=.5, b=.5), dict(x=.5, y=.5))
    assert P['T'] == .5 and P['dT'] == .5 and P['conjunta'] == .75
    with pytest.raises(ValueError):
        probabilidades(filas[:3], dict(a=.5, b=.5), dict(x=.5, y=.5))     # falta un escenario


def test_trayectorias_parten_el_ensemble_en_grupos_iguales_y_los_pesos_suman_uno():
    series = {f'member{i:02d}': serie_constante(10.0 + i, 24) for i in range(9)}
    for s in series.values():
        s.fuente = np.full(len(s.t_s), 'pronostico')
    elegidas, extra = elegir_trayectorias(series, 3)
    assert [e['miembro'] for e in elegidas] == ['member01', 'member04', 'member07']
    assert sum(e['peso'] for e in elegidas) == pytest.approx(1.0)
    assert extra['frio_extremo'] == 'member00' and extra['calido_extremo'] == 'member08'


def test_un_miembro_del_ensemble_se_vuelve_un_snapshot_normal():
    h = dict(time=['2026-11-19T00:00'], temperature_2m=[10], wind_speed_10m=[1], shortwave_radiation=[0],
             temperature_2m_member01=[12], wind_speed_10m_member01=[2], shortwave_radiation_member01=[5])
    snap = dict(nombre='ens', pronostico=dict(hourly=h), historico={})
    assert miembros_de(snap['pronostico']) == ['member00', 'member01']
    s1 = snapshot_de_miembro(snap, 'member01')
    assert s1['pronostico']['hourly']['temperature_2m'] == [12] and s1['nombre'] == 'ens:member01'
    assert snapshot_de_miembro(snap, 'member00')['pronostico']['hourly']['wind_speed_10m'] == [1]
