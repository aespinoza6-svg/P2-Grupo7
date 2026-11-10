"""E10: clima, snapshots, bordes ambientales, moldaje y suelo. Sin internet."""
import json
from datetime import datetime, timedelta, timezone
import numpy as np
import pytest
from mcplan.clima import descargar_snapshot, cargar_snapshot, ultimo_snapshot, serie_clima, serie_constante
from mcplan.bordes import BordeAmbiental, h_conveccion, h_radiacion, caras_ambientales
from mcplan.suelo import SUELO_ID, agregar_suelo, modelo_suelo, diferencial_por_vertida
from mcplan.voxelizacion import GrillaVoxeles
from mcplan.hidratacion import HidratacionExponencial
from mcplan.construccion import ProgramaVertidas
from mcplan.solver1d import Material, Neumann, Robin
from mcplan.solver3d import CARAS
from mcplan.solver_vertidas import SistemaActivo, resolver_vertidas


def respuesta(inicio_local, horas, T, V=2.0, G=0.0, desfase=-3 * 3600):
    t0 = datetime.fromisoformat(inicio_local)
    tiempos = [(t0 + timedelta(hours=i)).strftime('%Y-%m-%dT%H:%M') for i in range(horas)]
    f = (lambda x, i: x(i)) if callable(T) else (lambda x, i: x)
    return dict(utc_offset_seconds=desfase, hourly=dict(time=tiempos, temperature_2m=[f(T, i) for i in range(horas)],
                                                        wind_speed_10m=[V] * horas, shortwave_radiation=[G] * horas))


def pedir_falso(url, parametros):
    """Imita Open-Meteo: 3 dias de pronostico a 20 C; historico a 10 + ano - 2024 C."""
    if 'archive' in url:
        ano = int(parametros['start_date'][:4])
        return url, respuesta(parametros['start_date'] + 'T00:00', 24 * 6, 10.0 + (ano - 2024))
    return url, respuesta('2026-11-10T00:00', 72, 20.0)


@pytest.fixture
def snapshot(tmp_path):
    carpeta = descargar_snapshot(tmp_path, -37.6921, -72.624, dias=3, pedir=pedir_falso,
                                 historico=dict(desde='11-10', dias=5, anos=[2024, 2025]),
                                 ahora=datetime(2026, 11, 10, 15, tzinfo=timezone.utc))
    return carpeta


def test_snapshot_guarda_metadata_minima_y_no_se_sobrescribe(snapshot, tmp_path):
    meta = json.loads((snapshot / 'metadata.json').read_text(encoding='utf-8'))
    for campo in ('source', 'model', 'retrieved_at', 'latitude', 'longitude', 'timezone', 'variables', 'units'):
        assert campo in meta
    assert meta['units']['temperature_2m'] == 'degC' and snapshot.name == '2026-11-10T1500'
    assert len(meta['requests']) == 3
    with pytest.raises(FileExistsError):
        descargar_snapshot(tmp_path, -37.6921, -72.624, dias=3, pedir=pedir_falso,
                           ahora=datetime(2026, 11, 10, 15, tzinfo=timezone.utc))
    assert ultimo_snapshot(tmp_path) == snapshot


def test_la_cola_que_no_cubre_el_pronostico_se_completa_con_historico(snapshot):
    snap = cargar_snapshot(snapshot)
    s = serie_clima(snap, '2026-11-11T00:00:00-03:00', 96)
    # El pronostico llega hasta las 23:00 local del 12-nov: 48 h desde el inicio, inclusive.
    assert s.horas_pronostico == 48 and s.horas_relleno == 49
    assert np.all(s.T_aire_C[s.fuente == 'pronostico'] == 20.0)
    assert np.allclose(s.T_aire_C[s.fuente == 'historico'], 10.5)  # promedio de 10 y 11
    frio = serie_clima(snap, '2026-11-11T00:00:00-03:00', 96, 'historico_frio')
    assert np.allclose(frio.T_aire_C[frio.fuente == 'historico'], 10.0) and frio.metodo_cola == 'historico_frio:2024'


def test_sin_historico_y_sin_cobertura_avisa(tmp_path):
    carpeta = descargar_snapshot(tmp_path, -37.7, -72.6, dias=3, pedir=pedir_falso,
                                 ahora=datetime(2026, 11, 10, 15, tzinfo=timezone.utc))
    with pytest.raises(ValueError, match='no cubre'):
        serie_clima(cargar_snapshot(carpeta), '2026-11-11T00:00:00-03:00', 96)


def test_coeficientes_de_viento_y_radiacion():
    assert h_conveccion(0) == pytest.approx(5.7) and h_conveccion(5) == pytest.approx(24.7)
    assert h_radiacion(15.0) == pytest.approx(4.88, abs=0.01)
    assert h_conveccion(0) + h_radiacion(15.0) == pytest.approx(10.6, abs=0.05)


def programa(ids, inicios, temperaturas=20., d=(.1, .1, .1), elevacion=0., suelo=None):
    ids = np.asarray(ids, dtype=np.int32)
    props = {int(v): dict(Inicio_h=inicios.get(int(v), 0.), TemperaturaColocacion_C=temperaturas, TipoHormigon='M1',
                          RetiroMoldaje_h=1.0) for v in set(np.unique(ids)) - {0}}
    g = GrillaVoxeles(np.zeros(3), d, ids, props)
    m = Material('M1', 2., 2400., 900.)
    modelos, materiales = {'M1': HidratacionExponencial(elevacion, 3600., m.rho, m.cp)}, {'M1': m}
    if suelo:
        g = agregar_suelo(g, **suelo)
        modelos['SUELO'] = modelo_suelo(1900., 1000.)
        materiales['SUELO'] = Material('SUELO', 1.5, 1900., 1000.)
    return ProgramaVertidas(g, modelos), materiales


def test_borde_ambiental_constante_reproduce_robin():
    p, mat = programa(np.ones((6, 4, 3)), {}, elevacion=20.)
    clima = serie_constante(15.0, 10)
    a = resolver_vertidas(p, mat, {c: Robin(10., 15.) for c in CARAS}, 10 * 3600, max_dt_s=600.)
    b = resolver_vertidas(p, mat, {c: BordeAmbiental(clima, h_fijo=10.) for c in CARAS}, 10 * 3600, max_dt_s=600.)
    np.testing.assert_allclose(a.T_final, b.T_final, rtol=0, atol=1e-12)


def test_el_sol_entra_como_temperatura_sol_aire_sin_doble_conteo():
    p, mat = programa(np.ones((2, 2, 2)), {})
    s = SistemaActivo(p, mat, {c: Neumann(0.) for c in CARAS}); s.actualizar(np.ones((2, 2, 2), bool))
    clima = serie_constante(15.0, 2, viento_m_s=0., ghi_W_m2=800.)
    b = BordeAmbiental(clima, sol=True, absortividad=.6, h_fijo=10.)
    G, Tref, Q = b.coeficientes(0., 'z+', s.expuestas['z+'], s)
    assert np.allclose(Tref, 15 + .6 * 800 / 10) and np.all(Q == 0)
    # Con la celda a T_aire, el flujo entrante es alpha*G atenuado solo por la media celda.
    A, dz, k = s.areas[2], .1, 2.
    assert G[0] * (Tref[0] - 15.) == pytest.approx(A * .6 * 800 / (1 + 10 * dz / (2 * k)))


def test_desmoldar_quita_la_resistencia_en_la_hora_de_retiro_de_cada_vertida():
    p, mat = programa([[[1, 2]]], {1: 0., 2: 5 / 3600 * 3600})
    s = SistemaActivo(p, mat, {}); s.actualizar(np.ones((1, 1, 2), bool))
    b = BordeAmbiental(serie_constante(10., 20), h_fijo=10., resistencias={1: (1.43, 1.0), 2: (1.43, 1.0)})
    antes, _, _ = b.coeficientes(1800., 'x-', s.expuestas['x-'], s)
    despues, _, _ = b.coeficientes(5400., 'x-', s.expuestas['x-'], s)
    assert despues[0] > 5 * antes[0]          # V01 ya se desmoldo (1 h desde su inicio)
    assert despues[1] == pytest.approx(antes[1])  # V02 entro despues: aun tiene moldaje
    assert b.G_maximo('x-', s.expuestas['x-'], s)[0] == pytest.approx(despues[0])


def test_suelo_se_agrega_bajo_la_fundacion_y_cierra_el_balance():
    p, mat = programa(np.ones((4, 4, 2)), {}, elevacion=30., d=(.2, .2, .2),
                      suelo=dict(espesor_m=.6, margen_m=.4, T_inicial_C=12.))
    g = p.grilla
    assert g.vertida_id.shape == (8, 8, 5) and np.count_nonzero(g.vertida_id == SUELO_ID) == 8 * 8 * 3
    assert np.all(g.vertida_id[2:6, 2:6, 3:] == 1) and np.all(g.vertida_id[:2, :, 3:] == 0)
    clima = serie_constante(12.0, 24)
    bordes = caras_ambientales(clima, suelo=dict(id=SUELO_ID, T_fondo_C=12.), h_fijo=10., absortividad=0.)
    vistos = []
    r = resolver_vertidas(p, mat, bordes, 24 * 3600, max_dt_s=600., tiempos_guardado_s=[0, 12 * 3600, 24 * 3600],
                          al_guardar=lambda t, T, a, s: vistos.append(diferencial_por_vertida(T, a, s)))
    for f in r.balance:
        escala = max(1, abs(f['energia_colocacion_J']) + abs(f['energia_hidratacion_J']) + abs(f['energia_bordes_J']))
        assert abs(f['residuo_balance_J']) / escala < 1e-10
    suelo = g.vertida_id == SUELO_ID
    assert r.T_final[suelo].max() > 12.2          # el hormigon caliente calienta el suelo
    assert len(vistos) == 3 and vistos[-1][1] == 1 and vistos[-1][0] > 0


def test_costados_del_suelo_son_adiabaticos_y_el_fondo_a_temperatura_fija():
    p, mat = programa(np.ones((2, 2, 1)), {}, d=(.2, .2, .2), suelo=dict(espesor_m=.4, margen_m=.2, T_inicial_C=12.))
    s = SistemaActivo(p, mat, caras_ambientales(serie_constante(12., 2), suelo=dict(id=SUELO_ID, T_fondo_C=12.), h_fijo=10.))
    s.actualizar(np.ones(p.grilla.vertida_id.shape, bool) & (p.grilla.vertida_id > 0))
    G, _, _ = s.bordes['x-'].coeficientes(0., 'x-', s.expuestas['x-'], s)
    ids = p.grilla.vertida_id[s.expuestas['x-']]
    assert np.all(G[ids == SUELO_ID] == 0) and np.all(G[ids == 1] > 0)
    Gf, Tf, _ = s.bordes['z-'].coeficientes(0., 'z-', s.expuestas['z-'], s)
    assert np.all(Gf > 0) and np.all(Tf == 12.)
