"""visor_laboratorio.html: misma geometria que el IFC y mismos resultados que el solver Python."""
import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pytest
from mcplan import bim
from mcplan.construccion import ProgramaVertidas
from mcplan.hidratacion import modelos_desde_configuracion
from mcplan.solver1d import Dirichlet, Material, Neumann, Robin
from mcplan.solver3d import CARAS
from mcplan.solver_vertidas import resolver_vertidas
from mcplan.voxelizacion import extraer_mallas, voxelizar_mallas

RAIZ = Path(__file__).resolve().parents[1]
IFC = RAIZ / 'fundaciones' / 'fundacion.ifc'
TEXTO = (RAIZ / 'visor_laboratorio.html').read_text(encoding='utf-8')
espec = importlib.util.spec_from_file_location('visor_util', RAIZ / 'scripts' / 'actualizar_visor_laboratorio.py')
visor12 = importlib.util.module_from_spec(espec)
espec.loader.exec_module(visor12)

# Corre el motor del HTML (el Web Worker) en Node, sin navegador.
ARNES = r"""
const fs=require('fs');global.self={};global.postMessage=()=>{};
eval(fs.readFileSync(process.argv[2],'utf8')+';global.voxelizar=voxelizar;global.simular=simular;');
const G=JSON.parse(fs.readFileSync(process.argv[3],'utf8')),P=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
const v=voxelizar(G,P.d),r=simular(v,P,()=>{}),s=r.serie,i=s.t.length-1;
const ifc=v.tabla.reduce((a,x)=>a+x.ifc,0),vox=v.tabla.reduce((a,x)=>a+x.vox,0);
const media={};for(const id in s.pv)media[id]=s.pv[id].media[i];
console.log(JSON.stringify({error:100*(vox-ifc)/ifc,celdas:r.ncon,maximo:r.maximo,pasos:r.pasos,bordes:s.Ebor[i],media}));
"""


def test_la_geometria_del_visor_es_la_del_ifc():
    copia = visor12.geometria_del_visor(TEXTO)
    assert copia == json.loads(json.dumps(visor12.geometria_desde_ifc(IFC)))
    reales = {m.vertida_id: m for m in extraer_mallas(IFC)}
    assert {g['id'] for g in copia} == set(reales)
    for g in copia:
        assert len(g['t']) == 3 * len(reales[g['id']].triangulos)
        assert g['vol'] == pytest.approx(reales[g['id']].volumen_ifc_m3)


REPO = dict(arriba=dict(tipo='robin', h=10, R=0, T=15), lados=dict(tipo='robin', h=10, R=0, T=15), base=dict(tipo='neumann', q=0))
MEZCLA = dict(arriba=dict(tipo='dirichlet', T=20), lados=dict(tipo='neumann', q=50), base=dict(tipo='robin', h=5, R=0.1, T=12))


def borde_python(b):
    if b['tipo'] == 'robin':
        return Robin(b['h'], b['T'], R_extra=b['R'])
    if b['tipo'] == 'dirichlet':
        return Dirichlet(b['T'])
    return Neumann(b['q'])


@pytest.mark.skipif(shutil.which('node') is None, reason='Node no esta instalado')
@pytest.mark.parametrize('caras', [REPO, MEZCLA], ids=['bordes_del_repo', 'dirichlet_neumann_robin'])
def test_el_motor_del_visor_da_lo_mismo_que_python(tmp_path, caras):
    horas, paso = 12.0, (.4, .4, .05)
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    hc = bim.cargar_yaml(RAIZ / 'config' / 'hidratacion_e07.yaml')
    h = caso['hormigon']
    grilla = voxelizar_mallas(extraer_mallas(IFC), paso, max_celdas=4_000_000)
    programa = ProgramaVertidas(grilla, modelos_desde_configuracion(caso, hc))
    materiales = {'M1': Material('M1', h['conductividad_W_mK'], h['densidad_kg_m3'], h['calor_especifico_J_kgK'])}
    grupo = {'x-': 'lados', 'x+': 'lados', 'y-': 'lados', 'y+': 'lados', 'z-': 'base', 'z+': 'arriba'}
    bordes = {c: borde_python(caras[grupo[c]]) for c in CARAS}
    r = resolver_vertidas(programa, materiales, bordes, horas * 3600, max_dt_s=600.0,
                          tiempos_guardado_s=np.arange(0, horas + 1) * 3600)
    filas = grilla.resumen_volumenes()
    error_py = 100 * (sum(f['volumen_voxel_m3'] for f in filas) / sum(f['volumen_ifc_m3'] for f in filas) - 1)
    media_py = {f['vertida_id']: f['T_media_C'] for f in r.por_vertida if f['tiempo_h'] == horas}

    worker = re.search(r'id="worker-src">(.*?)</script>', TEXTO, re.S).group(1)
    (tmp_path / 'worker.js').write_text(worker, encoding='utf-8')
    (tmp_path / 'arnes.js').write_text(ARNES, encoding='utf-8')
    (tmp_path / 'geo.json').write_text(json.dumps(visor12.geometria_del_visor(TEXTO)), encoding='utf-8')
    P = dict(d=list(paso), k=h['conductividad_W_mK'], rho=h['densidad_kg_m3'], cp=h['calor_especifico_J_kgK'],
             theta=45, tau=30, horizonte=horas, maxdt=600, intervalo=1, **caras,
             vertidas=[dict(id=v, inicio=p['Inicio_h'], temp=p['TemperaturaColocacion_C'], calc=True) for v, p in grilla.vertidas.items()])
    (tmp_path / 'p.json').write_text(json.dumps(P), encoding='utf-8')
    salida = subprocess.run(['node', tmp_path / 'arnes.js', tmp_path / 'worker.js', tmp_path / 'geo.json', tmp_path / 'p.json'],
                            capture_output=True, text=True, check=True)
    js = json.loads(salida.stdout)

    assert js['celdas'] == int(np.count_nonzero(grilla.vertida_id))
    assert js['error'] == pytest.approx(error_py, abs=1e-9)
    assert js['pasos'] == r.n_pasos
    assert js['maximo'] == pytest.approx(r.maximo_C, abs=1e-6)
    assert js['bordes'] == pytest.approx(r.balance[-1]['energia_bordes_J'], rel=1e-9)
    for vid, media in media_py.items():
        if media is None:
            assert js['media'][str(vid)] is None
        else:
            assert js['media'][str(vid)] == pytest.approx(media, abs=1e-6)


# ---------------------------------------------------------------------------
# E10: clima del snapshot, moldaje con retiro y suelo
# ---------------------------------------------------------------------------
from datetime import datetime
from mcplan.clima import cargar_snapshot, serie_clima
from mcplan.bordes import caras_ambientales
from mcplan.suelo import SUELO_ID, agregar_suelo, modelo_suelo

ARNES_CLIMA = r"""
const fs=require('fs');const src=fs.readFileSync(process.argv[2],'utf8');
const ini=src.indexOf('function horasUTC'),fin=src.indexOf('function inicioClima');
eval(src.slice(ini,fin)+';global.serieClima=serieClima;');
const snap=JSON.parse(fs.readFileSync(process.argv[3],'utf8')),q=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
console.log(JSON.stringify(q.map(([inicio,H,m])=>{const s=serieClima(snap,inicio,H,m);return {T:s.T,V:s.V,G:s.G,fuente:s.fuente,metodo:s.metodo};})));
"""

ARNES_SUELO = r"""
const fs=require('fs');global.self={};global.postMessage=()=>{};
eval(fs.readFileSync(process.argv[2],'utf8')+';global.voxelizar=voxelizar;global.simular=simular;global.agregarSuelo=agregarSuelo;');
const G=JSON.parse(fs.readFileSync(process.argv[3],'utf8')),P=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
const v=agregarSuelo(voxelizar(G,P.d),P.suelo),r=simular(v,P,()=>{}),s=r.serie,i=s.t.length-1;
const media={};for(const id in s.pv)media[id]=s.pv[id].media[i];
console.log(JSON.stringify({forma:v.n,celdas:r.ncon,pasos:r.pasos,bordes:s.Ebor[i],gmax:s.gmax[i],media,relmax:r.relmax}));
"""


def _node(tmp_path, arnes, *archivos):
    (tmp_path / 'arnes.js').write_text(arnes, encoding='utf-8')
    salida = subprocess.run(['node', tmp_path / 'arnes.js', *archivos], capture_output=True, text=True, check=True)
    return json.loads(salida.stdout)


def test_el_visor_lleva_el_ultimo_snapshot_de_clima():
    copia = visor12.clima_del_visor(TEXTO)
    assert copia == json.loads(json.dumps(visor12.clima_desde_snapshot(RAIZ / 'data' / 'weather' / copia['nombre'])))
    assert copia['metadata']['retrieved_at'] and copia['historico']


@pytest.mark.skipif(shutil.which('node') is None, reason='Node no esta instalado')
def test_el_clima_del_visor_es_el_de_python(tmp_path):
    snap_visor = visor12.clima_del_visor(TEXTO)
    snap = cargar_snapshot(RAIZ / 'data' / 'weather' / snap_visor['nombre'])
    consultas = [('2026-11-11T08:00:00-03:00', 72, 'historico_promedio'), ('2026-11-27T08:00:00-03:00', 336, 'historico_promedio'),
                 ('2026-11-27T08:00:00-03:00', 336, 'historico_frio')]
    (tmp_path / 'main.js').write_text(re.findall(r'<script>(.*?)</script>', TEXTO, re.S)[-1], encoding='utf-8')
    (tmp_path / 'snap.json').write_text(json.dumps(snap_visor), encoding='utf-8')
    (tmp_path / 'q.json').write_text(json.dumps([[datetime.fromisoformat(i).timestamp(), h, m] for i, h, m in consultas]), encoding='utf-8')
    js = _node(tmp_path, ARNES_CLIMA, tmp_path / 'main.js', tmp_path / 'snap.json', tmp_path / 'q.json')
    for (inicio, H, metodo), j in zip(consultas, js):
        py = serie_clima(snap, inicio, H, metodo)
        assert j['metodo'] == py.metodo_cola and j['fuente'] == py.fuente.tolist()
        np.testing.assert_allclose(j['T'], py.T_aire_C, atol=1e-9)
        np.testing.assert_allclose(j['V'], py.viento_m_s, atol=1e-9)
        np.testing.assert_allclose(j['G'], py.ghi_W_m2, atol=1e-9)


@pytest.mark.skipif(shutil.which('node') is None, reason='Node no esta instalado')
def test_el_motor_del_visor_con_clima_moldaje_y_suelo_da_lo_mismo_que_python(tmp_path):
    horas, paso = 12.0, (.4, .4, .05)
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    hc = bim.cargar_yaml(RAIZ / 'config' / 'hidratacion_e07.yaml')
    h = caso['hormigon']
    snap = cargar_snapshot(RAIZ / 'data' / 'weather' / visor12.clima_del_visor(TEXTO)['nombre'])
    clima = serie_clima(snap, '2026-11-11T08:00:00-03:00', horas)
    suelo = dict(espesor=0.6, margen=0.8, k=1.5, rho=1900., cp=1000., T0=12., Tfondo=11.)
    g = agregar_suelo(voxelizar_mallas(extraer_mallas(IFC), paso, max_celdas=4_000_000),
                      espesor_m=suelo['espesor'], margen_m=suelo['margen'], T_inicial_C=suelo['T0'])
    modelos = modelos_desde_configuracion(caso, hc)
    modelos['SUELO'] = modelo_suelo(suelo['rho'], suelo['cp'])
    materiales = {'M1': Material('M1', h['conductividad_W_mK'], h['densidad_kg_m3'], h['calor_especifico_J_kgK']),
                  'SUELO': Material('SUELO', suelo['k'], suelo['rho'], suelo['cp'])}
    programa = ProgramaVertidas(g, modelos)
    hormigon = [v for v in programa.vertidas if v != SUELO_ID]
    bordes = caras_ambientales(clima, moldaje={v: (0.14, 3.0) for v in hormigon}, aislacion={v: (1.43, None) for v in hormigon},
                               suelo=dict(id=SUELO_ID, T_fondo_C=suelo['Tfondo']), absortividad=0.65, emisividad=0.9)
    r = resolver_vertidas(programa, materiales, bordes, horas * 3600, max_dt_s=600.0,
                          tiempos_guardado_s=np.arange(0, horas + 1) * 3600)
    activo = r.activa_final & (g.vertida_id != SUELO_ID)
    media_py = {f['vertida_id']: f['T_media_C'] for f in r.por_vertida if f['tiempo_h'] == horas}

    worker = re.search(r'id="worker-src">(.*?)</script>', TEXTO, re.S).group(1)
    (tmp_path / 'worker.js').write_text(worker, encoding='utf-8')
    (tmp_path / 'geo.json').write_text(json.dumps(visor12.geometria_del_visor(TEXTO)), encoding='utf-8')
    P = dict(d=list(paso), k=h['conductividad_W_mK'], rho=h['densidad_kg_m3'], cp=h['calor_especifico_J_kgK'],
             theta=45, tau=30, horizonte=horas, maxdt=600, intervalo=1, amb=dict(alfa=0.65, emis=0.9, hFijo=0),
             arriba=dict(tipo='ambiente', R=1.43, retiro=None), lados=dict(tipo='ambiente', R=0.14, retiro=3.0),
             base=dict(tipo='neumann', q=0), suelo=suelo,
             clima=dict(t=clima.t_s.tolist(), T=clima.T_aire_C.tolist(), V=clima.viento_m_s.tolist(), G=clima.ghi_W_m2.tolist()),
             vertidas=[dict(id=v, inicio=g.vertidas[v]['Inicio_h'], temp=g.vertidas[v]['TemperaturaColocacion_C'], calc=True) for v in hormigon])
    (tmp_path / 'p.json').write_text(json.dumps(P), encoding='utf-8')
    js = _node(tmp_path, ARNES_SUELO, tmp_path / 'worker.js', tmp_path / 'geo.json', tmp_path / 'p.json')

    assert js['forma'] == list(g.vertida_id.shape) and js['celdas'] == int(np.count_nonzero(g.vertida_id))
    assert js['pasos'] == r.n_pasos and js['relmax'] < 1e-8
    assert js['gmax'] == pytest.approx(float(r.T_final[activo].max()), abs=1e-6)
    assert js['bordes'] == pytest.approx(r.balance[-1]['energia_bordes_J'], rel=1e-9)
    for vid, media in media_py.items():
        if media is not None:
            assert js['media'][str(vid)] == pytest.approx(media, abs=1e-6)


# ---------------------------------------------------------------------------
# Geometría paramétrica: el navegador construye las mismas vertidas que el IFC
# ---------------------------------------------------------------------------
import copy

ARNES_GEO = r"""
const fs=require('fs');global.self={};global.postMessage=()=>{};
const main=fs.readFileSync(process.argv[2],'utf8');
const ini=main.indexOf('// --- Geometría paramétrica'),fin=main.indexOf('// --- fin geometría paramétrica');
eval(main.slice(ini,fin)+';global.generarGeometria=generarGeometria;global.validarGeometria=validarGeometria;');
eval(fs.readFileSync(process.argv[3],'utf8')+';global.voxelizar=voxelizar;');
const base=JSON.parse(fs.readFileSync(process.argv[4],'utf8')),casos=JSON.parse(fs.readFileSync(process.argv[5],'utf8'));
console.log(JSON.stringify(casos.map(g=>{const m=generarGeometria(g,base),v=voxelizar(m,[.4,.4,.05]);
  return {errores:validarGeometria(g),vol:m.map(x=>x.vol),celdas:v.tabla.map(x=>x.celdas),n:v.n};})));
"""


def _caso_y_plan(g):
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    plan = bim.cargar_yaml(RAIZ / 'config' / 'plan_v02.yaml')
    caso = copy.deepcopy(caso)
    caso['geometria'].update(lados_planta=g['lados'], distancia_caras_exteriores_m=g['caras'], espesor_borde_m=g['borde'],
                             espesor_centro_m=g['centro'], plataforma_distancia_caras_m=g['plataforma'],
                             pedestal_diametro_m=g['pedestal'], altura_pedestal_m=g['alto'],
                             altura_total_m=g['centro'] + g['alto'], lados_cilindro=g['ladoscil'])
    cotas = [(0, g['corte1']), (g['corte1'], g['borde']), (g['borde'], g['corte2']), (g['corte2'], g['centro']),
             (g['centro'], g['centro'] + g['alto'])]
    for v, (z0, z1) in zip(plan['vertidas'], cotas):
        v['z_inferior_m'], v['z_superior_m'] = z0, z1
    return caso, plan


@pytest.mark.skipif(shutil.which('node') is None, reason='Node no esta instalado')
def test_la_geometria_parametrica_del_visor_es_la_del_ifc(tmp_path):
    base = dict(lados=8, caras=18, borde=1.1, centro=3.3, plataforma=7, pedestal=7, alto=1, ladoscil=48, corte1=0.55, corte2=2.2)
    otra = dict(lados=6, caras=16, borde=1.2, centro=3.0, plataforma=6, pedestal=5.5, alto=1.4, ladoscil=32, corte1=0.6, corte2=2.0)
    casos = [base, otra]
    (tmp_path / 'main.js').write_text(re.findall(r'<script>(.*?)</script>', TEXTO, re.S)[-1], encoding='utf-8')
    (tmp_path / 'worker.js').write_text(re.search(r'id="worker-src">(.*?)</script>', TEXTO, re.S).group(1), encoding='utf-8')
    (tmp_path / 'base.json').write_text(json.dumps(visor12.geometria_del_visor(TEXTO)), encoding='utf-8')
    (tmp_path / 'casos.json').write_text(json.dumps(casos), encoding='utf-8')
    js = _node(tmp_path, ARNES_GEO, tmp_path / 'main.js', tmp_path / 'worker.js', tmp_path / 'base.json', tmp_path / 'casos.json')
    for g, j in zip(casos, js):
        caso, plan = _caso_y_plan(g)
        bim.validar_plan(caso['geometria'], plan)
        ruta = tmp_path / f"f{g['lados']}.ifc"
        bim.crear_modelo_fundacion(caso, plan).write(str(ruta))
        mallas = extraer_mallas(ruta)
        grilla = voxelizar_mallas(mallas, (.4, .4, .05), max_celdas=4_000_000)
        assert j['errores'] == []
        assert j['vol'] == pytest.approx([m.volumen_ifc_m3 for m in mallas], rel=1e-9)
        assert j['celdas'] == [f['celdas'] for f in grilla.resumen_volumenes()]
        assert j['n'] == list(grilla.vertida_id.shape)
