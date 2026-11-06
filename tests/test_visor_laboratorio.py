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
from mcplan.solver1d import Material, Neumann, Robin
from mcplan.solver3d import CARAS
from mcplan.solver_vertidas import resolver_vertidas
from mcplan.voxelizacion import extraer_mallas, voxelizar_mallas

RAIZ = Path(__file__).resolve().parents[1]
IFC = RAIZ / 'fundaciones' / 'fundacion.ifc'
TEXTO = (RAIZ / 'visor_laboratorio.html').read_text(encoding='utf-8')
espec = importlib.util.spec_from_file_location('visor12', RAIZ / 'scripts' / '12_actualizar_visor_laboratorio.py')
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


@pytest.mark.skipif(shutil.which('node') is None, reason='Node no esta instalado')
def test_el_motor_del_visor_da_lo_mismo_que_python(tmp_path):
    horas, paso = 12.0, (.4, .4, .05)
    caso = bim.cargar_yaml(RAIZ / 'config' / 'caso_base.yaml')
    hc = bim.cargar_yaml(RAIZ / 'config' / 'hidratacion_e07.yaml')
    h = caso['hormigon']
    grilla = voxelizar_mallas(extraer_mallas(IFC), paso, max_celdas=4_000_000)
    programa = ProgramaVertidas(grilla, modelos_desde_configuracion(caso, hc))
    materiales = {'M1': Material('M1', h['conductividad_W_mK'], h['densidad_kg_m3'], h['calor_especifico_J_kgK'])}
    bordes = {c: Robin(10.0, 15.0) for c in CARAS}
    bordes['z-'] = Neumann(0.0)
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
             theta=45, tau=30, horizonte=horas, maxdt=600, intervalo=1,
             arriba=dict(tipo='robin', h=10, R=0, T=15), lados=dict(tipo='robin', h=10, R=0, T=15), base=dict(tipo='adiabatica'),
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
