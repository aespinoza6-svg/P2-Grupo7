"""Comparar dos corridas con mismo dominio/fisica/horizonte y distinto dt.
No ejecuta nuevas simulaciones. Requiere salidas de script 09.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RAIZ=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,default=RAIZ/'resultados'/'e07_simulacion_acoplada')
    p.add_argument('--refinada',type=Path,default=RAIZ/'resultados'/'e07_sensibilidad_dt300')
    p.add_argument('--salida',type=Path,default=RAIZ/'resultados'/'e07_precision_temporal')
    a=p.parse_args()
    if not (a.base/'resumen.json').is_file() or not (a.refinada/'resumen.json').is_file():
        p.error('Faltan las corridas base y refinada. Ejecutelas primero con el script 09.')
    jb=json.loads((a.base/'resumen.json').read_text(encoding='utf-8'))
    jf=json.loads((a.refinada/'resumen.json').read_text(encoding='utf-8'))
    for campo in ['grilla_sha256','hidratacion','horizonte_h']:
        if jb[campo]!=jf[campo]:p.error(f'No se puede atribuir la diferencia solo al paso: cambio {campo}.')
    tb={k:v for k,v in jb['termica'].items() if k!='max_dt_s'}
    tf={k:v for k,v in jf['termica'].items() if k!='max_dt_s'}
    if tb!=tf:p.error('Las condiciones termicas de ambas corridas no coinciden.')
    with np.load(a.base/'resultado_termico.npz',allow_pickle=False) as b,np.load(a.refinada/'resultado_termico.npz',allow_pickle=False) as f:
        if not np.array_equal(b['tiempo_s'],f['tiempo_s']):p.error('No coinciden los tiempos de salida.')
        if not np.array_equal(b['activa_final'],f['activa_final']):p.error('No coinciden las celdas activas finales.')
        mask=b['activa_final'];error_final=np.abs(b['T_final_C'][mask]-f['T_final_C'][mask])
        tiempos=b['tiempo_s']/3600
        cortes_b=b['cortes_xz_C'].copy();cortes_f=f['cortes_xz_C'].copy()
    with (a.base/'balance_energia.csv').open(encoding='utf-8-sig') as f:hb=list(csv.DictReader(f))
    with (a.refinada/'balance_energia.csv').open(encoding='utf-8-sig') as f:hf=list(csv.DictReader(f))
    if len(hb)!=len(tiempos) or len(hf)!=len(tiempos):p.error('Historias incompletas.')
    max_b=np.array([float(x['T_max_C']) for x in hb]);max_f=np.array([float(x['T_max_C']) for x in hf])
    info=dict(alcance='Sensibilidad temporal a grilla fija; no prueba convergencia espacial ni validacion experimental.',
              horizonte_h=jb['horizonte_h'],max_dt_base_s=jb['termica']['max_dt_s'],
              max_dt_refinado_s=jf['termica']['max_dt_s'],
              Tmax_base_C=jb['maximo_global_C'],Tmax_refinada_C=jf['maximo_global_C'],
              diferencia_Tmax_C=abs(jb['maximo_global_C']-jf['maximo_global_C']),
              diferencia_maxima_campo_final_C=float(error_final.max()),
              diferencia_media_absoluta_campo_final_C=float(error_final.mean()),
              diferencia_maxima_historia_Tmax_C=float(np.max(abs(max_b-max_f))),
              residuo_balance_relativo_base=jb['residuo_balance_relativo_max'],
              residuo_balance_relativo_refinado=jf['residuo_balance_relativo_max'])
    a.salida.mkdir(parents=True,exist_ok=True)
    (a.salida/'comparacion.json').write_text(json.dumps(info,indent=2,ensure_ascii=False),encoding='utf-8')
    with (a.salida/'comparacion_por_hora.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['tiempo_h','Tmax_base_C','Tmax_refinada_C','diferencia_C'])
        w.writerows(zip(tiempos,max_b,max_f,max_b-max_f))
    fig,axs=plt.subplots(2,1,figsize=(10,7),layout='constrained')
    axs[0].plot(tiempos,max_b,label=f'Paso máximo {info["max_dt_base_s"]:g} s',color='#287bb5')
    axs[0].plot(tiempos,max_f,label=f'Paso máximo {info["max_dt_refinado_s"]:g} s',color='#d97922',ls='--')
    axs[0].set(ylabel='Máxima del dominio [°C]',title='Sensibilidad temporal del solver acoplado');axs[0].legend()
    axs[1].plot(tiempos,max_b-max_f,color='#985bad');axs[1].set(ylabel='Diferencia entre máximas [°C]',title=f'Diferencia máxima en el campo final: {error_final.max():.4f} °C')
    for ax in axs:ax.set_xlabel('Tiempo [h]');ax.grid(alpha=.2)
    fig.savefig(a.salida/'comparacion_temporal.png',dpi=150);plt.close(fig)
    print(json.dumps(info,indent=2,ensure_ascii=False),flush=True)
    print(f'Resultados: {a.salida}',flush=True)


if __name__=='__main__':main()
