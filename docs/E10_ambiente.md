# E10: el ambiente real en la física

Objetivo de la actividad: mostrar cómo el clima entra a la ecuación del calor y cómo se obtienen pronósticos trazables.

## Qué tiene el repo

| Lo que pide el curso | Dónde está | Qué hace |
|---|---|---|
| Módulo de condiciones de borde | `mcplan/bordes.py` | Arma el borde de cada cara expuesta con aire, viento, sol y moldaje o aislación |
| Tratamiento de suelo | `mcplan/suelo.py` | Agrega un bloque de suelo con capacidad térmica bajo la fundación |
| Moldaje y aislación | `mcplan/bordes.py` y `config/ambiente_e10.yaml` | Resistencia R en serie, que se retira en la hora `RetiroMoldaje_h` de cada vertida |
| Cliente o loader meteorológico | `mcplan/clima.py` y `scripts/13_descargar_pronostico.py` | Descarga Open-Meteo y lee lo guardado |
| Esquema de snapshots | `data/weather/<fecha de descarga>/` | `metadata.json`, `forecast.json` e `historico.json`, sin sobrescribir |
| Simulación con clima real | `scripts/14_simular_ambiente.py` | Corre el caso con todo lo anterior y deja figuras y `resumen.json` trazable |
| Tests | `tests/test_ambiente.py` | Corren sin internet, con un snapshot de ejemplo |

## Comandos

Desde `P2-Grupo7`, con el entorno activado:

```powershell
python -X utf8 scripts/13_descargar_pronostico.py
python -X utf8 scripts/14_simular_ambiente.py --constante --salida resultados/e10/constante
python -X utf8 scripts/14_simular_ambiente.py --comparar-con resultados/e10/constante
python -X utf8 scripts/14_simular_ambiente.py --inicio 2026-11-11T08:00:00-03:00 --salida resultados/e10/pronostico_11nov
python -X utf8 -m pytest tests/test_ambiente.py -v
```

Para cambiar el moldaje, la aislación, el suelo, la fecha o la forma de completar el clima, se edita `config/ambiente_e10.yaml`. No se toca el código.

## Cómo entra todo en una sola fórmula

Cada cara de celda que toca el aire recibe:

$$q = \frac{T_{ref}(t) - T_{celda}}{\dfrac{\Delta x}{2k} + R(t) + \dfrac{1}{h(t)}}\,A$$

| Término | Qué es | De dónde sale |
|---|---|---|
| $\Delta x/2k$ | Conducción desde el centro de la celda hasta su cara | Malla y hormigón |
| $R(t)$ | Moldaje o aislación; vale 0 después del retiro | Configuración e IFC |
| $1/h(t)$ | Película de aire | Viento y radiación |
| $T_{ref}(t)$ | Aire, o sol-aire en las caras que miran al cielo | Snapshot meteorológico |

## Qué cambia con el clima real

Son 14 días, con malla de 20 × 20 × 5 cm. El diferencial es la aproximación núcleo − superficie dentro de cada vertida.

| Caso | T máx | Cuándo | Diferencial máx. | Cuándo |
|---|---|---|---|---|
| Ambiente constante de E07: aire a 15 °C, h = 10, sin sol, sin suelo ni moldaje | 62,16 °C | 177 h | 46,1 °C | 249 h |
| Inicio el 27-nov, con relleno histórico, moldaje de madera retirado a las 72 h y suelo | 61,94 °C | 149 h | 49,3 °C | 143 h |
| Inicio el 11-nov, todo con pronóstico, con moldaje y suelo | 62,00 °C | 149 h | 47,7 °C | 189 h |

Figuras en `resultados/e10/<caso>/`: `clima.png`, `temperaturas_y_diferencial.png` y `cortes_con_suelo.png`.

- **El diferencial oscila cada día:** de noche el aire está a unos 11 °C y la superficie se enfría; de día el sol la calienta y el diferencial baja.
- **Con clima real, el diferencial máximo sube unos 3 °C:** las noches son más frías que los 15 °C constantes y hay viento.
- **El suelo saca calor por la base:** después del peak, la T máxima baja más rápido que con base adiabática.
- **El sol calienta la cara superior de V01 en las primeras horas:** la T máxima tiene un salto a las 6 h, que es la cara al sol, no el núcleo.
- **Ningún caso cumple 20 °C sin protección:** hay que aislar, como mostró la comparación de estrategias del visor.

## Las preguntas de la actividad

### ¿Qué término usa la temperatura ambiente?

El término **Robin**, o de convección: $q = h\,(T_{aire} - T_{superficie})$. La temperatura ambiente es $T_{ref}$ de la fórmula de arriba. Ahora cambia hora a hora según el snapshot, y antes era 15 °C fijo.

### ¿Cómo entra el viento?

A través de $h$. Más viento renueva más rápido el aire junto a la superficie, y el hormigón pierde calor más rápido. Usamos la correlación clásica $h_{conv} = 5{,}7 + 3{,}8\,V$ [W/m²K], con $V$ a 10 m en m/s. Es aproximada: la medición a 10 m no es el viento junto a la cara.

### ¿Cómo modelar la resistencia térmica de un moldaje?

Como una resistencia en serie: $R = e/k_{material}$. Por ejemplo, madera contrachapada de 18 mm: $R = 0{,}018/0{,}13 \approx 0{,}14$ m²K/W. Poliestireno de 5 cm: $0{,}05/0{,}035 \approx 1{,}43$ m²K/W. Se suma a $\Delta x/2k$ y a $1/h$, y por eso reduce el flujo. No tiene capacidad térmica: es una simplificación, porque la masa del moldaje es pequeña frente al hormigón.

### ¿Qué cambia al desmoldar?

$R$ pasa a 0 en la hora de retiro de **cada** vertida, que viene de `RetiroMoldaje_h` en el IFC. En ese instante la superficie queda expuesta al aire y se enfría más rápido. El diferencial núcleo–superficie puede **subir de golpe**: es el choque térmico. Por eso el retiro es una decisión del plan y no un detalle.

### ¿Qué es un flujo solar?

Es la energía del sol que llega a una superficie, en W/m². Open-Meteo la entrega como `shortwave_radiation` (GHI, sobre un plano horizontal). La superficie absorbe la fracción $\alpha$, que usamos como 0,65 para hormigón gris. La aplicamos en las caras que miran al cielo; las caras verticales no reciben sol en este modelo.

### ¿Cómo evitar doble-contar radiación y convección?

Hay que separar qué hace cada término:

- $h = h_{conv} + h_{rad}$ multiplica **diferencias de temperatura**: es la convección más la onda larga linealizada, $h_{rad} = 4\varepsilon\sigma T^3 \approx 4{,}9$ W/m²K.
- El sol es un **aporte aparte**, $\alpha G$, que no depende de la temperatura de la superficie. Lo metemos como **temperatura sol-aire**: $T_{sa} = T_{aire} + \alpha G/h$.

Así, $h(T_{sa} - T) = h(T_{aire} - T) + \alpha G$: cada término aparece una sola vez.

El error típico es usar un $h$ «combinado» que ya incluye la radiación y además sumarle $h_{rad}$, o meter el sol como Neumann y además como sol-aire. Con viento nulo, nuestro $h = 5{,}7 + 4{,}9 \approx 10{,}6$ W/m²K, coherente con el $h = 10$ constante que usábamos antes.

### ¿Qué diferencia hay entre forecast, reanalysis y observación?

| Tipo | Qué es | En el repo |
|---|---|---|
| **Observación** | Medición de una estación en el sitio | No la usamos todavía |
| **Reanalysis** | Reconstrucción del pasado que combina un modelo con observaciones (ERA5) | `historico.json`, de `archive-api` |
| **Forecast** | Predicción del futuro hecha por un modelo, que depende del día en que se descargó | `forecast.json` |

### ¿Qué ocurre si el forecast no cubre los 14 días después del 27?

Pasa de verdad. El snapshot del 10-nov llega solo hasta el 25-nov, así que **no cubre ninguna hora** desde el 27-nov. No es un error: la guía G07 dice que es parte del problema.

Lo resolvemos completando la cola con el histórico ERA5 de los **mismos días y horas** de 2021 a 2025. Cada hora queda marcada como `pronostico` o `historico`, y las figuras sombrean el relleno. Hay dos métodos:

- **`historico_promedio`:** el promedio de los cinco años.
- **`historico_frio`:** el año más frío, como escenario conservador para el diferencial.

A medida que se acerque la fecha, cada snapshot nuevo cubrirá más horas con pronóstico.

### ¿Por qué se deben conservar los forecasts antiguos?

1. **Trazabilidad:** cada resultado dice con qué snapshot se calculó (`resumen.json` guarda su nombre, su fecha de descarga y su huella SHA-256). Si se sobrescribe, ya no se puede repetir la figura.
2. **Defensa de una decisión:** un plan se decide con lo que se sabía ese día. Sirve para mostrar por qué se decidió así, aunque el clima después cambie.
3. **Medir el error del pronóstico:** comparar lo que se pronosticó con lo que pasó, según la anticipación, permite estimar la incertidumbre meteorológica para Monte Carlo.

Por eso `descargar_snapshot` falla si la carpeta ya existe: un snapshot no se sobrescribe.

## Límites declarados

- El cielo se supone a la temperatura del aire. De noche el cielo despejado está más frío, y el enfriamiento real sería mayor.
- Las caras verticales no reciben sol.
- $h(V)$ es una correlación aproximada.
- El suelo tiene valores típicos, sin ensayos del sitio. Su fondo está a 12 °C fijo a 2 m, y sus costados son adiabáticos.
- El diferencial núcleo–superficie se aproxima por celdas dentro de cada vertida, no a 25 mm.
- La hidratación sigue sin calibrar y no depende de la temperatura.
