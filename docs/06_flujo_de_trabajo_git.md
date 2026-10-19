# 06 — Flujo de trabajo con Git y GitHub

> Objetivo: que el repositorio esté **siempre ordenado y funcionando**, y que cada entrega
> apunte a un commit exacto (el curso evalúa sobre un SHA congelado).

## 1. Reglas del grupo

1. **`main` siempre funciona**: los tests pasan y los scripts corren.
2. Nadie trabaja directo en `main`. Cada tarea va en una **rama** y entra por **Pull Request (PR)**.
3. **Otro integrante revisa** cada PR antes de unirlo. Así todos conocen todo el código
   (en el oral se le puede preguntar a cualquiera por cualquier parte).
4. Antes de cada commit: `pytest`.
5. Los archivos generados (IFC v01, figuras) **no se editan a mano**: se regeneran con los scripts.

## 2. Nombres

- **Ramas**: `eXX-descripcion-corta` (p. ej. `e04-solver-2d`, `e08-voxelizacion`).
- **Commits**: `tipo: descripción en presente`, con tipos

| tipo | para |
|---|---|
| `feat` | funcionalidad nueva (`feat: agregar condición Robin al solver 1D`) |
| `fix` | corrección de un error |
| `test` | tests nuevos o modificados |
| `docs` | documentación |
| `refactor` | reorganizar código sin cambiar resultados |
| `chore` | configuración, dependencias, archivos auxiliares |

## 3. Ciclo de una tarea

```bash
git switch main
git pull                                  # traer lo último del grupo
git switch -c e04-solver-2d               # rama nueva
# ... trabajar ...
pytest                                    # todo debe pasar
git add <archivos>
git commit -m "feat: solver 2D por balance de celdas"
git push -u origin e04-solver-2d          # subir la rama
# en GitHub: abrir Pull Request -> revisión de otro integrante -> Merge
```

## 4. Entregas: commit evaluado

Cada actividad declara:

```
Repositorio:     https://github.com/aespinoza6-svg/P2-Grupo7
Commit evaluado: https://github.com/aespinoza6-svg/P2-Grupo7/commit/<SHA>
```

1. Unir a `main` todo lo de la entrega.
2. Obtener el SHA: `git switch main && git pull && git rev-parse HEAD`.
3. Marcarlo con una **etiqueta** (tag) para encontrarlo después:

```bash
git tag -a E01 -m "Entrega E01: BIM + solver 1D"
git push origin E01
```

4. Copiar el texto de `entregas/EXX.md` en Canvas con el SHA.

Una vez iniciada la revisión **ese SHA no cambia**: se puede seguir trabajando en `main`,
pero la evaluación es sobre el commit etiquetado.

## 5. Qué se versiona y qué no

| Sí | No (ver `.gitignore`) |
|---|---|
| código (`src/`, `scripts/`, `tests/`) | entornos virtuales (`.venv/`) |
| configuración (`config/`) | `__pycache__/`, `.pytest_cache/` |
| modelo IFC generado (`modelos/`) | respaldos de Blender (`*.blend1`) |
| figuras y tablas livianas (`resultados/`) | resultados pesados (Monte Carlo, campos 3D) — se regeneran |
| documentación (`docs/`, `entregas/`) | |
