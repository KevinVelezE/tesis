# Reproduccion numerica del capitulo 3

Desde la raiz de `tesis`, crear el entorno e instalar dependencias:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r code/requirements.txt
```

Ejecutar la reproduccion:

```bash
python code/generate_chapter3.py
```

El script es autocontenido dentro de `tesis/code` y no importa desde `../radon` ni `../radon_intr`.
Crea figuras, tablas, metadatos y auditoria de normalizacion en `Figures/Chapter3`, `Tables/Chapter3` y `Results/Chapter3`.
