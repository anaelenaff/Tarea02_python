# Análisis de RNA-seq con PyDESeq2 (Python) 
**Autora:** Ana Elena Fernández Fernández, Cuarto curso del Grado en Bioquímica, Asignatura de Bioinformática y Análisis Genómico  
**Lenguaje elegido:** Python  
**Enlace al repositorio de GitHub del proyecto:** https://github.com/anaelenaff/Tarea02_python.git

## Estructura del proyecto

```
tarea02_python/
|-- README.md
|-- environment.yml  
|-- main.sh
|-- informe.pdf
|-- data/raw/           # counts.csv.gz, metadata.csv.gz (del aula virtual)
|-- src/01_run_deseq2.py 
|-- results/            # tablas .tsv  y figuras .png               
|-- logs/               # run_summary.txt y sessionInfo.txt
|-- notebook/pruebas_01_run_deseq2.ipynb
```

## Instrucciones
### 1. Creación del entorno de Conda
Desde la carpeta del proyecto (tarea02_python):
```bash
conda env create -p ./.venv -f environment.yml
```
Se creará así un entorno de conda que contendrá los paquetes necesarios para el análisis y sus dependencias. Fundamentalmente:
- python=3.12 
- pandas 
- numpy 
- matplotlib 
- pydeseq2 
- gseapy  

Además, el entorno incluye el paquete jupyterlab, puesto que el proyecto contiene un notebook de JupyterLab con el que se probaron los bloques de código por separado antes de lanzar el script.

Para activar el entorno:
```bash
conda activate ./.venv
```
**NOTA:** A la hora de ejecutar el análisis, el archivo main.sh ya contiene la instrucción para activar el entorno, pero si se desea consultar el notebook, el entorno virtual debe activarse previamente con el comando anterior.

### 2. Ejecución del análisis 
```bash
bash main.sh
```
Requiere conexión a internet, pues el enriquecimiento con gseapy consulta Enrichr.

### 3. Resultados
- `results/`: tablas de DESeq2, conteos normalizados, enriquecimiento GO y Reactome (.tsv) y figuras (.png).
- `logs/`: versiones de las librerías (sessionInfo.txt) y `run_summary.txt` con las cifras clave.

## Comentarios
- El script se desarrolló en JupyterLab (`notebook/`) para probar cada bloque.
- Git se inicializó después de terminar el script, por lo que su elaboración no queda versionada.
- Las cifras obtenidas difieren ligeramente de las de control de la guía. Tras el filtro se obtienen 17615 genes (guía: 17606), y 4133 genes DE (guía: 4054), de los cuales 2166 están sobreexpresados y 1967 subexpresados (guía: 2115 y 1939). Se probaron valores diferentes de `min_samples` pero ninguno reprodujo las cifras de la guía, por lo que el análisis se realizó con `min_samples = 8` . No obstante, las cifras del run_summary son prácticamente idénticas en el análisis con R y con Python, por lo que la diferencia no se debe a un fallo de traducción de R a Python.



