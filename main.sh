# Situarse en la carpeta del proyecto:
cd "$(dirname "$0")" || exit 1
# $0 obtiene la ruta con la que se ejecuta el main.sh, y dirname extrae
# la ruta hacia la carpeta que contiene main.sh para fijar el directorio
# de trabajo.
# || exit 1 significa que, si no se puede cambiar de carpeta, main.sh 
# termina indicando un error.

# Inicializar Conda para poder usarlo en este fichero:
eval "$(conda shell.bash hook)"

# Activar el entorno virtual de Conda del proyecto:
conda activate ./.venv

# Ejecutar el script de Python que realiza el análisis RNA-Seq:
python src/01_run_deseq2.py

