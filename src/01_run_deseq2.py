### ANÁLISIS DE RNA-SEQ CON pyDESEQ2 EN PYTHON 

# Carga de paquetes necesarios
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Para guardar figuras sin abrir ventanas
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import dotplot
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
# Para la trazabilidad
from importlib.metadata import version 
import sys 
from datetime import datetime

# Creación de las carpeta results y logs.
tarea02_python = Path.cwd() # Devuelve el directorio de trabajo actual (Current Working Directory)
(tarea02_python/"results").mkdir(parents=True, exist_ok=True)
(tarea02_python/"logs").mkdir(parents=True, exist_ok=True)
# parents = True para que se creen las carpeta que contiene results y logs si no está creada.
# exist_ok = True para que no devuelva error si la carpeta ya está creada.

# Definición de variables 
counts_path = Path("data/raw/counts.csv.gz")
metadata_path = Path("data/raw/metadata.csv.gz")

condition_column = "site"  
plate_column = "plate"  
condition_test = "Primary_Solid_Tumor"
condition_ref = "Solid_Tissue_Normal"

alpha = 0.05
lfc_threshold = 1

# Análisis de expresión diferencial
print("[1/8] Leyendo datos de entrada (TCGA-BRCA)...")
counts_df = pd.read_csv(counts_path, index_col =0) 
# filas = genes, columnas = muestras
metadata = pd.read_csv(metadata_path, index_col=0)
# filas = muestras, columnas= "site" y "plate"

print("[2/8] Alineando muestras entre conteos y metadata...")
metadata = metadata.loc[list(counts_df.columns)]
metadata[plate_column] = metadata[plate_column].astype("category")

# Establecimiento del orden de las categorías en metadata para especificar la condición de referencia.
metadata[condition_column] = pd.Categorical(
    metadata[condition_column],
    categories=[condition_ref, condition_test]
)

print("[3/8] Creando objeto DESeq2 y filtrando genes de baja cuenta...")
# Filtrado de genes
min_count = 20
min_samples = round(len(counts_df.columns)* 0.10) # 10% del total de muestras.

genes_to_keep = ((counts_df >= min_count).sum(axis=1)) >= min_samples
# counts_df >= min_count devuelve una matriz booleana que indica en qué pares (muestra,gen),
# las lecturas del gen son mayor que el número mínimo de lecturas. 
# (counts_df >= min_count).sum(axis=1) suma los "True" de la matriz booleana anterior, que
# será el número de muestras para cada gen que supera el número mínimo de lecturas.
# Con >= min_samples restringimos que el número mínimo de muestras anterior.
counts_df = counts_df.loc[genes_to_keep]

# pydeseq2 espera muestras en filas y genes en columnas, luego usamos la transposición.
dds = DeseqDataSet(
    counts = counts_df.T, 
    metadata = metadata,
    design = f"~ {condition_column} + {plate_column}",
)

print("[4/8] Ejecutando DESeq...")
dds.deseq2()

normalized_counts = pd.DataFrame(dds.layers["normed_counts"], 
                                 index= dds.obs_names, 
                                 columns = dds.var_names)
# filas = muestras, columnas = genes

normalized_counts.to_csv("results/normalized_counts.tsv", 
                         sep="\t")

dds.plot_dispersions(save_path=str("results/plot_dispersion.png"))
plt.close("all")

print("[5/8] Extrayendo resultados... ")
res = DeseqStats(dds, 
                 contrast = [condition_column, condition_test, condition_ref],
                 alpha = alpha)

res.summary()
res_df = res.results_df.copy()
# Al filtrar un DataFrame, a veces pandas devuelve una vista (una referencia al DataFrame original) en lugar de un nuevo objeto independiente.
# Por tanto, puede no guardar los cambios correctamente en el objeto filtrado. Al usar .copy(), se le fuerza a crear una copia independiente 
# en memoria para ese nuevo DataFrame.

res_df["gene_id"] = res_df.index
res_df = res_df.sort_values("padj")

has_padj = res_df["padj"].notna()
sig_df = res_df[has_padj & (res_df["padj"] < alpha)].copy()
de_df = res_df[has_padj & (res_df["padj"] < alpha) & 
    (res_df["log2FoldChange"].abs() > lfc_threshold)].copy()
up_df = de_df[de_df["log2FoldChange"] > 0]
down_df = de_df[de_df["log2FoldChange"] < 0]

top10_up = up_df.head(10)
top10_down = down_df.head(10)

de_df["top10_label"] = "other"
de_df.loc[de_df["gene_id"].isin(top10_up["gene_id"]),"top10_label"] = "top10_up"
# Acceso por [fila,columna]
de_df.loc[de_df["gene_id"].isin(top10_down["gene_id"]),"top10_label"] = "top10_down"

# Guardado de las tablas en ficheros TSV.
res_df.to_csv("results/deseq2_all_results.tsv", sep="\t")
sig_df.to_csv("results/deseq2_sig_padj_0.05.tsv", sep="\t")
de_df.to_csv("results/deseq2_de_padj_0.05_lfc_1.tsv", sep="\t")
top10_up.to_csv("results/top10_up_genes.tsv", sep="\t")
top10_down.to_csv("results/top10_down_genes.tsv", sep="\t")

print("[6/8] Generando gráficos... ")

# TOP 10
# Se ordena top10_up y top10_down por log2FoldChange para representarlos.
top10_up = top10_up.sort_values("log2FoldChange", ascending=False)
# ascending = False -> De mayor a menor
top10_down = top10_down.sort_values("log2FoldChange", ascending=True)

fig, axs = plt.subplots(1, 2, figsize=(10, 4)) 
axs[0].barh(top10_up["gene_id"],top10_up["log2FoldChange"], color= "darkred") 
axs[0].set_title("Up (Tumor > Normal)")
axs[0].set_xlabel("log2FoldChange")
axs[0].invert_yaxis()  # Para que aparezca el gen con mayor Log2FC arriba

axs[1].barh(top10_down["gene_id"],top10_down["log2FoldChange"], color= "darkblue") 
axs[1].set_title("Down (Tumor < Normal)")
axs[1].set_xlabel("log2FoldChange")
axs[1].invert_yaxis()  

fig.tight_layout()
fig.savefig("results/plot_top10_genes.png", dpi=150)
plt.close(fig)

# VOLCANO PLOT
volcano_df = res_df[has_padj].copy()

volcano_df["status"] = "No significativo" 
volcano_df.loc[volcano_df["gene_id"].isin(up_df["gene_id"]), "status"] = "Up"
volcano_df.loc[volcano_df["gene_id"].isin(down_df["gene_id"]), "status"] = "Down"

volcano_df["-log10padj"] = -np.log10(volcano_df["padj"].clip(lower=np.finfo(float).tiny))
# np.finfo(float).tiny devuelve el menor número positivo normal representable con el tipo float de numpy.
# .clip(lowe=...) establece un límite para que cualquier valor de padj que sea menor que el número 
# definido anteriormente, se sustituya por ese mismo. 
# Se utiliza porque si algún padj es cero, daría infinito y no podría representarse correctamente.

colors = {"Up": "darkred", "Down": "darkblue", "No significativo": "gray"}

fig, ax = plt.subplots(figsize=(7,6))

for status, color in colors.items():
    subset = volcano_df[volcano_df["status"] == status]
    ax.scatter(
        subset["log2FoldChange"],
        subset["-log10padj"],
        c = color,
        label = status,
        alpha = 0.7, # Define la opacidad.
        s = 4
    )

ax.axvline(-lfc_threshold, linestyle="--", color="black", linewidth=1)
ax.axvline(lfc_threshold, linestyle="--", color="black", linewidth=1)
ax.axhline(-np.log10(alpha), linestyle="--", color="black", linewidth=1)

labelled = volcano_df[volcano_df["gene_id"].isin(pd.concat([top10_up, top10_down])["gene_id"])]
# Se concatenan los DataFrames top10_up y top10_down, para definir que los gene_id en ambos 
# Data Frames serán los que se etiquetarán.
for _, row in labelled.iterrows():
    # iterrows() permite recorrer las filas del DataFrame una por una.
    ax.text(row["log2FoldChange"],  # Posición eje X
            row["-log10padj"],      # Posición eje Y
            row["gene_id"],         # Texto 
            fontsize=6)

ax.set_xlabel("log2FoldChange")
ax.set_ylabel("-log10(padj)")
ax.legend(title="status", markerscale=3)
# markerscale=3 aumenta tres veces el tamaño de los símbolos de la leyenda respecto a los puntos originales.

fig.tight_layout()
fig.savefig("results/plot_volcano.png", dpi=150)
plt.close(fig)

# MA PLOT
ma_df = res_df[res_df["baseMean"] > 0]
# baseMean contiene la expresión media normalizada, y seleccionamos aquellos genes en los 
# que esta sea mayor que cero.
is_sig = ma_df["padj"].notna() & (ma_df["padj"] < alpha)  # Genes significativos
no_sig = ma_df["padj"].notna() & (ma_df["padj"] >= alpha) # Genes no significativos

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(ma_df.loc[no_sig, "baseMean"],        # x = expresión media
           ma_df.loc[no_sig, "log2FoldChange"],  # y = log2FoldChange
           s=4, c="grey") # s=4 establece un tamaño pequeño para el punto

ax.scatter(ma_df.loc[is_sig, "baseMean"],
           ma_df.loc[is_sig, "log2FoldChange"],
           s=4, c="blue")

ax.axhline(0, color="black", linewidth=1)
ax.set_xscale("log")

ax.set_xlabel("mean of normalized counts")
ax.set_ylabel("log fold change")
ax.set_title("pyDESeq2 MA plot")

fig.tight_layout()
fig.savefig("results/plot_ma.png", dpi=150)
plt.close(fig)

print("[7/8] Contextualizacion biologica: realizando enriquecimiento GO y Reactome...")

# Gen destacado
gen_destacado = de_df.iloc[[0]][["gene_id", "log2FoldChange", "padj"]]
# Se usan dobles corchetes para seleccionar el subconjunto de columnas que queremos,
# manteniendo el formato de DataFrame, sin que se convierta en una Series de Pandas.
gen_destacado.to_csv("results/top_gene_context.txt", sep="\t")

# Listas de genes para el enriquecimiento.
genes_de = list(de_df["gene_id"].unique())  # Universo
genes_up = list(up_df["gene_id"].unique())
genes_down = list(down_df["gene_id"].unique())

# Enriquecimiento en términos GO - Biological Process
# gp.get_library_name("human") -> Para consultar las posibles librerías que podemos poner en gene_sets.

def enrichment(genes, library, output_name):
    """
    Función para realizar el enriquecimiento funcional de genes.
    Recibe una lista de genes, el nombre de una librería concreta 
    de Enrichr y el nombre de salida (ambos entre comillas); para 
    generar los archivos: <output_name>.tsv y plot_<output_name>.png.
    """
    enr = gp.enrichr(
        gene_list = genes,
        gene_sets = library,
        background = genes_de, # El universo queda fijo => Genes diferencialmente expresados.
        outdir = None,  # Guarda los resultados solo en la memoria.
        cutoff = alpha
    )

    table = enr.res2d[enr.res2d["Adjusted P-value"] < alpha] # Términos significativos
    table.to_csv(f"results/{output_name}.tsv", sep= "\t", index = False)
    # .res2d = Atributo del objeto de enriquecimiento de gseapy para obtener la tabla de resultados.
    
    dotplot(table,   
            column = "Adjusted P-value", 
            title = output_name, 
            size=5,
            top_term = 10,  # showCategory en R
            figsize=(3, 5),
            ofname= f"results/plot_{output_name}.png"
            )
    plt.close("all")
    
enrichment(genes_up, "GO_Biological_Process_2026", "go_bp_up")
enrichment(genes_down, "GO_Biological_Process_2026", "go_bp_down")
enrichment(genes_up, "Reactome_2022", "reactome_up")
enrichment(genes_down, "Reactome_2022", "reactome_down")

print("[8/8] Generando informes de trazabilidad...")

packages = ["pydeseq2", "pandas", "numpy", "matplotlib", "gseapy"]
versions = [f"{package}:{version(package)}" for package in packages] # Lista de f-strings

session_info = [f"Version de Python: {sys.version}"]
session_info.extend(versions) # Añade cada uno de los elementos de versions a la lista session_info.

Path("logs/sessionInfo.txt").write_text("\n".join(session_info)); 
# "\n".join(run_summary) para poner cada elemento de la lista run_summary
# en una línea, porque \n representa el salto de linea.
# ; para que no devuelva el número de caracteres escritos.

run_summary = [
    "RESUMEN DEL ANALISIS DE RNA-SEQ",
    "Proyecto: TCGA-BRCA",
    f"Fecha de ejecucion: {datetime.now()}",
    f"Version de Python: {sys.version}",
    f"Entrada total de genes: {len(counts_df)}",
    f"Genes analizados: {dds.n_vars}",
    f"Muestras analizadas: {dds.n_obs}",
    f"Genes significativos (padj < {alpha}): {len(sig_df)}",
    f"Genes DE (padj < {alpha} y |log2FC| > {lfc_threshold}): {len(de_df)}",
    f"Genes sobreexpresados: {len(up_df)}",
    f"Genes subexpresados: {len(down_df)}",
    f"Condicion de interes: {condition_test}",
    f"Condicion de referencia: {condition_ref}",
    f"Columna de condicion: {condition_column}",
    f"Columna de lote: {plate_column}",
    f"Umbral de significacion: {alpha}",
    f"Umbral absoluto de log2FC: {lfc_threshold}",
    "Metodo: PyDESeq2",
]

Path("logs/run_summary.txt").write_text("\n".join(run_summary));
