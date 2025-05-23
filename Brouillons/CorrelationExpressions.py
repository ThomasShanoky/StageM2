import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt



expression_kmers_file = "ScriptsPrincipaux/ScriptFinal/ExpressionsWithKmers.csv"
expression_beataml_file = "ScriptsPrincipaux/ScriptFinal/BEATAML_Expressions.csv"

expression_kmers = pd.read_csv(expression_kmers_file, sep=",")
expression_kmers.set_index("Gene", inplace=True)
expression_beataml = pd.read_csv(expression_beataml_file, sep=",")
expression_beataml.set_index("Gene", inplace=True)


##### Intersection des deux listes d'échantillons #####

for sample, col in expression_kmers.items():
    if sample not in expression_beataml.columns:
        expression_kmers.drop(sample, axis=1, inplace=True)

for sample, col in expression_beataml.items():
    if sample not in expression_kmers.columns:
        expression_beataml.drop(sample, axis=1, inplace=True)

# print(expression_kmers.head())
# print(expression_beataml.head())


##### Choix d'un gène et plot des deux expressions #####

def plot_expressions(gene, expr_ba, expr_km, plot=True):
    gene_expr_ba = expr_ba.loc[gene]
    gene_expr_km = expr_km.loc[gene]

    gene_expr_ba = gene_expr_ba.sort_values(ascending=True)
    gene_expr_km = gene_expr_km.reindex(gene_expr_ba.index)

    # gene_expr_ba = np.exp(gene_expr_ba)

    pente, origine, r, pval, std = linregress(gene_expr_ba, gene_expr_km)
    regre_line = pente*gene_expr_ba + origine

    if plot:
        plt.scatter(gene_expr_ba, gene_expr_km, alpha=0.5)
        plt.plot(gene_expr_ba, regre_line, color="red", 
                 label=f"y = {pente:.2f}x + {origine:.2f}\nr = {r:.2f}, p = {pval:.2e}")
        plt.xlabel("Expression BEATAML")
        plt.ylabel("Expression kmers")
        plt.title(f"Expressions de {gene}")
        plt.grid("--")
        plt.legend()
        plt.show()
    
    return r, pval, std

Genes = list(expression_beataml.index)
Genes.sort()


R = []
p_values = []
stds = []

for gene in Genes:
    r, p, std = plot_expressions(gene, expression_beataml, expression_kmers, True)
    R.append(r)
    p_values.append(p)
    stds.append(std)

# print(Genes)
# print("r")
# print([f"{r:.5f}" for r in R_carres])
# print("p")
# print([f"{p:.5e}" for p in p_values])

print(f"Moyenne des coefficients de corrélation : {np.mean(R):.5f}")