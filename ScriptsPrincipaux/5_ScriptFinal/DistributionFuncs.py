import os
import numpy as np
import pandas as pd
import scipy.stats as stats



def get_number_for_bar(self, cat_name, ind_beataml, ind_list, feature):
    number_for_bar = [0 for _ in range(len(cat_name))]
    
    for ind in ind_list:
        index_ind = ind_beataml.index(ind)
        OneCat = self.data_beat_aml.iloc[index_ind][feature]

        for i, cat in enumerate(cat_name):
            if OneCat == cat:
                number_for_bar[i] += 1
    return number_for_bar


def get_number_for_bar_spe(self, cat_name:list[str], ind_beataml, ind_list:list[str], feature:str)->list[int]:
    number_for_bar = [0, 0]
    
    for ind in ind_list:
        index_ind = ind_beataml.index(ind)
        OneCat = self.data_beat_aml.iloc[index_ind][feature]

        if pd.isna(OneCat) or OneCat == 'nan':
            number_for_bar[1] += 1
        else:
            number_for_bar[0] += 1

    return number_for_bar


def Chi2Test(self, number_for_plot, number_for_plot_nonMut):
    table = np.array([number_for_plot, number_for_plot_nonMut])
    table = 100 * table / np.sum(table)

    if 0 in np.sum(table, axis=1):
        return 1, 0

    _, p, _, expected = stats.chi2_contingency(table)
    residus = (table - expected) / np.sqrt(expected)

    return p, residus


def rearrangeZeros(self, number_for_plot, number_for_plot_nonMut):
    """S'il y a une colonne remplie de 0, on la supprime, car cela pose problème pour le test du chi2 (en plus de ne pas donner d'informations supplémentaires)"""

    L_ind = [] #liste des indices dont la colonne est remplie de 0

    for i in range(len(number_for_plot)):
        if number_for_plot[i] == 0 and number_for_plot_nonMut[i] == 0:
            L_ind.append(i)

    return L_ind


def CreateFileResFeat(self, ind_beataml, gene, feature, gene_cat, number_for_plot, number_for_plot_nonMut, p, ind_geneMut, ind_geneNonMut):

    if not os.path.exists("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes"):
        os.mkdir("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes")

    output_file = f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/ResBarplots_{gene}_{feature}/Resultats_Barplot.txt"

    if not os.path.exists(f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/ResBarplots_{gene}_{feature}"):
        os.mkdir(f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/ResBarplots_{gene}_{feature}")

    with open(output_file, 'w') as f:
        f.write(f"Résultats de barplot pour {gene} et {feature}\n")
        f.write(f"p-value = {p:.5f} selon le test du χ²\n")
        f.write(f"Résultats bruts : \n")
        Res = pd.DataFrame({
            f"{feature}": gene_cat,
            f"{gene} muté": number_for_plot,
            f"{gene} non muté": number_for_plot_nonMut
        })
        f.write(Res.to_string())
        f.write("\n")

        f.write("Liste des ID Sample par catégorie :\n")
        f.write(f"{gene} muté :\n")
        for i, cat in enumerate(gene_cat):
            if feature in ["TP53", "RUNX1", "ASXL1"]:
                if cat == "Positive":
                    ids = [ind for ind in ind_geneMut if not pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                else:
                    ids = [ind for ind in ind_geneMut if pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
            else:
                ids = [ind for ind in ind_geneMut if self.data_beat_aml.loc[ind_beataml.index(ind), feature] == cat]
            f.write(f"{cat} : {', '.join(ids)}\n")

        f.write(f"{gene} non muté :\n")
        for i, cat in enumerate(gene_cat):
            if feature in ["TP53", "RUNX1", "ASXL1"]:
                if cat == "Positive":
                    ids = [ind for ind in ind_geneNonMut if not pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                else:
                    ids = [ind for ind in ind_geneNonMut if pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
            else:
                ids = [ind for ind in ind_geneNonMut if self.data_beat_aml.loc[ind_beataml.index(ind), feature] == cat]
            f.write(f"{cat} : {', '.join(ids)}\n")

    self.fig.savefig(f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/ResBarplots_{gene}_{feature}/Barplot_plot.png")

    print(f"Les résultats bruts et la figure ont été sauvegardés dans les fichiers {output_file} et Barplot_plot.png")

    return


def plot_graph_without_abundance(self, cat_name, number_for_plot, number_for_plot_nonMut, gene, feature):
    bar_width = 0.35

    r1 = np.arange(len(cat_name))
    r2 = [x + bar_width for x in r1]

    self.fig.clear()
    self.ax = self.fig.add_subplot(111)
    self.ax.bar(r1, number_for_plot, color="blue", width=bar_width, label=f"{gene} mutated")
    self.ax.bar(r2, number_for_plot_nonMut, color="green", width=bar_width, label=f"{gene} non mutated")

    self.ax.set_xticks([r + bar_width / 2 for r in range(len(cat_name))])
    self.ax.set_xticklabels(cat_name, rotation=90)
    self.ax.set_ylabel("Nombre d'échantillons")
    self.ax.legend()
    self.ax.set_title(f"Distribution de {feature} selon la mutation (ou non) de {gene}")

    self.fig.tight_layout()
    self.canvas.draw()

    return