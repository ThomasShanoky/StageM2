import os
import numpy as np
import pandas as pd
import scipy.stats as stats



def get_number_for_bar(data_beat_aml, cat_name, ind_beataml, ind_list, feature):
    """Etant donné l'état du gène (muté ou non muté), renvoie en liste, le nombre d'échantillons pour chaque catégorie de la feature"""
    number_for_bar = [0 for _ in range(len(cat_name))]
    
    for ind in ind_list:
        index_ind = ind_beataml.index(ind)
        OneCat = data_beat_aml.iloc[index_ind][feature] # On récupère la catégorie de la feature pour cet échantillon

        for i, cat in enumerate(cat_name):
            if OneCat == cat:
                number_for_bar[i] += 1

    return number_for_bar


def Chi2Test(number_for_plot, number_for_plot_nonMut):
    """Effectue un test du chi2, renvoyant la p-value"""
    table = np.array([number_for_plot, number_for_plot_nonMut]) #Table de contingence

    if 0 in np.sum(table, axis=1): #si une ligne est remplie de 0, on ne peut pas faire le test du chi2, on fixe la p-value à 1 (non-significatif)
        return 0, 1

    chi2stat, p, _, expected = stats.chi2_contingency(table)
    # residus = (table - expected) / np.sqrt(expected)

    return chi2stat, p


def rearrangeZeros(number_for_plot, number_for_plot_nonMut):
    """S'il y a une colonne remplie de 0, on la supprime, car cela pose problème pour le test du chi2 (en plus de ne pas donner d'informations supplémentaires)"""

    L_ind = [] #liste des indices dont la colonne est remplie de 0

    for i in range(len(number_for_plot)):
        if number_for_plot[i] == 0 and number_for_plot_nonMut[i] == 0:
            L_ind.append(i)

    return L_ind


def CreateFileResFeat(data_beat_aml, gene, feature, gene_cat, number_for_plot, number_for_plot_nonMut, p, ind_geneMut, ind_geneNonMut, fig, SaveAll, directory):
    """Créer le dossier résultats (avec le graphe + les résultats bruts)"""

    if not os.path.exists(f"{directory}/DossierRes"):
        os.mkdir(f"{directory}/DossierRes")

    output_file = f"{directory}/DossierRes/1_ResBarplots_{gene}_{feature}"

    if not os.path.exists(f"{directory}/DossierRes/1_ResBarplots_{gene}_{feature}"):
        os.mkdir(f"{directory}/DossierRes/1_ResBarplots_{gene}_{feature}")

    with open(f"{output_file}/Resultat_Brut.txt", 'w') as f:
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
        for cat in gene_cat:
            ids = [ind for ind in ind_geneMut if data_beat_aml.loc[ind, feature] == cat]
            f.write(f"{cat} : {', '.join(ids)}\n")

        f.write(f"{gene} non muté :\n")
        for i, cat in enumerate(gene_cat):
            ids = [ind for ind in ind_geneNonMut if data_beat_aml.loc[ind, feature] == cat]
            f.write(f"{cat} : {', '.join(ids)}\n")

    fig.savefig(f"{output_file}/Resultat_Plot.png")

    if not(SaveAll):
        print(f"Les résultats ont été sauvegardés dans le dossier {output_file}")

    return


def plot_graph_without_abundance(fig, canvas, cat_name, number_for_plot, number_for_plot_nonMut, gene, feature, p_value):
    """Trace le graphe de la distribution des échantillons selon la mutation (ou non) du gène"""
    bar_width = 0.35

    r1 = np.arange(len(cat_name))
    r2 = [x + bar_width for x in r1]

    fig.clear()
    ax = fig.add_subplot(111)
    ax.bar(r1, number_for_plot, color="blue", width=bar_width, label=f"{gene} mutated")
    ax.bar(r2, number_for_plot_nonMut, color="green", width=bar_width, label=f"{gene} non mutated")

    ax.set_xticks([r + bar_width / 2 for r in range(len(cat_name))])
    ax.set_xticklabels(cat_name, rotation=90)
    ax.set_ylabel("Nombre d'échantillons")
    ax.legend()
    ax.set_title(f"Distribution de {feature} selon la mutation (ou non) de {gene}")
    ax.set_title(f"Répartition des échantillons {gene} mutés et non {gene} mutés selon la catégorie \"{feature}\"")

    # if p_value < 0.05:
    #     if p_value < 0.0001:
    #         stars = "****"
    #     elif p_value < 0.001:
    #         stars = "***"
    #     elif p_value < 0.01:
    #         stars = "**"
    #     else:
    #         stars = "*"

        # max_height = max(max(number_for_plot), max(number_for_plot_nonMut))
        # y = max_height + 1
        # x1 = r1[0]
        # x2 = r2[-1]

        # ax.plot([x1, x1, x2, x2], [y, y + 0.5, y + 0.5, y], color="black", lw=1.5) # Ligne en forme de crochet
        # ax.text((x1 + x2) / 2, y + 0.7, stars, ha="center", va="bottom", color="black", fontsize=12) #étoiles


    fig.tight_layout()
    canvas.draw()

    return canvas, fig, ax