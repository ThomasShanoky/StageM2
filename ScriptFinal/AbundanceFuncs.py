import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns



def NormalizationByTotKmers(ind_list, data_abundance, tot_kmers_file):
        """Normalise les valeurs d'abondance par le nombre total de kmers"""

        tot_kmers = pd.read_csv(tot_kmers_file, sep=",")
        tot_kmers.set_index("Sample", inplace=True) 

        NormalizedExpression = pd.DataFrame({
            'ID Sample':ind_list,
            'Abundance': data_abundance,
        })

        for index, row in NormalizedExpression.iterrows():
            sample = row["ID Sample"]
            tot_kmers_value = tot_kmers.loc[sample]["TotalKmers"]
            NormalizedExpression.at[index, "TotalKmers"] = tot_kmers_value
            NormalizedExpression.at[index, "NormalizedExpression"] = row["Abundance"] * 10**(9) / tot_kmers_value #multiplication par 10^9 pour avoir des valeurs aux ordres de grandeur de 10 à 100
        
        return NormalizedExpression


def getFeatForPlotAbundance(data_beat_aml, NormalizedExpression, feature):
    """Récupère la valeur de la feature pour chaque échantillon dans le tableau d'expression normalisée"""
    NormalizedExpression[feature] = "" #création d'une nouvelle colonne vide

    for ind, _ in NormalizedExpression.iterrows():
        featValue = data_beat_aml.loc[ind, feature]
        NormalizedExpression.at[ind, feature] = featValue

    return NormalizedExpression


def MannWhitneyUTest(Df, feature):
    """Effectue un test de Mann-Whitney U pour deux groupes"""
    cat_name = np.unique(list(Df[feature])) #Liste de valeurs uniques de la feature (2 valeurs)
    group1 = Df[Df[feature] == cat_name[0]]["NormalizedExpression"].values #valeurs d'abondance (normalisées) pour le premier groupe
    group2 = Df[Df[feature] == cat_name[1]]["NormalizedExpression"].values

    if len(group1) > 2 and len(group2) > 2: #si on a strictement plus de 2 échantillons dans chaque groupe
        u_stat, p = stats.mannwhitneyu(group1, group2)
        return u_stat, p
    else:
        return 1, 1


def ANOVATest(Df, feature):
    """Effectue un test ANOVA pour plus de deux groupes"""
    groups = [Df[Df[feature] == cat]["NormalizedExpression"].values for cat in Df[feature].unique()] #liste de listes de valeurs d'abondance (normalisées) pour chaque groupe

    long = []
    for group in groups:
        long.append(len(group))

    if 0 in long or 1 in long or 2 in long: #s'il existe un groupe avec seulement 0, 1 ou 2 échantillons, on ne peut pas faire le test ANOVA
        return 1, 1
    
    f_stat, p = stats.f_oneway(*groups)
    return f_stat, p


def plot_graph_with_abundance(canvas, fig, NormalizedExpressionAndFeat, gene, feature, p_value):
    """Trace le graphe de l'abondance normalisée selon la feature"""

    fig.clear()
    ax = fig.add_subplot(111)

    NormalizedExpressionAndFeat = NormalizedExpressionAndFeat.reset_index(drop=True)

    sns.boxplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, ax=ax)
    sns.stripplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, color='black', alpha=0.5, jitter=True, ax=ax)

    ax.set_title(f"Distribution de l'expression normalisée selon {feature} pour le gène {gene}")
    ax.set_xlabel(feature)
    ax.set_ylabel("Expression normalisée par TBP")

    if p_value < 0.05:
        if p_value < 0.0001:
            stars = "****"
        elif p_value < 0.001:
            stars = "***"
        elif p_value < 0.01:
            stars = "**"
        else:
            stars = "*"

        x1, x2 = 0, 1
        y, h, col = NormalizedExpressionAndFeat['NormalizedExpression'].max() + 1, 1, "black"
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        ax.text((x1 + x2) * 0.5, y + h, stars, ha="center", va="bottom", color=col)


    fig.tight_layout()
    canvas.draw()

    return canvas, fig, ax


def CreateFileResAbund(fig, NormalizedExpressionAndFeat, Test, p, gene, feature, SaveAll, directory):
    """Créer le dossier résultats (avec le graphe + les résultats bruts)"""

    if not os.path.exists(f"{directory}/DossierRes"):
        os.mkdir(f"{directory}/DossierRes")

    output_file = f"{directory}/DossierRes/2_ResAbundance_{gene}_{feature}/"

    if not os.path.exists(f"{directory}/DossierRes/2_ResAbundance_{gene}_{feature}"):
        os.mkdir(f"{directory}/DossierRes/2_ResAbundance_{gene}_{feature}")
    
    with open(f"{output_file}/Resultat_Bruts.txt", 'w') as f:
        f.write(f"#Résultats bruts de l'abondance de mutation du gène {gene} selon la feature {feature}, normalisée par le nombre total de kmers par échantillon\n")
        f.write(f"#Test d'hypothèse de différence significative entre les moyennes des groupes\n#p-value = {p:.6f} avec le test {Test}\n")

    NormalizedExpressionAndFeat.to_csv(f"{output_file}/Resultat_Table.tsv", sep="\t", index=True)

    fig.savefig(f"{output_file}/Resultat_Plot.png")

    if not(SaveAll):
        print(f"Les résultats ont été sauvegardés dans le dossier {output_file}")

    return