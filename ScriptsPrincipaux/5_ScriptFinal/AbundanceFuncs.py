import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns



def NormalizationByTotKmers(ind_list, data_abundance):

        tot_kmers_file = "Documents/ScriptsPrincipaux/TotalKmersPerSample.csv"
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
            NormalizedExpression.at[index, "NormalizedExpression"] = row["Abundance"] * 10**(9) / tot_kmers_value
        
        return NormalizedExpression


def getFeatForPlotAbundance(data_beat_aml, NormalizedExpression, feature, ind_beataml):
    NormalizedExpression[feature] = ""

    for ind, _ in NormalizedExpression.iterrows():
        print(ind)
        featValue = data_beat_aml.loc[ind_beataml.index(ind)][feature]
        NormalizedExpression.at[ind, feature] = featValue

    return NormalizedExpression


def MannWhitneyUTest(Df, feature):
    cat_name = np.unique(list(Df[feature]))
    group1 = Df[Df[feature] == cat_name[0]]["NormalizedExpression"].values
    group2 = Df[Df[feature] == cat_name[1]]["NormalizedExpression"].values

    if len(group1) > 2 and len(group2) > 2:
        _, p = stats.mannwhitneyu(group1, group2)
        return p
    else:
        return 1


def ANOVATest(Df, feature):
    groups = [Df[Df[feature] == cat]["NormalizedExpression"].values for cat in Df[feature].unique()]

    long = []
    for group in groups:
        long.append(len(group))

    if 0 in long or 1 in long or 2 in long:
        return 1
    _, p = stats.f_oneway(*groups)
    return p


def plot_graph_with_abundance(canvas, fig, NormalizedExpressionAndFeat, gene, feature):

    fig.clear()
    ax = fig.add_subplot(111)

    NormalizedExpressionAndFeat = NormalizedExpressionAndFeat.reset_index(drop=True)

    sns.boxplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, ax=ax) #showfliers=False
    sns.stripplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, color='black', size=3, jitter=True, ax=ax)

    ax.set_title(f"Distribution de l'expression normalisée selon {feature} pour le gène {gene}")
    ax.set_xlabel(feature)
    ax.set_ylabel("Expression normalisée par TBP")

    fig.tight_layout()
    canvas.draw()

    return canvas, fig, ax


def CreateFileResAbund(fig, NormalizedExpressionAndFeat, p, gene, feature):
    if not os.path.exists("Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes"):
        os.mkdir("Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes")

    output_file = f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/ResAbundance_{gene}_{feature}/Resultats_Abundance.txt"

    if not os.path.exists(f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/ResAbundance_{gene}_{feature}"):
        os.mkdir(f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/ResAbundance_{gene}_{feature}")

    if len(np.unique(NormalizedExpressionAndFeat[feature])) == 2:
        Test = "Mann-Whitney U"
    elif len(np.unique(NormalizedExpressionAndFeat[feature])) > 2:
        Test = "ANOVA"
    
    with open(output_file, 'w') as f:
        f.write(f"#Résultats bruts de l'abondance de mutation du gène {gene} selon la feature {feature}, normalisée par le gène de ménage TBP.\n")
        f.write(f"#Test d'hypothèse de différence significative entre les moyennes des groupes\n#p-value = {p:.6f} avec le test {Test}\n")
        f.write(NormalizedExpressionAndFeat.to_string(index=False))

    fig.savefig(f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/ResAbundance_{gene}_{feature}/Abundance_plot.png")

    print(f"Les résultats bruts ont été sauvegardés dans le fichier {output_file}")

    return