import pandas as pd
import os



def getSamplesAndFeatures(data_beat_aml, feature, featureValues):
    samples = data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample", feature]]
    samples = samples[samples[feature].isin(featureValues)]  # ne prendre que les échantillons qui ont une valeur de feature présente dans featureValues
    ind_beataml = data_beat_aml.index.tolist() #liste d'échantillons 
    ind_feature = [] #liste des valeurs de la feature associée à l'échantillon dans ind_beataml au même indice

    for ind in ind_beataml:
        ind_feature.append(data_beat_aml.loc[ind, feature])

    SamplesAndFeatures = pd.DataFrame({"Sample": ind_beataml, "Feature": ind_feature})

    return SamplesAndFeatures


def CreateFileRes(fig, Tableau, stats_res_test, premiereLigne, fileName):
    if not os.path.exists("Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes"):
        os.mkdir("Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes")

    output_brut_file = f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/Res_{fileName}/Resultats.txt"
    if not(os.path.exists(f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/Res_{fileName}")):
        os.mkdir(f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/Res_{fileName}")

    with open(output_brut_file, "w") as f:
        f.write(premiereLigne)
        f.write(f"#Résultat du test statistique {stats_res_test}\n")
        f.write("#Tableau des résultats :\n")
        f.write(Tableau.to_string(index=False))

    output_graph_file = f"Documents/ScriptsPrincipaux/5_ScriptFinal/DossierRes/Res_{fileName}/Plot.png"
    fig.savefig(output_graph_file)

    print(f"Les résultats bruts et le graphique ont été sauvegardés respectivement dans le fichier {output_brut_file} et dans le fichier {output_graph_file}")

    return