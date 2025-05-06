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


def CreateFileRes(directory, fig, Tableau, stats_res_test, premiereLigne, fileName, SaveAll):
    if not os.path.exists(f"{directory}/DossierRes"):
        os.mkdir(f"{directory}/DossierRes")

    output_file = f"{directory}/DossierRes/Res_{fileName}"
    if not(os.path.exists(f"{directory}/DossierRes/Res_{fileName}")):
        os.mkdir(f"{directory}/DossierRes/Res_{fileName}")

    with open(f"{output_file}/Resultat_Brut.txt", "w") as f:
        f.write(premiereLigne)
        f.write(f"#Résultat du test statistique {stats_res_test}\n")

    Tableau.to_csv(f"{output_file}/Resultat_Table.tsv", sep="\t", index=False)
    fig.savefig(f"{output_file}/Resultat_Plot.png")

    if not(SaveAll):
        print(f"Les résultats ont été sauvegardés dans le dossier {output_file}")

    return