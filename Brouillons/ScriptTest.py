import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from pprint import pprint



###############################################################################
##### Mettre l'ID sample directement dans le fichier BEATAML_Clinique.csv #####
################################################################################

# file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv"
# output_file = "Documents/ScriptsPrincipaux/BEATAMLdata/new_BEATAML_Cliniques.csv"

# data = pd.read_csv(file, sep=",")

# data.insert(3, "ID Sample", "")
# # print(data.columns)

# for index, row in data.iterrows():
#     var1 = row["dbgap_dnaseq_sample"]
#     var2 = row["dbgap_rnaseq_sample"]
#     if var1 == "nan" or pd.isna(var1) or var1 == float('nan'):
#         sample = var2[:6]
#     else:
#         sample = var1[:6]

#     data.at[index, "ID Sample"] = sample
# # print(data.head())

# data.to_csv(output_file, sep=",", index=False)


###################################################
##### Reformater les métadonnées de Leucegene #####
###################################################

# file = "Documents/ScriptsPrincipaux/Brouillons/GSE67040_series_matrix.txt"
# output_file = "Documents/ScriptsPrincipaux/Brouillons/Leucegene_data.csv"

# Categories = []
# Valeurs = []

# with open(file, 'r') as f:
#     for i, line in enumerate(f.readlines()):
#         line_list = line.strip().split("\t")
#         if i == 0:
#             Cat = line_list[0][1:]
#             Categories.append(Cat)
#             Val = line_list[1:][0].split(" ")
#             Valeurs.append(Val)
#         else:
#             Cat = line_list[0][1:]
#             Categories.append(Cat)
#             Val = line_list[1:]
#             Valeurs.append(Val)


# print(f"Nombre de catégories : {len(Categories)}")  # Devrait être 39
# print(f"Nombre de listes dans Valeurs : {len(Valeurs)}")  # Devrait être 39
# print(f"Longueur de chaque liste dans Valeurs : {[len(v) for v in Valeurs]}")  # Devrait être 452 pour chaque liste

# data = pd.DataFrame(Valeurs, index=Categories)
# data = data.T
# # print(data)
# data.to_csv(output_file, sep=",", index=False)


#####################################################
##### Combiner les tables cliniques de Beat-AML #####
#####################################################

#Dans un second temps, rajouter les features non présentes dans le premier fichier

# file1 = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv"
# file2 = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_CliniquesShort.csv"

# df1 = pd.read_csv(file1, sep=",")
# df2 = pd.read_csv(file2, sep=",", comment='#')
# df1.set_index("ID Sample", inplace=True)
# df2.set_index("ID Sample", inplace=True)


# common_features = list(set(df1.columns) & set(df2.columns))
# common_features.sort()
# print(f"Common features: {common_features}") #Mettre les deux features d'âge (ageAtDiagnosis et AgeCategory)
# print(len(common_features))

# data_comp = pd.DataFrame()

# data_comp["ID Sample"] = df1.index
# # data_comp.set_index("ID Sample", inplace=True)
# # print(data_comp)

# for feat in common_features:
#     col1 = df1[feat]
#     col2 = df2[feat]
#     data_comp[f"{feat}_1"] = data_comp["ID Sample"].map(col1)
#     data_comp[f"{feat}_2"] = data_comp["ID Sample"].map(col2)

# data_comp.set_index("ID Sample", inplace=True)
# # print(data_comp)
# # print(len(data_comp.columns))
# # data_comp.to_csv("Documents/ScriptsPrincipaux/Brouillons/BEATAML_Compare.csv", sep=",", index=True)

# CompMismatches = {}
# for i in range(int(len(data_comp.columns)/2)):
#     c = 0
#     Diff = []
#     feat1 = data_comp.columns[2*i]
#     feat2 = data_comp.columns[2*i+1]
#     for index, row in data_comp.iterrows():
#         val1 = row[feat1]
#         val2 = row[feat2]
#         if val1 != val2:
#             c += 1
#             Diff.append((index, val1, val2))
#     CompMismatches[common_features[i]] = [c, Diff]

# # pprint(CompMismatches)

# for feat, (c, Diff) in CompMismatches.items():
#     if not("nan" in str(Diff)):
#         print(f"{feat} : {Diff}")
#         #Rien => Les différences ne sont que des valeurs manquantes entre deux versions des fichiers


# #Ajouter : isTherapy

# df1.insert(20, "isTherapy", "")

# df1["isTherapy"] = df1.index.map(df2["isTherapy"])


# #Ajouter : Categoriser les âges

# df1.insert(17, "AgeCategory", "")
# # print(df1.columns)

# for index, row in df1.iterrows():
#     age = row["ageAtDiagnosis"]
#     if age <= 45:
#         df1.at[index, "AgeCategory"] = "young"
#     elif age <= 60:
#         df1.at[index, "AgeCategory"] = "middle"
#     elif age <= 75:
#         df1.at[index, "AgeCategory"] = "older"
#     elif age <= 120:
#         df1.at[index, "AgeCategory"] = "oldest"

# print(df1)
# df1.to_csv("Documents/ScriptsPrincipaux/BEATAMLdata/test.csv", sep=",", index=True)


#################################################################
##### Refaire la distribution de CEBPA avec CEBPA_Biallelic #####
#################################################################

# file = "ScriptsPrincipaux/ScriptFinal/BEATAML_Cliniques.csv"
# df = pd.read_csv(file, sep=",", comment='#')[["ID Sample", "CEBPA_Biallelic"]]


# file_mut = "ScriptsPrincipaux/ScriptFinal/MUTdata/CEBPA_alt_perso.csv"
# df_mut = pd.read_csv(file_mut, sep=",", comment='#')[["sampleID"]]
# sampleMutated = df_mut["sampleID"].tolist()

# for index, row in df.iterrows():
#     samp = row["ID Sample"]
#     if samp in sampleMutated:
#         df.at[index, "CEBPAmut"] = 1
#     else:
#         df.at[index, "CEBPAmut"] = 0

# df["CEBPA_Biallelic"] = df["CEBPA_Biallelic"].fillna("NaN")
# print(df)

# plt.figure(figsize=(10,6))
# sns.countplot(
#     data=df,
#     x="CEBPA_Biallelic",
#     hue="CEBPAmut",
#     palette={0:"green", 1:"blue"}
# )

# plt.tight_layout()
# plt.show()

# contingency_table = pd.crosstab(df["CEBPA_Biallelic"], df["CEBPAmut"])
# chi2, p, dof, expected = chi2_contingency(contingency_table)
# print(f"Chi2: {chi2}, p-value: {p}, dof: {dof}")


#########################################
##### Discrétiser la survie globale #####
#########################################

# file = "ScriptsPrincipaux/ScriptFinal/BEATAML_Cliniques.csv"
# df = pd.read_csv(file, sep=",", comment='#')

# overall_survival = list(df["overallSurvival"])
# # print(overall_survival)

# # plt.hist(overall_survival, bins=80)
# # plt.show()

# for index, row in df.iterrows():
#     surv = row["overallSurvival"]
#     if surv <= 365:
#         df.at[index, "overallSurvivalDiscretized"] = "LessThanAYear"
#     else:
#         df.at[index, "overallSurvivalDiscretized"] = "MoreThanAYear"

# # print(df)
# df.to_csv("ScriptsPrincipaux/Brouillons/BEATAML_Cliniques.csv", sep=",", index=False)


############################################################
##### Trouver une ligne décrivant un échantillon donné #####
############################################################

# file = "ScriptFinal/BEATAML_Cliniques.csv"
# df = pd.read_csv(file, sep=",", comment='#')

# sample_id = "BA3216"

# pd.set_option('display.max_columns', None)
# print(df[df["ID Sample"] == sample_id])


#################################################################
##### Taux d'annotation pour les gènes ASXL1, RUNX1 et TP53 #####
#################################################################

# file = "ScriptFinal/BEATAML_Cliniques.csv"
# df = pd.read_csv(file, sep=",", comment='#')[["ID Sample", "ASXL1", "RUNX1", "TP53"]]

# # print(df)

# nb_ech = 0
# nb_asxl1_annoted = 0
# nb_runx1_annoted = 0
# nb_tp53_annoted = 0

# for index, row in df.iterrows():
#     if pd.isna(row["ASXL1"]):
#         nb_asxl1_annoted += 1
#     if pd.isna(row["RUNX1"]):
#         nb_runx1_annoted += 1
#     if pd.isna(row["TP53"]):
#         nb_tp53_annoted += 1
#     nb_ech += 1

# print(100*(nb_ech-nb_asxl1_annoted)/nb_ech)
# print(100*(nb_ech-nb_runx1_annoted)/nb_ech)
# print(100*(nb_ech-nb_tp53_annoted)/nb_ech)

# print(f"{0.3+0.5:.20f}")


######################################################################################################################
##### Vérifier que la liste des échantillons indexé et celle des échantillons possédant un RNAseq sont les mêmes #####
######################################################################################################################

file_index = "ScriptFinal/BEATAML_index.tsv"

ech_index = list(pd.read_csv(file_index, sep="\t", comment='#'))
ech_index = [ech[:6] for ech in ech_index]

file_metadata = "ScriptFinal/BEATAML_Cliniques2.csv"
data = pd.read_csv(file_metadata, sep=",", comment='#')
ech_rnaseq = list(data["dbgap_rnaseq_sample"])
ech_rnaseq = [ech[:6] for ech in ech_rnaseq if not pd.isna(ech)]


c = 0
for ech in ech_index:
    if ech not in ech_rnaseq:
        c += 1

print(c) # c = 0 + les deux listes sont de même longueurs => Les deux listes sont les mêmes


data = data.dropna(subset=["dbgap_rnaseq_sample"])
print(data)
print(len(data)) 
data.to_csv("ScriptFinal/BEATAML_Cliniques.csv", sep=",", index=False)
