import numpy as np
import pandas as pd
from pprint import pprint



##### Mettre l'ID sample directement dans le fichier BEATAML_Clinique.csv

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


##### Reformater les métadonnées de Leucegene

file = "Documents/ScriptsPrincipaux/Brouillons/GSE67040_series_matrix.txt"
output_file = "Documents/ScriptsPrincipaux/Brouillons/Leucegene_data.csv"

Categories = []
Valeurs = []

with open(file, 'r') as f:
    for i, line in enumerate(f.readlines()):
        line_list = line.strip().split("\t")
        if i == 0:
            Cat = line_list[0][1:]
            Categories.append(Cat)
            Val = line_list[1:][0].split(" ")
            Valeurs.append(Val)
        else:
            Cat = line_list[0][1:]
            Categories.append(Cat)
            Val = line_list[1:]
            Valeurs.append(Val)


print(f"Nombre de catégories : {len(Categories)}")  # Devrait être 39
print(f"Nombre de listes dans Valeurs : {len(Valeurs)}")  # Devrait être 39
print(f"Longueur de chaque liste dans Valeurs : {[len(v) for v in Valeurs]}")  # Devrait être 452 pour chaque liste

data = pd.DataFrame(Valeurs, index=Categories)
data = data.T
# print(data)
data.to_csv(output_file, sep=",", index=False)