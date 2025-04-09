import pandas as pd



# Fusioner les fichiers

S_file = "SRSF1.xlsx"
S_Table = pd.ExcelFile(S_file)
SRSF1_Table = S_Table.parse("Feuil1")
# print(SRSF1_Table)

M_file = "MNK.xlsx"
M_Table = pd.ExcelFile(M_file)
MNK_Table = M_Table.parse("Feuil1")
# print(MNK_Table)

new_Table = pd.DataFrame({"DiseaseStageAtSpecimenCollection":[], "Individus":[], "SRSF1":[], "ATF4":[], "MNK1-E11E12":[], "MNK1-E11E13":[], "MNK2-E13E14a":[], "MNK2-E13E14b":[]})
# print(new_Table)

for _, ligneS in SRSF1_Table.iterrows():
    IndS = ligneS["Individus"]
    for _, ligneM in MNK_Table.iterrows():
        indM = ligneM["seq_name"]
        if indM == IndS:
            new_Table = new_Table._append({"DiseaseStageAtSpecimenCollection":ligneS["DiseaseStageAtSpecimenCollection"], "Individus":ligneS["Individus"], "SRSF1":ligneS["SRSF1"], "ATF4":ligneS["ATF4"], "MNK1-E11E12":ligneM["MNK1-E11E12"], "MNK1-E11E13":ligneM["MNK1-E11E13"], "MNK2-E13E14a":ligneM["MNK2-E13E14a"], "MNK2-E13E14b":ligneM["MNK2-E13E14b"]}, ignore_index=True)

print(new_Table)
new_Table.to_excel("SRSF1_MNK.xlsx", index=False)



# Ajouter 'BA' aux samples

file = "Documents/Classeur1.xlsx"
Table = pd.ExcelFile(file)
Table = Table.parse("Feuil1")
Table["Individus"] = "BA" + Table["Individus"].astype(str)
Table.to_excel("Documents/Classeur1_new.xlsx", index=False)


