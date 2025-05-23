import pandas as pd



main_file = "ScriptFinal/BEATAML_Cliniques.csv"
secc_file = "Brouillons/ClinicalSummary.csv"

main_df = pd.read_csv(main_file, comment='#')
secc_df = pd.read_csv(secc_file, comment='#')
secc_df.set_index("ID Sample", inplace=True)



features_main = list(main_df.columns)
features_secc = list(secc_df.columns)

print(len(features_main))
print(len(features_secc))

FeaturesNotFound = []
for feat in features_secc:
    if feat not in features_main:
        FeaturesNotFound.append(feat)

print(FeaturesNotFound)
print(len(FeaturesNotFound)) #nombre de colonne à ajouter


df_main_copy = main_df.copy()
df_main_copy.set_index("ID Sample", inplace=True)

for feat in FeaturesNotFound:
    df_main_copy[feat] = ""
    for index, row in secc_df.iterrows():
        val = row[feat]
        df_main_copy.at[index, feat] = val

print(df_main_copy)
df_main_copy.to_csv("ScriptFinal/BEATAML_Cliniques2.csv", index=True)
        