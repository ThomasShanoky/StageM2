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

# file_index = "ScriptFinal/BEATAML_index.tsv"

# ech_index = list(pd.read_csv(file_index, sep="\t", comment='#'))
# ech_index = [ech[:6] for ech in ech_index]

# file_metadata = "ScriptFinal/BEATAML_Cliniques2.csv"
# data = pd.read_csv(file_metadata, sep=",", comment='#')
# ech_rnaseq = list(data["dbgap_rnaseq_sample"])
# ech_rnaseq = [ech[:6] for ech in ech_rnaseq if not pd.isna(ech)]


# c = 0
# for ech in ech_index:
#     if ech not in ech_rnaseq:
#         c += 1

# print(c) # c = 0 + les deux listes sont de même longueurs => Les deux listes sont les mêmes


# data = data.dropna(subset=["dbgap_rnaseq_sample"])
# print(data)
# print(len(data)) 
# data.to_csv("ScriptFinal/BEATAML_Cliniques.csv", sep=",", index=False)


############################################################################
##### Normaliser les données d'expression obtenues par approches k-mer #####
############################################################################

# file_totkmer = "ScriptFinal/TotalKmersPerSample.csv"
# file_expr = "ScriptFinal/ExpressionsWithKmers.csv"

# totkmers_dict = {}
# with open(file_totkmer, 'r') as f:
#     for i, line in enumerate(f.readlines()):
#         if i != 0:
#             line_list = line.strip().split(",")
#             samp = line_list[0]
#             tot_kmers = int(line_list[1])
#             totkmers_dict[samp] = tot_kmers

# # print(totkmers_dict)

# data = pd.read_csv(file_expr, sep=",", comment='#')

# for col in data.columns:
#     if col in totkmers_dict:
#         data[col] = (data[col] * 10**9) / totkmers_dict[col]


# data.to_csv("ScriptFinal/ExpressionsWithKmers2.csv", sep=",", index=False)


##############################################################
##### Avoir la liste des 141 gènes impliquées dans l'AML #####
##############################################################

# file = "ScriptFinal/data_mutations.vcf"
# df = pd.read_csv(file, sep="\t", comment='#')

# genes_list = list(np.unique(df["GENE"]))
# print(len(genes_list)) #125
# print(genes_list)


# import pandas as pd


# df[['CHROM', 'POS']] = df['localisation'].str.extract(r'(Chr[\w]+):(\d+)', expand=True)
# df['POS'] = df['POS'].astype(int)
# df['ID'] = df['ID_sample'].astype(str).str[:6]
# df['GENE'] = df['Gene']
# df['REF'] = df['seq_ref']
# df['ALT'] = df['seq_alt']
# df['ABUNDANCE'] = df['count_alt'].astype(str)
# df['RATIO'] = df['% ratio'].str.replace('%', '', regex=False).astype(float)
# df['QUAL'] = '.'
# df['FILTER'] = 'PASS'

# df['INFO'] = (
#     'Gene=' + df['Gene'] +
#     ';Ensembl_ID=' + df['Ensembl_ID'] +
#     ';type_alt=' + df['type_alt'] +
#     ';abundance=' + df['ABUNDANCE'].astype(str) +
#     ';ratio=' + df['RATIO'].astype(str) +
#     ';count_ref=' + df['count_ref'].astype(str) +
#     ';count_alt=' + df['count_alt'].astype(str)
# )

# vcf_df = df[['CHROM', 'POS', 'ID', 'GENE', 'REF', 'ALT', 'ABUNDANCE', 'RATIO', 'QUAL', 'FILTER', 'INFO']]
# vcf_df.columns = ['CHROM', 'POS', 'ID', 'GENE', 'REF', 'ALT', 'ABUNDANCE', 'RATIO', 'QUAL', 'FILTER', 'INFO']
# vcf_df.to_csv("output.vcf", sep="\t", index=False)


#######################################
##### Formatage table expressions #####
#######################################

# import pandas as pd


# file = "beataml_waves1to4_norm_exp_dbgap.txt"
# df = pd.read_csv(file, sep="\t", comment='#')

# df = df.drop(["stable_id", "description", "biotype"], axis=1)
# df.columns = [col[:6] for col in df.columns]

# df.to_csv("Expressions.csv", sep=",", index=False)



############################################################
##### Formatage nombre total de k-mers par échantillon #####
############################################################

# import pandas as pd


# file = "kmerstot_per_sample.txt"
# df = pd.DataFrame(columns=["Sample", "TotalKmers"])

# with open(file, 'r') as f:
#     for line in f.readlines():
#         line_list = line.split(":")
#         sample = line_list[0][:6]
#         tot_kmers = line_list[1].strip()
#         df = df._append({"Sample": sample, "TotalKmers": tot_kmers}, ignore_index=True)

# df.to_csv("TotalKmersPerSample.csv", sep=",", index=False)



import pandas as pd
import os

# file = "Elmer/ExpressionKmers/CountKmers.tsv"
# df = pd.read_csv(file, sep="\t", comment='#')

# df.columns = [col[:6] for col in df.columns]

# df.to_csv("Elmer/ExpressionKmers/CountKmersFibn.tsv", sep="\t", index=False)



# data_beat_aml_file = f"ScriptFinal/Metadata.csv"
# data_beat_aml = pd.read_csv(data_beat_aml_file, sep=",", comment="#")
# data_beat_aml = data_beat_aml[data_beat_aml["diseaseStageAtSpecimenCollection"] == "Initial Diagnosis"] 
# data_beat_aml.to_csv(f"ScriptFinal/MetadataInitialDiagnosis.csv", sep=",", index=False)



directory = '/'.join(os.path.abspath(__file__).split('/')[:-1])


# vcf_file = f"{directory}/data_mutations.vcf"
vcf_file = f"ScriptFinal/gene_test_thomas.vcf"

#Pré traitement VCF

with open(vcf_file, "r") as f:
    for i, line in enumerate(f):
        if line.startswith("#CHROM"):
            header_line = i 
            break

vcf_df = pd.read_csv(vcf_file, sep="\t", skiprows=header_line)
vcf_df.columns = vcf_df.columns.str.replace("#", "")  
# print(vcf_df.columns)


vcf_df["GENE"] = ""
vcf_df["SAMPLE"] = ""
vcf_df["ABUNDANCE"] = 0.0

new_rows = []

for index, row in vcf_df.iterrows():
    gene = row["INFO"].split("GENE=")[1].split(";")[0]
    nb_sample = int(row["INFO"].split("NUMBER_SAMPLES=")[1].split(";")[0])
    all_samples = row["INFO"].split("SAMPLES=")[2]
    for i in range(nb_sample):
        sample = all_samples.split("(")[1+i].split(",")[0].strip("'\"")[:6]
        abundance = float(all_samples.split("(")[1+i].split(",")[1].strip(")"))

        if i == 0: #réécrire la ligne
            vcf_df.at[index, "GENE"] = gene
            vcf_df.at[index, "SAMPLE"] = sample
            vcf_df.at[index, "ABUNDANCE"] = abundance

        else: #ajouter une ligne
            new_row = row.to_dict()
            new_row["GENE"] = gene
            new_row["SAMPLE"] = sample
            new_row["ABUNDANCE"] = abundance
            new_rows.append(new_row)

if new_rows:
    vcf_df = pd.concat([vcf_df, pd.DataFrame(new_rows)], ignore_index=True)

# print(vcf_df)
vcf_df.to_csv(f"{directory}/data_mutations_processed.csv", sep=",", index=False)
