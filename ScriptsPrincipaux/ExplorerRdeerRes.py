# à partir d'un gène et d'un individu, trouve les résultats bruts de rdeer
from pprint import pprint
import pandas as pd
import os


gene = "NPM1"
# samples = ["BA2525", "BA2069", "BA2631", "BA2455", "BA2377", "BA2036", "BA2970", "BA2024", 
#  "BA2521", "BA2400", "BA2762", "BA3076", "BA2228", "BA2656", "BA2764", "BA2456", 
#  "BA2218", "BA3109", "BA3173", "BA3436", "BA3162", "BA3169", "BA3255", "BA3335", 
#  "BA3227", "BA3376"]
samples = ["BA2879", "BA3070", "BA3216"]

variants_file = f"Documents/RecreateData/Full_AML_GeneVar/Full_AML_{gene}_variants.csv"
Variants_gene = pd.read_csv(variants_file)[["in_cosmic"]]
# print(Variants_gene.iloc[5])

files_list = os.listdir("Documents/RecreateData/RdeerOutput/")
files_list = ["Documents/RecreateData/RdeerOutput/"+file for file in files_list if gene in file]

def extract_variant_number(file_name):
    # Extraire le numéro entre "var" et "_output"
    start = file_name.find("_var") + 4
    end = file_name.find("_output")
    return int(file_name[start:end])

# Trier la liste en fonction du numéro du variant
files_list_sorted = sorted(files_list, key=extract_variant_number)
pprint(files_list_sorted)
# print(int(len(files_list)/2))


AllDF = []
for i in range(int(len(files_list_sorted)/2)):
    filePart1 = files_list_sorted[2*i]
    filePart2 = files_list_sorted[2*i+1]
    df1 = pd.read_csv(filePart1, sep="\t")
    df2 = pd.read_csv(filePart2, sep="\t")

    df1.set_index(df1.columns[0], inplace=True)
    df2.set_index(df2.columns[0], inplace=True)

    df = pd.concat([df1, df2], axis=1)
    df.columns = [col[:6] for col in df.columns]
    AllDF.append(df)

for sample in samples:
    maxi = 0
    col_maxi = None
    long = AllDF[0].shape[0]
    in_cos = None
    var = None
    df_t = None
    for i, df in enumerate(AllDF):
        # print(df[sample])
        somme = sum(df[sample])
        if somme > maxi:
            maxi = somme
            col_maxi = df[sample]
            long = df.shape[0]
            in_cos = Variants_gene.iloc[i]["in_cosmic"]
            var = i
            df_t = df

    print("\n")
    print(f"Somme des kmers pour {sample} : {maxi}")
    seuil = 10*long
    print(f"Seuil : {seuil}")
    print(f"Variant : {var}")
    # print(f"DF : {df_t}")
    print(f"Dans cosmic ? : {in_cos}")



