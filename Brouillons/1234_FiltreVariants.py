#on prend le résultats de ScanR'N, on compare avec les variants trouvés dans la base BEAT AML et on filtre

import pandas as pd



Nucleotides = ["A", "T", "G", "C"]
r = 40 #nucléotide de décalage qu'on est prêt à accepter



# file_variants = "BEATAMLdata/BEATAMLvariants/Full_AML_NPM1_variants.csv"
# file_variants = "BEATAMLdata/BEATAMLvariants/Full_AML_KRAS_variants.csv"
file_variants = "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000133703_2025_03_20_11_35_05.csv"
variants = pd.read_csv(file_variants, sep=',')

# file_variantsAVerifier = "MUTdata/NPM1_alt.tsv"
file_variantsAVerifier = "MUTdata/KRAS_alt.tsv"
variantsAVerifier = pd.read_csv(file_variantsAVerifier, sep='\t')


FilteredVariants = pd.DataFrame(columns=variantsAVerifier.columns)

for _, row in variantsAVerifier.iterrows():
    pos = int(row['localisation'].split(":")[1])
    pos = [n for n in range(pos-r, pos+r+1)]
    alt = row['seq_alt']
    ref = row['seq_ref']
    nat = row['type_alt']
    found = False

    for posComp in pos:

        if nat == "Insertion" or nat == "Deletion":
            for i in range(4):
                altComp = Nucleotides[i] + alt
                refComp = Nucleotides[i] + ref
                refComp = refComp[:-1] #on enlève le tiret

                for _, row2 in variants.iterrows():
                    if row2['Reference'] == refComp and row2['Alternate'] == altComp and row2['gnomAD ID'].split("-")[1] == posComp:
                        # print("OK")
                        FilteredVariants = FilteredVariants._append(row, ignore_index=True)
                        found = True
                        break
                if found:
                    break
        else:
            for _, row2 in variants.iterrows():
                if row2['Reference'] == ref and row2['Alternate'] == alt and row2['gnomAD ID'].split("-")[1] == posComp:
                    FilteredVariants = FilteredVariants._append(row, ignore_index=True)
                    found = True
                    break

        if found:
            break
    


print(FilteredVariants)


170837543

171410543

171410539