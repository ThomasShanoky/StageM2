import numpy as np
import pandas as pd
from collections import Counter
from itertools import combinations



Genes = ['ABL1', 'ADA', 'ANKRD26', 'ASXL1', 'ASXL2', 'ATM', 'ATRX', 'BCL6', 'BCOR', 'BCORL1', 'BCR', 'BIRC3', 'BLM', 'BRAF', 'BRCA1', 'BRCA2', 'CALR', 'CARD11', 'CBL', 'CBLB', 'CDKN2A', 'CEBPA', 'CHEK2', 'CREBBP', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CUX1', 'DAXX', 'DDX41', 'DNM2', 'DNMT1', 'DNMT3A', 'EED', 'EP300', 'ETNK1', 'ETV6', 'EZH2', 'FAS', 'FBXW7', 'FLRT2', 'FLT3', 'GATA1', 'GATA2', 'GNAS', 'HNRNPK', 'HRAS', 'IDH1', 'IDH2', 'IKZF1', 'IKZF3', 'IL7R', 'JAK1', 'JAK2', 'JAK3', 'KAT6A', 'KDM6A', 'KDR', 'KIT', 'KLHDC8B', 'KLHL6', 'KMT2A', 'KMT2C', 'KRAS', 'LRRC4', 'MAP2K1', 'MLH1', 'MPL', 'MSH2', 'MSH6', 'MYC', 'MYD88', 'NBN', 'NF1', 'NOTCH1', 'NPAT', 'NPM1', 'NRAS', 'NSD1', 'NTRK3', 'P2RY2', 'PAX5', 'PDGFRA', 'PHF6', 'PML', 'PMS2', 'PRF1', 'PRPF40B', 'PRPF8', 'PTPN11', 'RAD21', 'RB1', 'RELN', 'RUNX1', 'SETBP1', 'SF1', 'SF3A1', 'SF3B1', 'SH2B3', 'SH2D1A', 'SMARCB1', 'SMC1A', 'SMC3', 'SRP72', 'SRSF2', 'STAG2', 'STAT3', 'STXBP2', 'SUZ12', 'TAL1', 'TERT', 'TET2', 'TNFRSF13B', 'TP53', 'TPMT', 'U2AF1', 'U2AF2', 'WAS', 'WRN', 'WT1', 'XPO1', 'ZRSR2']

data_mutation_files = [f"MUTdata/{gene}_alt.tsv" for gene in Genes]
data_mutation = [pd.read_csv(file, sep="\t", comment="#")[["ID_sample"]] for file in data_mutation_files]
data_mutation = [list(df["ID_sample"]) for df in data_mutation]

for i in range(len(data_mutation)):
    data_mutation[i] = [Id[:6] for Id in data_mutation[i]] #on tronque les ID pour ne garder que les 6 premiers caractères



################################
#   Nombre de gènes illimité   #
################################

def getPatientsAndMutatedGenes(data_mutation, Genes):
    PatientsAndMutatedGenes = {} #dictionnaire qui va contenir les patients en clefs et les gènes mutés en valeurs
    for ind_gene, IndListPerGene in enumerate(data_mutation):
        for IdPatient in IndListPerGene:
            if IdPatient not in PatientsAndMutatedGenes:
                PatientsAndMutatedGenes[IdPatient] = []
            PatientsAndMutatedGenes[IdPatient].append(Genes[ind_gene])
    return PatientsAndMutatedGenes


def UnicityGenes(PatientsAndMutatedGenes):
    """Parfois, des gènes se retrouvent plusieurs fois (car plusieurs fois mutés pour un même patient). On les enlève"""
    for IdPatient in PatientsAndMutatedGenes:
        PatientsAndMutatedGenes[IdPatient] = list(set(PatientsAndMutatedGenes[IdPatient]))
    return PatientsAndMutatedGenes


def getProfilMutationnel(PatientsAndMutatedGenes):
    """Renvoie un dictionnaire avec en clefs les profils mutationnels des patients et en valeurs le nombre de patients avec ce profil"""
    
    AllProfil = list(PatientsAndMutatedGenes.values())
    AllProfil = [str(sorted(geneList)) for geneList in AllProfil]
    # print(AllProfil)
    ProfilMutationnel = Counter(AllProfil)

    return dict(sorted(ProfilMutationnel.items(), key=lambda item: item[1], reverse=False))


PatientsAndMutatedGenes = getPatientsAndMutatedGenes(data_mutation, Genes)
PatientsAndMutatedGenes = UnicityGenes(PatientsAndMutatedGenes)
# print(PatientsAndMutatedGenes)
ProfilMutationnel = getProfilMutationnel(PatientsAndMutatedGenes)
# print(ProfilMutationnel)


###############################################
# Avoir la moyenne de gènes mutés par patient #
###############################################

NbGenesMutated = []
for pat in PatientsAndMutatedGenes:
    NbGenesMutated.append(len(PatientsAndMutatedGenes[pat]))

# print(f"Nombre moyen de gènes mutés par patient: {np.mean(NbGenesMutated)}")
# = 10

#############################################
# Voir les gènes les plus fréquemment mutés #
#############################################

#Attention : on ne considère qu'une seule fois un gène par patient, si un gène est muté plusieurs fois pour un patient, on ne le compte qu'une seule fois

m = 20 #voir les m gènes les plus fréquemment mutés

GenesNbMutations = {gene:0 for gene in Genes}

for pat in PatientsAndMutatedGenes:
    profil = PatientsAndMutatedGenes[pat]
    for gene in profil:
        GenesNbMutations[gene] += 1

GenesNbMutations = dict(sorted(GenesNbMutations.items(), key=lambda item: item[1], reverse=True))
GenesNbMutations = dict(list(GenesNbMutations.items())[:m])
# print(GenesNbMutations)


####################################
# Nombre de gènes limités (n fixé) #
####################################

def getCombiInOnePatient(PatientsAndMutatedGenes, IdPatient, n):
    GenesMutatedInThePatient = PatientsAndMutatedGenes[IdPatient]
    if len(GenesMutatedInThePatient) < n:
        return []
    
    Allcombi = []
    for combi in combinations(GenesMutatedInThePatient, n):
        Allcombi.append(combi)

    return Allcombi



n = 3 #nombre de gène muté à considérer
m = 20 #voir les m combinaisons les plus fréquentes

MutatedCombinations = {}
for pat in PatientsAndMutatedGenes:
    CombiOnePatient = getCombiInOnePatient(PatientsAndMutatedGenes, pat, n)
    for combi in CombiOnePatient:
        if combi not in MutatedCombinations:
            MutatedCombinations[combi] = 0
        MutatedCombinations[combi] += 1

MutatedCombinations = dict(sorted(MutatedCombinations.items(), key=lambda item: item[1], reverse=True))
MutatedCombinations = dict(list(MutatedCombinations.items())[:m])
print(MutatedCombinations)




