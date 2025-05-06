import pandas as pd
from pprint import pprint


NuclPossible = ["A", "T", "C", "G"]

r = 25 #erreur sur la position qu'on est prêt à accepter


# Etant donnés la nature de la mutation (insertion, substitution ...), sa position, sa séquence référente et sa séquence alternative, on peut identifier les mutations trouvées par ScanR'N


def readFiles(Gene, MutationAnnotatedFile):
    Gene_file = f"MUTdata/{Gene}_alt.tsv"
    Gene_mut = pd.read_csv(Gene_file, sep="\t", comment="#")[["ID_sample", "localisation", "type_alt", "seq_ref", "seq_alt"]]
    MutationAnnotated = pd.read_csv(MutationAnnotatedFile, sep=",")[["gnomAD ID", "Position", "rsIDs", "Reference", "Alternate", "HGVS Consequence", "Protein Consequence", "VEP Annotation"]]
    return Gene_mut, MutationAnnotated


def AnnotateMutation(NatureVar, SeqRefVarInit, SeqAltVarInit, RangePosVar, MutationAnnotated):
    """Pour une mutation trouvée par ScanR'N, retourne la ligne correspondante dans le fichier de mutation annotée"""
    VariantsRetenus = []
    for _, row in MutationAnnotated.iterrows():
        pos = row["Position"]
        SeqRef = row["Reference"]
        SeqAlt = row["Alternate"]

        if NatureVar == "Mutation":
            SeqRefVar = SeqRefVarInit
            SeqAltVar = SeqAltVarInit
            if pos in RangePosVar and SeqRef == SeqRefVar and SeqAlt == SeqAltVar:
                VariantsRetenus.append(row)
        else:
            for nucl in NuclPossible: #dans ScanR'N, la lettre de départ de la mutation n'est pas écrite, mais dans le fichier de mutation si
                SeqRefVar = nucl + SeqRefVarInit[:-1]
                SeqAltVar = nucl + SeqAltVarInit
                if pos in RangePosVar and SeqRef == SeqRefVar and SeqAlt == SeqAltVar:
                    VariantsRetenus.append(row)
    return VariantsRetenus


def AnnotateAllMutations(Gene_mut, MutationAnnotated):
    Variants = []
    for _, row in Gene_mut.iterrows():
        NatureVar = row["type_alt"]
        SeqRefVarInit = row["seq_ref"]
        SeqAltVarInit = row["seq_alt"]
        posVarInit = row["localisation"]
        posVarInit = int(posVarInit.split(":")[1])

        RangePosVarInf = [i for i in range(posVarInit-r, posVarInit)]
        RangePosVarSup = [i for i in range(posVarInit, posVarInit+r)]
        RangePosVar = RangePosVarInf + RangePosVarSup

        VariantsRetenus = AnnotateMutation(NatureVar, SeqRefVarInit, SeqAltVarInit, RangePosVar, MutationAnnotated)
        Variants.append(VariantsRetenus)
    return Variants


def WhereAreAnnotatedMutation(Gene_var):
    new_NPM1_var = []
    for i, lst in enumerate(Gene_var):
        if lst == []:
            new_NPM1_var.append(0)
        else:
            new_NPM1_var.append(i+2) #+2 car header + liste en python commence à 0
    return new_NPM1_var


def PercentOfAnnotatedMutations(Gene_var):
    per = 0
    for n in Gene_var:
        if n != 0:
            per += 1
    return per/len(Gene_var)



########## NPM1 ##########

# NPM1_mut, MutationAnnotatedNPM1 = readFiles("NPM1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000181163_2025_03_14_15_51_31.csv")

# NPM1_var = AnnotateAllMutations(NPM1_mut, MutationAnnotatedNPM1)
# new_NPM1_var = WhereAreAnnotatedMutation(NPM1_var)
# print(f"{PercentOfAnnotatedMutations(new_NPM1_var)*100:.2f}% des mutations trouvées de NPM1 ont été annotées")

# print(new_NPM1_var)

# AnnotatedMutations_dico = {}
# for nMut in new_NPM1_var:
#     if nMut != 0:
#         AnnotatedMutations_dico[nMut] = [NPM1_var[nMut-2][0]["gnomAD ID"], 
#                                          NPM1_var[nMut-2][0]["Position"],
#                                          NPM1_var[nMut-2][0]["rsIDs"],
#                                          NPM1_var[nMut-2][0]["Reference"], 
#                                          NPM1_var[nMut-2][0]["Alternate"],
#                                          NPM1_var[nMut-2][0]["HGVS Consequence"], 
#                                          NPM1_var[nMut-2][0]["Protein Consequence"],
#                                          NPM1_var[nMut-2][0]["VEP Annotation"]]
# pprint(AnnotatedMutations_dico)



rnaseq = "BEATAMLdata/BEATAML_Cliniques.csv"

n = 0
df = pd.read_csv(rnaseq, sep=",")[["dbgap_rnaseq_sample"]]
# print(df)
for _, row in df.iterrows():
    samp = row["dbgap_rnaseq_sample"]
    if pd.isna(samp):
        n += 1
print(942-n)


########## TP53 ##########

TP53_mut, MutationAnnotatedTP53 = readFiles("TP53", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000141510_2025_03_20_09_25_55.csv")

TP53_var = AnnotateAllMutations(TP53_mut, MutationAnnotatedTP53)
new_TP53_var = WhereAreAnnotatedMutation(TP53_var)
print(f"{PercentOfAnnotatedMutations(new_TP53_var)*100:.2f}% des mutations trouvées de TP53 ont été annotées")


########## DNMT3A ##########

DNMT3A_mut, MutationAnnotatedDNMT3A = readFiles("DNMT3A", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000119772_2025_03_20_09_32_34.csv")

DNMT3A_var = AnnotateAllMutations(DNMT3A_mut, MutationAnnotatedDNMT3A)
new_DNMT3A_var = WhereAreAnnotatedMutation(DNMT3A_var)
print(f"{PercentOfAnnotatedMutations(new_DNMT3A_var)*100:.2f}% des mutations trouvées de DNMT3A ont été annotées")


# ########## FLT3 ##########

# FLT3_mut, MutationAnnotatedFLT3 = readFiles("FLT3", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000122025_2025_03_20_09_33_57.csv")

# FLT3_var = AnnotateAllMutations(FLT3_mut, MutationAnnotatedFLT3)
# new_FLT3_var = WhereAreAnnotatedMutation(FLT3_var)
# print(f"{PercentOfAnnotatedMutations(new_FLT3_var)*100:.2f}% des mutations trouvées de FLT3 ont été annotées")


# ########## TET2 ##########

# TET2_mut, MutationAnnotatedTET2 = readFiles("TET2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000168769_2025_03_20_11_28_48.csv")

# TET2_var = AnnotateAllMutations(TET2_mut, MutationAnnotatedTET2)
# new_TET2_var = WhereAreAnnotatedMutation(TET2_var)
# print(f"{PercentOfAnnotatedMutations(new_TET2_var)*100:.2f}% des mutations trouvées de TET2 ont été annotées")


# ########## NRAS ##########

# NRAS_mut, MutationAnnotatedNRAS = readFiles("NRAS", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000213281_2025_03_20_11_30_36.csv")

# NRAS_var = AnnotateAllMutations(NRAS_mut, MutationAnnotatedNRAS)
# new_NRAS_var = WhereAreAnnotatedMutation(NRAS_var)
# print(f"{PercentOfAnnotatedMutations(new_NRAS_var)*100:.2f}% des mutations trouvées de NRAS ont été annotées")


# ########## RUNX1 ##########

# RUNX1_mut, MutationAnnotatedRUNX1 = readFiles("RUNX1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000159216_2025_03_20_11_31_41.csv")

# RUNX1_var = AnnotateAllMutations(RUNX1_mut, MutationAnnotatedRUNX1)
# new_RUNX1_var = WhereAreAnnotatedMutation(RUNX1_var)
# print(f"{PercentOfAnnotatedMutations(new_RUNX1_var)*100:.2f}% des mutations trouvées de RUNX1 ont été annotées")


# ########## IDH2 ##########

# IDH2_mut, MutationAnnotatedIDH2 = readFiles("IDH2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000182054_2025_03_20_11_32_29.csv")

# IDH2_var = AnnotateAllMutations(IDH2_mut, MutationAnnotatedIDH2)
# new_IDH2_var = WhereAreAnnotatedMutation(IDH2_var)
# print(f"{PercentOfAnnotatedMutations(new_IDH2_var)*100:.2f}% des mutations trouvées de IDH2 ont été annotées")


# ########## ASXL1 ##########

# ASXL1_mut, MutationAnnotatedASXL1 = readFiles("ASXL1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000171456_2025_03_20_11_33_17.csv")

# ASXL1_var = AnnotateAllMutations(ASXL1_mut, MutationAnnotatedASXL1)
# new_ASXL1_var = WhereAreAnnotatedMutation(ASXL1_var)
# print(f"{PercentOfAnnotatedMutations(new_ASXL1_var)*100:.2f}% des mutations trouvées de ASXL1 ont été annotées")


# ########## WT1 ##########

# WT1_mut, MutationAnnotatedWT1 = readFiles("WT1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000184937_2025_03_20_11_34_13.csv")

# WT1_var = AnnotateAllMutations(WT1_mut, MutationAnnotatedWT1)
# new_WT1_var = WhereAreAnnotatedMutation(WT1_var)
# print(f"{PercentOfAnnotatedMutations(new_WT1_var)*100:.2f}% des mutations trouvées de WT1 ont été annotées")


# ########## KRAS ##########

# KRAS_mut, MutationAnnotatedKRAS = readFiles("KRAS", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000133703_2025_03_20_11_35_05.csv")

# KRAS_var = AnnotateAllMutations(KRAS_mut, MutationAnnotatedKRAS)
# new_KRAS_var = WhereAreAnnotatedMutation(KRAS_var)
# print(f"{PercentOfAnnotatedMutations(new_KRAS_var)*100:.2f}% des mutations trouvées de KRAS ont été annotées")


# ########## IDH1 ##########

# IDH1_mut, MutationAnnotatedIDH1 = readFiles("IDH1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000138413_2025_03_20_11_35_53.csv")

# IDH1_var = AnnotateAllMutations(IDH1_mut, MutationAnnotatedIDH1)
# new_IDH1_var = WhereAreAnnotatedMutation(IDH1_var)
# print(f"{PercentOfAnnotatedMutations(new_IDH1_var)*100:.2f}% des mutations trouvées de IDH1 ont été annotées")


# ########## PTPN11 ##########

# PTPN11_mut, MutationAnnotatedPTPN11 = readFiles("PTPN11", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000179295_2025_03_20_11_37_00.csv")

# PTPN11_var = AnnotateAllMutations(PTPN11_mut, MutationAnnotatedPTPN11)
# new_PTPN11_var = WhereAreAnnotatedMutation(PTPN11_var)
# print(f"{PercentOfAnnotatedMutations(new_PTPN11_var)*100:.2f}% des mutations trouvées de PTPN11 ont été annotées")


# ########## SRSF2 ##########

# SRSF2_mut, MutationAnnotatedSRSF2 = readFiles("SRSF2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000161547_2025_03_20_11_37_49.csv")

# SRSF2_var = AnnotateAllMutations(SRSF2_mut, MutationAnnotatedSRSF2)
# new_SRSF2_var = WhereAreAnnotatedMutation(SRSF2_var)
# print(f"{PercentOfAnnotatedMutations(new_SRSF2_var)*100:.2f}% des mutations trouvées de SRSF2 ont été annotées")


# ########## CEBPA ##########

# CEBPA_mut, MutationAnnotatedCEBPA = readFiles("CEBPA", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000245848_2025_03_20_11_39_05.csv")

# CEBPA_var = AnnotateAllMutations(CEBPA_mut, MutationAnnotatedCEBPA)
# new_CEBPA_var = WhereAreAnnotatedMutation(CEBPA_var)
# print(f"{PercentOfAnnotatedMutations(new_CEBPA_var)*100:.2f}% des mutations trouvées de CEBPA ont été annotées")


# ########## KIT ##########

# KIT_mut, MutationAnnotatedKIT = readFiles("KIT", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000157404_2025_03_20_11_39_52.csv")

# KIT_var = AnnotateAllMutations(KIT_mut, MutationAnnotatedKIT)
# new_KIT_var = WhereAreAnnotatedMutation(KIT_var)
# print(f"{PercentOfAnnotatedMutations(new_KIT_var)*100:.2f}% des mutations trouvées de KIT ont été annotées")


# ########## NF1 ##########

# NF1_mut, MutationAnnotatedNF1 = readFiles("NF1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000196712_2025_03_20_11_40_40.csv")

# NF1_var = AnnotateAllMutations(NF1_mut, MutationAnnotatedNF1)
# new_NF1_var = WhereAreAnnotatedMutation(NF1_var)
# print(f"{PercentOfAnnotatedMutations(new_NF1_var)*100:.2f}% des mutations trouvées de NF1 ont été annotées")


# ########## STAG2 ##########

# STAG2_mut, MutationAnnotatedSTAG2 = readFiles("STAG2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000101972_2025_03_20_11_41_28.csv")

# STAG2_var = AnnotateAllMutations(STAG2_mut, MutationAnnotatedSTAG2)
# new_STAG2_var = WhereAreAnnotatedMutation(STAG2_var)
# print(f"{PercentOfAnnotatedMutations(new_STAG2_var)*100:.2f}% des mutations trouvées de STAG2 ont été annotées")


# ########## GATA2 ##########

# GATA2_mut, MutationAnnotatedGATA2 = readFiles("GATA2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000179348_2025_03_20_11_45_47.csv")

# GATA2_var = AnnotateAllMutations(GATA2_mut, MutationAnnotatedGATA2)
# new_GATA2_var = WhereAreAnnotatedMutation(GATA2_var)
# print(f"{PercentOfAnnotatedMutations(new_GATA2_var)*100:.2f}% des mutations trouvées de GATA2 ont été annotées")


# ########## EZH2 ##########

# EZH2_mut, MutationAnnotatedEZH2 = readFiles("EZH2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000106462_2025_03_20_11_46_31.csv")

# EZH2_var = AnnotateAllMutations(EZH2_mut, MutationAnnotatedEZH2)
# new_EZH2_var = WhereAreAnnotatedMutation(EZH2_var)
# print(f"{PercentOfAnnotatedMutations(new_EZH2_var)*100:.2f}% des mutations trouvées de EZH2 ont été annotées")


# ########## BCOR ##########

# BCOR_mut, MutationAnnotatedBCOR = readFiles("BCOR", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000183337_2025_03_20_11_48_45.csv")

# BCOR_var = AnnotateAllMutations(BCOR_mut, MutationAnnotatedBCOR)
# new_BCOR_var = WhereAreAnnotatedMutation(BCOR_var)
# print(f"{PercentOfAnnotatedMutations(new_BCOR_var)*100:.2f}% des mutations trouvées de BCOR ont été annotées")


# ########## JAK2 ##########

# JAK2_mut, MutationAnnotatedJAK2 = readFiles("JAK2", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000096968_2025_03_20_11_49_45.csv")

# JAK2_var = AnnotateAllMutations(JAK2_mut, MutationAnnotatedJAK2)
# new_JAK2_var = WhereAreAnnotatedMutation(JAK2_var)
# print(f"{PercentOfAnnotatedMutations(new_JAK2_var)*100:.2f}% des mutations trouvées de JAK2 ont été annotées")


# ########## SMC1A ##########

# SMC1A_mut, MutationAnnotatedSMC1A = readFiles("SMC1A", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000072501_2025_03_20_11_50_50.csv")

# SMC1A_var = AnnotateAllMutations(SMC1A_mut, MutationAnnotatedSMC1A)
# new_SMC1A_var = WhereAreAnnotatedMutation(SMC1A_var)
# print(f"{PercentOfAnnotatedMutations(new_SMC1A_var)*100:.2f}% des mutations trouvées de SMC1A ont été annotées")


# ########## RAD21 ##########

# RAD21_mut, MutationAnnotatedRAD21 = readFiles("RAD21", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000164754_2025_03_20_11_51_44.csv")

# RAD21_var = AnnotateAllMutations(RAD21_mut, MutationAnnotatedRAD21)
# new_RAD21_var = WhereAreAnnotatedMutation(RAD21_var)
# print(f"{PercentOfAnnotatedMutations(new_RAD21_var)*100:.2f}% des mutations trouvées de RAD21 ont été annotées")


# ########## SF3B1 ##########

# SF3B1_mut, MutationAnnotatedSF3B1 = readFiles("SF3B1", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000115524_2025_03_20_11_52_27.csv")

# SF3B1_var = AnnotateAllMutations(SF3B1_mut, MutationAnnotatedSF3B1)
# new_SF3B1_var = WhereAreAnnotatedMutation(SF3B1_var)
# print(f"{PercentOfAnnotatedMutations(new_SF3B1_var)*100:.2f}% des mutations trouvées de SF3B1 ont été annotées")


# ########## CBL ##########

# CBL_mut, MutationAnnotatedCBL = readFiles("CBL", "AnnotatedMutations/gnomAD_v4.1.0_ENSG00000110395_2025_03_20_11_53_11.csv")

# CBL_var = AnnotateAllMutations(CBL_mut, MutationAnnotatedCBL)
# new_CBL_var = WhereAreAnnotatedMutation(CBL_var)
# print(f"{PercentOfAnnotatedMutations(new_CBL_var)*100:.2f}% des mutations trouvées de CBL ont été annotées")