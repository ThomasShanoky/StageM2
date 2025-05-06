#kmerator -s GAPDH -o resultsGAPDH => NIVEAU GENIQUE
#requête des kmers sur Transipédia avec les deux cohortes BEAT AML de patients (atteints d'AML)
#Count method = raw


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#Ici, on veut juste voir l'expression chez tous les patients de gène de ménage, tel que GAPDH. Cette expression doit être la même pour tous

def getMeanByPatient(file):
    ExpressionsKmers = pd.read_csv(file, sep="\t")

    #faire une moyenne de l'expression de GAPDH pour chaque patient
    ExpressionGeneMoy = ExpressionsKmers.drop(columns=["seq_name"])
    ExpressionGeneMoy = ExpressionGeneMoy.mean()

    return ExpressionGeneMoy, ExpressionGeneMoy.std()




file_HPRT1 = "Documents/DeuxiemePartie/GeneMenage/query_results_HPRT1.tsv"
file_RPLP0 = "Documents/DeuxiemePartie/GeneMenage/query_results_RPLP0.tsv"
file_GUSB = "Documents/DeuxiemePartie/GeneMenage/query_results_GUSB.tsv"
file_TBP = "Documents/DeuxiemePartie/GeneMenage/query_results_TBP.tsv"
file_PPIA = "Documents/DeuxiemePartie/GeneMenage/query_results_PPIA.tsv"
file_ABL1 = "Documents/DeuxiemePartie/GeneMenage/query_results_ABL1.tsv"
file_GAPDH = "Documents/DeuxiemePartie/GeneMenage/query_results_GAPDH.tsv"
file_ACTB = "Documents/DeuxiemePartie/GeneMenage/query_results_ACTB.tsv"


ExpressionMoyHPRT1, EcartTypeHPRT1 = getMeanByPatient(file_HPRT1)
ExpressionMoyRPLP0, EcartTypeRPLP0 = getMeanByPatient(file_RPLP0)
ExpressionMoyGUSB, EcartTypeGUSB = getMeanByPatient(file_GUSB)
ExpressionMoyTBP, EcartTypeTBP = getMeanByPatient(file_TBP)
ExpressionMoyPPIA, EcartTypePPIA = getMeanByPatient(file_PPIA)
ExpressionMoyABL1, EcartTypeABL1 = getMeanByPatient(file_ABL1)



all_means = pd.concat([ExpressionMoyHPRT1, ExpressionMoyGUSB, ExpressionMoyTBP, ExpressionMoyPPIA, ExpressionMoyABL1])
# all_means = pd.concat([ExpressionMoyNPM1, ExpressionMoyFLT3, ExpressionMoyHPRT1, ExpressionMoyGUSB, ExpressionMoyTBP, ExpressionMoyPPIA, ExpressionMoyABL1])
y_min = all_means.min()
y_max = all_means.max()


fig, axes = plt.subplots(1, 8, figsize=(20, 6))

# sns.boxplot(data=ExpressionMoyNPM1, ax=axes[0])
# sns.stripplot(data=ExpressionMoyNPM1, ax=axes[0], color=".25", alpha=0.3)
# axes[0].set_title("NPM1")
# # axes[0].set_yscale("log")
# axes[0].set_ylim(y_min, y_max)
# axes[0].annotate(rf"$\bar{{\sigma}}$: {EcartTypeNPM1:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

# sns.boxplot(data=ExpressionMoyFLT3, ax=axes[1])
# sns.stripplot(data=ExpressionMoyFLT3, ax=axes[1], color=".25", alpha=0.3)
# axes[1].set_title("FLT3")
# # axes[1].set_yscale("log")
# axes[1].set_ylim(y_min, y_max)
# axes[1].annotate(rf"$\bar{{\sigma}}$: {EcartTypeFLT3:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

sns.boxplot(data=ExpressionMoyHPRT1, ax=axes[3])
sns.stripplot(data=ExpressionMoyHPRT1, ax=axes[3], color=".25", alpha=0.3)
axes[3].set_title("HPRT1")
# axes[3].set_yscale("log")
axes[3].set_ylim(y_min, y_max)
axes[3].annotate(rf"$\bar{{\sigma}}$: {EcartTypeHPRT1:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

sns.boxplot(data=ExpressionMoyGUSB, ax=axes[4])
sns.stripplot(data=ExpressionMoyGUSB, ax=axes[4], color=".25", alpha=0.3)
axes[4].set_title("GUSB")
# axes[4].set_yscale("log")
axes[4].set_ylim(y_min, y_max)
axes[4].annotate(rf"$\bar{{\sigma}}$: {EcartTypeGUSB:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

sns.boxplot(data=ExpressionMoyTBP, ax=axes[5])
sns.stripplot(data=ExpressionMoyTBP, ax=axes[5], color=".25", alpha=0.3)
axes[5].set_title("TBP")
# axes[5].set_yscale("log")
axes[5].set_ylim(y_min, y_max)
axes[5].annotate(rf"$\bar{{\sigma}}$: {EcartTypeTBP:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

sns.boxplot(data=ExpressionMoyPPIA, ax=axes[6])
sns.stripplot(data=ExpressionMoyPPIA, ax=axes[6], color=".25", alpha=0.3)
axes[6].set_title("PPIA")
# axes[6].set_yscale("log")
axes[6].set_ylim(y_min, y_max)
axes[6].annotate(rf"$\bar{{\sigma}}$: {EcartTypePPIA:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

sns.boxplot(data=ExpressionMoyABL1, ax=axes[7])
sns.stripplot(data=ExpressionMoyABL1, ax=axes[7], color=".25", alpha=0.3)
axes[7].set_title("ABL1")
# axes[7].set_yscale("log")
axes[7].set_ylim(y_min, y_max)
axes[7].annotate(rf"$\bar{{\sigma}}$: {EcartTypeABL1:.2f}", xy=(0.5, -0.1), xycoords='axes fraction', ha='center', fontsize=12)

plt.tight_layout()
plt.show()


print(min(ExpressionMoyTBP))