import numpy as np
import pandas as pd
import os
from tabulate import tabulate
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import RandomOverSampler

np.random.seed(42)


beataml_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv"
data_beat_aml = pd.read_csv(beataml_file)
data_beat_aml["sample"] = data_beat_aml.apply(
    lambda row: row["dbgap_dnaseq_sample"][:6] if not pd.isna(row["dbgap_dnaseq_sample"]) else row["dbgap_rnaseq_sample"][:6],
    axis=1
)
data_beat_aml.set_index("sample", inplace=True)
# data_beat_aml = data_beat_aml[data_beat_aml["diseaseStageAtSpecimenCollection"] == "Initial Diagnosis"]
# print(len(data_beat_aml))

beataml_expressions_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_NormalizedExpression.tsv"
beataml_expressions = pd.read_csv(beataml_expressions_file, sep="\t")
beataml_expressions.columns = list(beataml_expressions.columns[:3]) + [col[:6] for col in beataml_expressions.columns[3:]] 
beataml_expressions.set_index("display_label", inplace=True)
expr_ind = beataml_expressions.columns[3:].tolist()

index_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_index.tsv"


usable_cat = ["consensus_sex", "reportedRace", "reportedEthnicity", "CEBPA_Biallelic", "consensusAMLFusions", "isRelapse", "isDenovo", "isTransformed", "specificDxAtAcquisition_MDSMPN", "nonAML_MDSMPN_specificDxAtAcquisition", "cumulativeChemo", "priorMalignancyRadiationTx", "priorMDS", "priorMDSMoreThanTwoMths", "priorMDSMPN", "priorMDSMPNMoreThanTwoMths", "priorMPN", "priorMPNMoreThanTwoMths", "ELN2017", "diseaseStageAtSpecimenCollection", "specimenType", "totalDrug", "cumulativeTreatmentRegimenCount", "cumulativeTreatmentStageCount", "responseToInductionTx", "typeInductionTx", "mostRecentTreatmentType", "vitalStatus", "causeOfDeath", "FLT3-ITD", "NPM1", "RUNX1", "ASXL1", "TP53"] #liste de feature utilisable pour la classification

# GenesBA = ['ABL1', 'ADA', 'ANKRD26', 'ASXL1', 'ASXL2', 'ATM', 'ATRX', 'BCL6', 'BCOR', 'BCORL1', 'BCR', 'BIRC3', 'BLM', 'BRAF', 'BRCA1', 'BRCA2', 'CALR', 'CARD11', 'CBL', 'CBLB', 'CDKN2A', 'CEBPA', 'CHEK2', 'CREBBP', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CUX1', 'DAXX', 'DDX41', 'DNM2', 'DNMT1', 'DNMT3A', 'EED', 'EP300', 'ETNK1', 'ETV6', 'EZH2', 'FAS', 'FBXW7', 'FLRT2', 'FLT3', 'GATA1', 'GATA2', 'GNAS', 'HNRNPK', 'HRAS', 'IDH1', 'IDH2', 'IKZF1', 'IKZF3', 'IL7R', 'JAK1', 'JAK2', 'JAK3', 'KAT6A', 'KDM6A', 'KDR', 'KIT', 'KLHDC8B', 'KLHL6', 'KMT2A', 'KMT2C', 'KRAS', 'LRRC4', 'LUC7L2', 'MAP2K1', 'MLH1', 'MPL', 'MSH2', 'MSH6', 'MYC', 'MYD88', 'NBN', 'NF1', 'NOTCH1', 'NPAT', 'NPM1', 'NRAS', 'NSD1', 'NTRK3', 'P2RY2', 'PAX5', 'PDGFRA', 'PHF6', 'PML', 'PMS2', 'PRF1', 'PRPF40B', 'PRPF8', 'PTPN11', 'RAD21', 'RB1', 'RELN', 'RUNX1', 'SETBP1', 'SF1', 'SF3A1', 'SF3B1', 'SH2B3', 'SH2D1A', 'SMARCB1', 'SMC1A', 'SMC3', 'SRP72', 'SRSF2', 'STAG2', 'STAT3', 'STXBP2', 'SUZ12', 'TAL1', 'TERT', 'TET2', 'TNFRSF13B', 'TP53', 'TPMT', 'TUBA3C', 'U2AF1', 'U2AF2', 'WAS', 'WRN', 'WT1', 'XPO1', 'ZRSR2'] 
Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] #prevenant de Leucegene : on prend les plus abondantes pour travailler sur un nombre limité de gènes (Sample count >= 10). Seul KMT2D n'est pas dans la liste des 140 gènes



###########################
# Préparation des données #
###########################

##### Récupérer la liste d'individus (et les abunds) dont on a détecté une mutation pour chaque gène (Script de Sandra) #####

def getIndMut(Genes:list[str]) -> list[list[str]]:
    IndMut = [] #liste de liste d'individus qui ont le gène i muté (i étant entre 0 et n-1, n étant le nombre de gènes)
    for gene in Genes:
        file = f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv"
        IndMut.append(pd.read_csv(file, sep=",")["sampleID"].tolist())
        IndMut[-1] = [name[:6] for name in IndMut[-1]] #on ne garde que les 6 premiers caractères de chaque Id Sample d'individu
    return IndMut


def getAbund(IndMut:list[list[str]], Genes:list[str]) -> list[list[float]]:
    Abundances = []
    for gene_nb, lst in enumerate(IndMut):
        file = f"Documents/ScriptsPrincipaux/newMUTdata/{Genes[gene_nb]}_alt_perso.csv"
        Abundances.append(pd.read_csv(file, sep=",")["mean_count_kmer_alt"].tolist())
        # Abundances[-1] = [float(Abund[:-1]) for Abund in Abundances[-1]]
    return Abundances
    

##### Récupérer la liste d'individus de BEATAML (dont on va extraire les valeurs de la feature choisie) #####

def getIndFeatAndBeatAml(feature):
    samples = data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample", feature]]
    ind_beataml = []
    ind_feature = []

    for row, sample in samples.iterrows():
        if pd.isna(sample["dbgap_dnaseq_sample"]):
            sample_id = sample["dbgap_rnaseq_sample"][:6]
        else:
            sample_id = sample["dbgap_dnaseq_sample"][:6]
        ind_beataml.append(sample_id)
        ind_feature.append(sample[feature])
    return ind_beataml, ind_feature


##### Récupérer la liste des individus indexés (dont la mutation du gène i PEUT être détectée par l'outil de Sandra) #####

with open(index_file, "r") as f:
    index_list = f.readlines()[0].split("\t")
index_list = [ind[:6] for ind in index_list]


##### Construire le tableau #####


def creteTable(index, Genes, IndMut, SamplesAndFeatures, featureValues, DataBool, Abundances, Expressions):
    Tableau = pd.DataFrame(index=index, columns=Genes, dtype=float)
    Tableau = Tableau.fillna(0)
    Tableau["Feature"] = ""

    for lst in IndMut:
        for ind in lst:
            if ind in index:
                if DataBool == 0:
                    Tableau.at[ind, Genes[IndMut.index(lst)]] = Abundances[IndMut.index(lst)][lst.index(ind)]
                elif DataBool == 1:
                    Tableau.at[ind, Genes[IndMut.index(lst)]] = 1
                elif DataBool == 2:
                    Tableau.at[ind, Genes[IndMut.index(lst)]] = Expressions.at[Genes[IndMut.index(lst)], ind]
                Tableau.at[ind, "Feature"] = SamplesAndFeatures[SamplesAndFeatures["Sample"] == ind]["Feature"].values[0]

    Tableau = Tableau[Tableau["Feature"].isin(featureValues)] #filtrer les individus qui n'ont pas les valeurs enregistrées pour la feature choisie

    return Tableau


# pd.set_option('display.max_rows', None)  # Afficher toutes les lignes
# pd.set_option('display.max_columns', None)  # Afficher toutes les colonnes
# pd.set_option('display.width', None)  # Ajuster la largeur de l'affichage
# pd.set_option('display.max_colwidth', None)  # Ajuster la largeur maximale des colonnes



##### Filtrer les individus qui n'ont pas les valeurs enregistrées pour la feature choisie #####

import sklearn

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV



##################################
# Modèle de classification : SVM #
##################################

from sklearn.svm import SVC
from sklearn.inspection import permutation_importance


def TrainSVM(X_train, y_train, X_test, y_test, gridSearchBool):
    if gridSearchBool:
        param_grid = {
            'C': [0.1, 1, 10, 100],
            'kernel': ['linear', 'rbf', 'poly'],
            'probability': [True, False]
        }
        grid_search = GridSearchCV(SVC(random_state=42), param_grid, cv=5)
        grid_search.fit(X_train, y_train)

    # print(f"Meilleurs paramètres : {grid_search.best_params_}")

        clf = grid_search.best_estimator_
    else:
        clf = SVC(random_state=42)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    #afficher les features les plus importantes

    result = permutation_importance(clf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=-1)
    importances = np.abs(result.importances_mean)
    indices = np.argsort(importances)[::-1]

    # print("Ordre des feature les plus importantes (du plus au moins) :")
    features = Genes + ["Âge"]
    res = ""
    for f in range(7):
        res += f"{f + 1}. feature {features[indices[f]]} ({importances[indices[f]]:.4f})\n"

    return acc, res


#################################################
# Modèle de classification : Decision Tree (DT) #
#################################################

from sklearn.tree import DecisionTreeClassifier

def TrainDT(X_train, y_train, X_test, y_test, gridSearchBool):
    if gridSearchBool:
        param_grid = {
            'criterion': ['gini', 'entropy'],
            'splitter': ['best', 'random'],
            'max_depth': [None, 10, 20, 30, 40, 50],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
            'max_features': [None, 'sqrt', 'log2']
        }
        grid_search = GridSearchCV(DecisionTreeClassifier(random_state=42), param_grid, cv=5)
        grid_search.fit(X_train, y_train)

        # print(f"Meilleurs paramètres : {grid_search.best_params_}")

        clf = grid_search.best_estimator_
    else:
        clf = DecisionTreeClassifier(random_state=42)
        
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    importances = clf.feature_importances_
    indices = np.argsort(importances)[::-1]

    # print("Ordre des feature les plus importantes (du plus au moins) :")
    features = Genes + ["Âge"]
    res = ""
    for f in range(7):
        res += f"{f + 1}. feature {features[indices[f]]} ({importances[indices[f]]:.4f})\n"
    
    return acc, res


############################################
# Modèle de classification : Random Forest #
############################################

from sklearn.ensemble import RandomForestClassifier

def TrainRF(X_train, y_train, X_test, y_test, GridSearchBool):
    if GridSearchBool:
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_features': ['sqrt', 'log2'],
            'max_depth': [None, 10, 20, 30],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
            'bootstrap': [True, False]
        }
        grid_search = GridSearchCV(RandomForestClassifier(random_state=42), param_grid, cv=5)
        grid_search.fit(X_train, y_train)

        # print(f"Meilleurs paramètres : {grid_search.best_params_}")

        clf = grid_search.best_estimator_
    else:
        clf = RandomForestClassifier(random_state=42)

    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    #afficher les features les plus importantes

    importances = clf.feature_importances_
    indices = np.argsort(importances)[::-1]

    # print("Ordre des feature les plus importantes (du plus au moins) :")
    features = Genes + ["Âge"]
    res = ""
    for f in range(7):
        res += f"{f + 1}. feature {features[indices[f]]} ({importances[indices[f]]:.4f})\n"

    return acc, res


######################################
# Modèle de classification : XGBoost #
######################################

from sklearn.ensemble import GradientBoostingClassifier

def TrainXGBoost(X_train, y_train, X_test, y_test, GridSearchBool):
    if GridSearchBool:
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [3, 5, 7, 10],
            'learning_rate': [0.01, 0.1, 0.2],
            'subsample': [0.6, 0.8, 1.0]
        }
        grid_search = GridSearchCV(GradientBoostingClassifier(random_state=42), param_grid, cv=5)
        grid_search.fit(X_train, y_train)

        # print(f"Meilleurs paramètres : {grid_search.best_params_}")

        clf = grid_search.best_estimator_
    else:
        clf = GradientBoostingClassifier(random_state=42)

    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    #afficher les features les plus importantes

    importances = clf.feature_importances_
    indices = np.argsort(importances)[::-1]

    # print("Ordre des feature les plus importantes (du plus au moins) :")
    features = Genes + ["Âge"]
    res = ""
    for f in range(7):
        res += f"{f + 1}. feature {features[indices[f]]} ({importances[indices[f]]:.4f})\n"
    
    return acc, res


#########################
# Régression logistique #
#########################

from sklearn.linear_model import LogisticRegression

def TrainLogReg(X_train, y_train, X_test, y_test, GridSearchBool):
    if GridSearchBool:
        param_grid = {
            'C': [0.1, 1, 10, 100],
            'penalty': ['l1', 'l2'],
            'solver': ['liblinear']
        }
        grid_search = GridSearchCV(LogisticRegression(random_state=42), param_grid, cv=5)
        grid_search.fit(X_train, y_train)

        # print(f"Meilleurs paramètres : {grid_search.best_params_}")

        clf = grid_search.best_estimator_
    else:
        clf = LogisticRegression(random_state=42)
    
    clf.fit(X_train, y_train)
    
    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    importances = clf.coef_[0]
    indices = np.argsort(importances)[::-1]

    #ordre des features les plus importantes
    features = Genes + ["Âge"]
    res = ""
    for f in range(7):
        res += f"{f + 1}. feature {features[indices[f]]} ({importances[indices[f]]:.4f})\n"

    return acc, res


#####################################################################





IndMut = getIndMut(Genes)
Abundances = getAbund(IndMut, Genes)

AbundBools = [0, 1, 2] #Si on veut utiliser les Abondances à la place de simples "1" pour les mutations
# UseGridSearchBool = False
Features = ["ELN2017", "vitalStatus", "isRelapse"] 
FeaturesValuesChosen = [["Favorable", "Adverse"], ["Alive", "Dead"], ["FALSE", "TRUE"]]

DonneesUtilisees = ["0/1", "Abondances", "Expressions"]


for i, feature in enumerate(Features):
    ResTableau = pd.DataFrame({"Données utilisées":DonneesUtilisees, "SVM":["", "", ""], "Decision Tree":["", "", ""], "Random Forest":["", "", ""], "XGBoost":["", "", ""], "Régression Logistique":["", "", ""]})
    ResTableau.set_index("Données utilisées", inplace=True)
    # print(ResTableau)
    featureValuesPossible = np.unique(data_beat_aml[feature].tolist())
    featureValues = FeaturesValuesChosen[i]
    for feat in featureValues:
        if feat not in featureValuesPossible:
            print(data_beat_aml["isRelapse"])
            raise ValueError(f"Les valeurs des features entrées ne sont pas dans les valeurs possibles de la feature {feature} : {featureValuesPossible}")
            
    for DataBool in AbundBools:
        # for UseGridSearchBool in [False, True]: 

        ind_beataml, ind_feature = getIndFeatAndBeatAml(feature)
        SamplesAndFeatures = pd.DataFrame({"Sample": ind_beataml, "Feature": ind_feature})

        indToRemove = [] #on vérifie que tous les individus de BEATAML sont bien dans la liste des individus indexés
        for ind in ind_beataml:
            if ind not in index_list or ind not in expr_ind:
                indToRemove.append(ind)

        ind_beataml = [ind for ind in ind_beataml if ind not in indToRemove]

        Tableau = creteTable(ind_beataml, Genes, IndMut, SamplesAndFeatures, featureValues, DataBool, Abundances, beataml_expressions)
        Tableau.insert(Tableau.shape[1]-1, "Âge", [0.0 for _ in range(Tableau.shape[0])])
        for index, row in Tableau.iterrows():
            if not pd.isna(data_beat_aml.at[index, "ageAtDiagnosis"]):
                Tableau.at[index, "Âge"] = data_beat_aml.at[index, "ageAtDiagnosis"]
            else:
                Tableau.at[index, "Âge"] = 0
        Tableau["Âge"] = [np.mean(Tableau["Âge"].dropna()) if age == 0 else float(age) for age in Tableau["Âge"]]
        print(Tableau)

        X = Tableau[Genes+["Âge"]]
        y = Tableau["Feature"]
        # print(feature)
        lst_temp = [val if val == FeaturesValuesChosen[i][0] else 0 for val in y.tolist()]
        # print(lst_temp.count(FeaturesValuesChosen[i][0]))
        # print(lst_temp.count(0))
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # ros = RandomOverSampler(random_state=42)
        # X_train, y_train = ros.fit_resample(X_train, y_train)

        scaler = StandardScaler()
        X_train, X_test = scaler.fit_transform(X_train), scaler.transform(X_test) #on normalise l'ensemble des données d'entrainement et on applique la même normalisation à l'ensemble de test

        accSVM, genesSVM = TrainSVM(X_train, y_train, X_test, y_test, False)
        resSVM = f"acc = {accSVM:.3f}\navec {genesSVM.split()[2]} et {genesSVM.split()[6]}"
        accDT, genesDT = TrainDT(X_train, y_train, X_test, y_test, False)
        resDT = f"acc = {accDT:.3f}\navec {genesDT.split()[2]} et {genesDT.split()[6]}"
        accRF, genesRF = TrainRF(X_train, y_train, X_test, y_test, False)
        resRF = f"acc = {accRF:.3f}\navec {genesRF.split()[2]} et {genesRF.split()[6]}"
        accXGB, genesXGB = TrainXGBoost(X_train, y_train, X_test, y_test, False)
        resXGB = f"acc = {accXGB:.3f}\navec {genesXGB.split()[2]} et {genesXGB.split()[6]}"
        accLogReg, genesLogReg = TrainLogReg(X_train, y_train, X_test, y_test, False)
        resLogReg = f"acc = {accLogReg:.3f}\navec {genesLogReg.split()[2]} et {genesLogReg.split()[6]}"

        output_file = f"Documents/ScriptsPrincipaux/3_Classification/DossierRes/output_{feature}.txt"

        if not(os.path.exists("Documents/ScriptsPrincipaux/3_Classification/DossierRes/")):
            os.makedirs("Documents/ScriptsPrincipaux/3_Classification/DossierRes/")

        ResTableau.at[DonneesUtilisees[AbundBools.index(DataBool)], "SVM"] = resSVM
        ResTableau.at[DonneesUtilisees[AbundBools.index(DataBool)], "Decision Tree"] = resDT
        ResTableau.at[DonneesUtilisees[AbundBools.index(DataBool)], "Random Forest"] = resRF
        ResTableau.at[DonneesUtilisees[AbundBools.index(DataBool)], "XGBoost"] = resXGB
        ResTableau.at[DonneesUtilisees[AbundBools.index(DataBool)], "Régression Logistique"] = resLogReg

        with open(output_file, 'w') as fw:
            fw.write(f"Feature {feature} with values : {featureValues}\n")
            fw.write(f"Using Abundances : {DataBool}\n")
            fw.write(tabulate(ResTableau, headers="keys", tablefmt="grid", stralign="center"))

        print(f"Feature {feature} avec {featureValues}; Données utilisées : {DonneesUtilisees[DataBool]} | terminé")

