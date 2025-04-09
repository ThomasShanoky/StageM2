import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Charger les données
data_beat_aml = pd.read_csv("Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv", comment="#")
ExpressionsFile = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_NormalizedExpression.tsv"
data_expressions = pd.read_csv(ExpressionsFile, sep="\t", comment="#")

index_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_index.tsv"
with open(index_file) as f:
    index = f.readlines()[0]
index_list = index.split("\t")
index_list = [ind[:6] for ind in index_list] #liste des individus indexés sur Transipédia

usable_cat = ["consensus_sex", "reportedRace", "reportedEthnicity", "CEBPA_Biallelic", "consensusAMLFusions", "isRelapse", "isDenovo", "isTransformed", "specificDxAtAcquisition_MDSMPN", "nonAML_MDSMPN_specificDxAtAcquisition", "cumulativeChemo", "priorMalignancyRadiationTx", "priorMDS", "priorMDSMoreThanTwoMths", "priorMDSMPN", "priorMDSMPNMoreThanTwoMths", "priorMPN", "priorMPNMoreThanTwoMths", "ELN2017", "diseaseStageAtSpecimenCollection", "specimenType", "totalDrug", "cumulativeTreatmentRegimenCount", "cumulativeTreatmentStageCount", "responseToInductionTx", "typeInductionTx", "mostRecentTreatmentType", "vitalStatus", "causeOfDeath", "FLT3-ITD", "NPM1", "RUNX1", "ASXL1", "TP53"]

# Genes = ['ABL1', 'ADA', 'ANKRD26', 'ASXL1', 'ASXL2', 'ATM', 'ATRX', 'BCL6', 'BCOR', 'BCORL1', 'BCR', 'BIRC3', 'BLM', 'BRAF', 'BRCA1', 'BRCA2', 'CALR', 'CARD11', 'CBL', 'CBLB', 'CDKN2A', 'CEBPA', 'CHEK2', 'CREBBP', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CUX1', 'DAXX', 'DDX41', 'DNM2', 'DNMT1', 'DNMT3A', 'EED', 'EP300', 'ETNK1', 'ETV6', 'EZH2', 'FAS', 'FBXW7', 'FLRT2', 'FLT3', 'GATA1', 'GATA2', 'GNAS', 'HNRNPK', 'HRAS', 'IDH1', 'IDH2', 'IKZF1', 'IKZF3', 'IL7R', 'JAK1', 'JAK2', 'JAK3', 'KAT6A', 'KDM6A', 'KDR', 'KIT', 'KLHDC8B', 'KLHL6', 'KMT2A', 'KMT2C', 'KRAS', 'LRRC4', 'MAP2K1', 'MLH1', 'MPL', 'MSH2', 'MSH6', 'MYC', 'MYD88', 'NBN', 'NF1', 'NOTCH1', 'NPAT', 'NPM1', 'NRAS', 'NSD1', 'NTRK3', 'P2RY2', 'PAX5', 'PDGFRA', 'PHF6', 'PML', 'PMS2', 'PRF1', 'PRPF40B', 'PRPF8', 'PTPN11', 'RAD21', 'RB1', 'RELN', 'RUNX1', 'SETBP1', 'SF1', 'SF3A1', 'SF3B1', 'SH2B3', 'SH2D1A', 'SMARCB1', 'SMC1A', 'SMC3', 'SRP72', 'SRSF2', 'STAG2', 'STAT3', 'STXBP2', 'SUZ12', 'TAL1', 'TERT', 'TET2', 'TNFRSF13B', 'TP53', 'TPMT', 'U2AF1', 'U2AF2', 'WAS', 'WRN', 'WT1', 'XPO1', 'ZRSR2']
Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] #prevenant de Leucegene : on prend les plus abondantes pour travailler sur un nombre limité de gènes (Sample count >= 10). Seul KMT2D n'est pas dans la liste des 140 gènes
Genes.sort()

data_mutation_files = [f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv" for gene in Genes]
data_mutation = [pd.read_csv(file, sep=",", comment="#")[["sampleID", "localisation", "ref", "alt"]] for file in data_mutation_files]
# print(data_mutation)



class GUI:

    def __init__(self, minIndThreshold, data_beat_aml=data_beat_aml, data_expressions=data_expressions, data_mutation=data_mutation, Genes=Genes, usable_cat=usable_cat, index_list=index_list):

        self.minIndThreshold = minIndThreshold

        # Importation des données dans l'objet
        self.data_beat_aml = data_beat_aml
        self.data_expressions = data_expressions
        self.data_mutation = data_mutation
        self.Genes = Genes
        self.usable_cat = usable_cat
        self.index_list = index_list

        self.window = tk.Tk()  # Création d'une fenêtre
        self.window.geometry("1100x675")
        self.window.config(bg='#87CEEB')  # Fond bleu ciel
        self.window.title("Analyse des Expressions Génétiques")

        # self.menu = tk.Menu(self.window)  # Création d'une barre de menu
        # self.filemenu = tk.Menu(self.menu, tearoff=0)  # Création d'un menu déroulant
        # self.menu.add_cascade(label="Options", menu=self.filemenu)  # Ajouter le menu déroulant à la barre de menu
        # self.menu.config(bg='#87CEEB', fg='#000000')
        # self.window.config(menu=self.menu)  # Ajouter la barre de menu à la fenêtre

        # Variables pour les listes déroulantes
        self.gene_var = tk.StringVar(value=self.Genes[0])
        self.feature_var = tk.StringVar(value=self.usable_cat[0])

        # Listes déroulantes pour le gène et la feature et la mutation (si sélectionnée)
        self.gene_label = tk.Label(self.window, text="Sélectionnez un gène (1, 2 & 3) :")
        self.gene_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.gene_label.place(relx=0.03, rely=0.05)
        self.gene_dropdown = ttk.Combobox(self.window, textvariable=self.gene_var, values=self.Genes)
        self.gene_dropdown.place(relx=0.03, rely=0.10)
        self.gene_dropdown.bind("<<ComboboxSelected>>", self.update_mutations)

        self.dico_IndAndMut, self.dico_mut = self.getTypesOfMutationsAndInd(self.Genes[0])
        self.mut_var = tk.StringVar(value=self.format_mutations(self.dico_mut)[0])

        self.mut_label = tk.Label(self.window, text="Sélectionnez une mutation (3) :")
        self.mut_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.mut_label.place(relx=0.03, rely=0.18)
        self.mut_dropdown = ttk.Combobox(self.window, textvariable=self.mut_var, values=self.format_mutations(self.dico_mut), width=45)
        self.mut_dropdown.place(relx=0.03, rely=0.23)

        self.feature_label = tk.Label(self.window, text="Sélectionnez une feature (1 & 3) :")
        self.feature_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feature_label.place(relx=0.03, rely=0.31)
        self.feature_dropdown = ttk.Combobox(self.window, textvariable=self.feature_var, values=self.usable_cat)
        self.feature_dropdown.place(relx=0.03, rely=0.36)

        # Bouton interrupteur (switch on/off)
        self.switch_state = False
        self.switch_button_expr = tk.Button(self.window, text="Expressions de BEAT AML", command=self.toggle_switch, font=("DejaVu Serif", 13), bg="#FF6347", fg="#FFFFFF", width=20)
        self.switch_button_expr.place(relx=0.03, rely=0.40)

        #Bouton pour activer ou non l'analyse d'expression des gènes selon la feature choisie
        self.switch_feat = tk.BooleanVar(value=False)
        self.switch_feat_button = tk.Checkbutton(self.window, text="1. Expression selon les features", variable=self.switch_feat)
        self.switch_feat_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.switch_feat_button.place(relx=0.03, rely=0.47)

        # Bouton pour activer ou non l'analyse d'expression des gènes selon les mutations/features
        self.switch_var  = tk.BooleanVar(value=False)
        self.switch_button = tk.Checkbutton(self.window, text="2. Expression selon les mutations", variable=self.switch_var)
        self.switch_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.switch_button.place(relx=0.03, rely=0.52)

        # Bouton pour activer ou non l'analyse d'expression du gène ayant la mutation sélectionnée selon la feature choisie
        self.switch_mut = tk.BooleanVar(value=False)
        self.switch_mut_button = tk.Checkbutton(self.window, text="3. Expression de cette mutation selon\nla feature", variable=self.switch_mut)
        self.switch_mut_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.switch_mut_button.place(relx=0.03, rely=0.57)

        # Bouton pour générer le graphe
        self.generate_button = tk.Button(self.window, text="Générer le graphe", command=self.generate_plot)
        self.generate_button.config(width=15, font=("DejaVu Serif", 20, "bold"), highlightbackground="#370028", bg="#87CEEB", fg="#000000")
        self.generate_button.place(relx=0.03, rely=0.65)

        # Résultats des tests statistiques
        self.stats_label = tk.Label(self.window, text="")
        self.stats_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.stats_label.place(relx=0.03, rely=0.75)

        # Choix de sauvegarder les résultats
        self.SaveBool = tk.BooleanVar()
        self.save_button = tk.Checkbutton(self.window, text="Sauvegarder les résultats", variable=self.SaveBool)
        self.save_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.save_button.place(relx=0.03, rely=0.84)

        # Canvas pour afficher le graphe
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(relx=0.42, rely=0.05)

        self.window.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.window.mainloop()


    def toggle_switch(self):
        """Change l'état du switch et met à jour le texte et la couleur."""
        self.switch_state = not self.switch_state
        if self.switch_state:
            self.switch_button_expr.config(text="Expression kmers", bg="#32CD32")
        else:
            self.switch_button_expr.config(text="Expressions de BEAT AML", bg="#FF6347")
        return 
    

    def on_closing(self):
        self.window.quit()
        self.window.destroy()
    

    def generate_plot(self):
        if (self.switch_feat.get() and self.switch_var.get()) or (self.switch_feat.get() and self.switch_mut.get()) or (self.switch_var.get() and self.switch_mut.get()) or (self.switch_feat.get() and self.switch_var.get() and self.switch_mut.get()):
            messagebox.showwarning("Choix invalide", "Vous ne pouvez pas choisir plus d'une option")
            return

        if self.switch_feat.get() :
            self.generate_plot_feature()
        elif self.switch_var.get():
            self.generate_plot_mutations()
        elif self.switch_mut.get():
            self.generate_plot_mutations_features()
        return
    

    ##### Générer les plots avec les features #####

    def getSamplesAndFeatures(self, feature, featureValues):
        samples = self.data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample", feature]]
        samples = samples[samples[feature].isin(featureValues)]  # ne prendre que les lignes qui ont une valeur présente dans featureValues
        ind_beataml = []
        ind_feature = []

        for _, sample in samples.iterrows():
            if pd.isna(sample["dbgap_dnaseq_sample"]):
                sample_id = sample["dbgap_rnaseq_sample"][:6]
            else:
                sample_id = sample["dbgap_dnaseq_sample"][:6]
            ind_beataml.append(sample_id)
            ind_feature.append(sample[feature])

        SamplesAndFeatures = pd.DataFrame({"Sample": ind_beataml, "Feature": ind_feature})

        return SamplesAndFeatures


    def generate_plot_feature(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrer un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrer une feature valide", f"La feature \"{feature}\" n'est pas valide")

        featureValues = self.data_beat_aml[feature].dropna().unique().tolist()
        featureValues = [val for val in featureValues if val not in [np.nan, "Unknown", "unknown", "UNKNOWN"]]

        ##### Avoir les individus dont on peut extraire les features #####
        SamplesAndFeatures = self.getSamplesAndFeatures(feature, featureValues)

        ##### Avoir les individus dont on peut extraire les expressions #####
        ind_expr = list(self.data_expressions.columns)[4:]
        ind_expr = [ind[:6] for ind in ind_expr]
        # print(len(ind_expr))

        ##### Prendre l'intersection des deux listes #####
        inter_ind_pd = pd.DataFrame({"Sample": "", "Feature": "", "ExpressionGene": 0}, index=[0])

        for _, row in SamplesAndFeatures.iterrows():
            Id = row["Sample"]
            featureVal = row["Feature"]
            if Id in ind_expr:
                Expression = self.data_expressions.loc[self.data_expressions['display_label'] == gene, Id + "R"].values[0]
                inter_ind_pd = inter_ind_pd._append({"Sample": Id, "Feature": featureVal, "ExpressionGene": Expression}, ignore_index=True)

        inter_ind_pd = inter_ind_pd.drop(0).reset_index(drop=True)
        # print(inter_ind_pd)

        ##### Effectuer les tests statistiques #####
        if len(inter_ind_pd["Feature"].unique()) > 2:
            #ANOVA
            groups = [inter_ind_pd[inter_ind_pd["Feature"] == featureVal]["ExpressionGene"] for featureVal in inter_ind_pd["Feature"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.3f}"
        elif len(inter_ind_pd["Feature"].unique()) == 2:
            #Mann-Whitney U test
            group1 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[0]]["ExpressionGene"]
            group2 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[1]]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney U : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.3f}"
        else:
            messagebox.showwarning("Pas assez de features", "Il n'y a pas assez de valeurs de features pour effectuer un test statistique")

        if p_val < 0.05:
            significance = "La différence est significative"
        else:
            significance = "La différence n'est pas significative"

        self.stats_label.config(text=stats_res_test+"\n"+significance)

        ##### Créer le graphe #####
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, ax=self.ax)
        sns.stripplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction de la feature {feature}")
        self.ax.set_ylabel("Expression du gène")
        self.ax.set_xticklabels(self.ax.get_xticklabels(), rotation=45, horizontalalignment='right')
        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get():
            premiereLigne = f"#Boxplots d'expressions du gène {gene} en fonction de la feature {feature}\n"
            fileName = f"FeatureOnly_{gene}_{feature}"
            self.CreateFileRes(inter_ind_pd, stats_res_test, premiereLigne, fileName)

        return


    ##### Générer les plots avec les mutations #####

    def getTypesOfMutationsAndInd(self, gene):
        """Retourne un dictionnaire avec les types de mutations et les individus associés"""
        
        geneMutations = self.data_mutation[self.Genes.index(gene)]
        dico_mut = {}
        dico_IndAndMut = {}
        for _, row in geneMutations.iterrows():
            Id = row["sampleID"]
            local = row["localisation"]
            alt = row["alt"]
            ref = row["ref"]

            if [local, alt, ref] not in dico_mut.values():
                type_mut = "mut" + str(len(dico_mut)+1)
                dico_mut[type_mut] = [local, alt, ref]
                dico_IndAndMut[type_mut] = []

            dico_IndAndMut[type_mut].append(Id[:6])

        return dico_IndAndMut, dico_mut
    

    def checkIfPatientsAreInExpressions(self, dico_IndAndMut):
        """Retourne la liste d'individus à supprimer (qui ne sont pas dans les données d'expressions)"""
        Ind = dico_IndAndMut.values()
        Ind = [item for sublist in Ind for item in sublist]

        IndToRemove = []
        for ind in Ind:
            if ind+"R" not in self.data_expressions.columns:
                IndToRemove.append(ind)

        return IndToRemove
    

    def filterDicoIndAndMut(self, dico_IndAndMut, minInd):
        """Supprime les mutations qui n'ont pas assez d'individus associés"""

        mutToDel = []
        for mut in dico_IndAndMut:
            if len(dico_IndAndMut[mut]) < minInd:
                mutToDel.append(mut)

        dico_IndAndMut = {key: value for key, value in dico_IndAndMut.items() if key not in mutToDel}

        return dico_IndAndMut


    def getTable(self, gene, dico_IndAndMut):
        """Retourne un tableau avec 3 colonnes : l'ID sample, le type de mutation et l'expression du gène"""
        Tableau = pd.DataFrame(columns=["ID_sample", "TypeMutation", "ExpressionGene"])

        for mut in dico_IndAndMut.keys():
            for Id in dico_IndAndMut[mut]:
                Expression = self.data_expressions.loc[self.data_expressions['display_label'] == gene, Id + "R"].values[0]
                Tableau = Tableau._append({"ID_sample": Id, "TypeMutation": mut, "ExpressionGene": Expression}, ignore_index=True)

        return Tableau
    

    def getNonMutatedInd(self, AllIndMutated):
        NonMutatedInd = []
        AllInd = list(self.data_expressions.columns[4:]) #Tous les individus présents dans le fichier d'expressions géniques
        AllInd = [ind[:6] for ind in AllInd]
        
        for ind in AllInd:
            if ind not in AllIndMutated:
                NonMutatedInd.append(ind)
        return NonMutatedInd
    

    def checkIfIndExpressionAreInIndex(self, NonMutatedInd):
        IndIndex = self.index_list
        IndToDel = []

        for ind in NonMutatedInd:
            if ind not in IndIndex:
                IndToDel.append(ind)

        new_NonMutatedInd = [ind for ind in NonMutatedInd if ind not in IndToDel]

        return new_NonMutatedInd
    

    def addNonMutatedIndExpression(self, Tableau, gene, nonMutatedInd):

        for ind in nonMutatedInd:
            Expression = self.data_expressions.loc[self.data_expressions['display_label'] == gene, ind+"R"].values[0]
            Tableau = Tableau._append({"ID_sample": ind, "TypeMutation": "NonMut", "ExpressionGene": Expression}, ignore_index=True)

        return Tableau


    def generate_plot_mutations(self):
    
        gene = self.gene_var.get()
        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrer un gène valide", f"Le gène \"{gene}\" n'est pas valide")

        dico_IndAndMut, dico_mut = self.getTypesOfMutationsAndInd(gene)

        AllIndMutated = list(dico_IndAndMut.values())
        AllIndMutated = [item for sublist in AllIndMutated for item in sublist] #à utiliser pour avoir les individus NON mutés

        IndToRemove = self.checkIfPatientsAreInExpressions(dico_IndAndMut)

        if len(IndToRemove) > 0:
            messagebox.showwarning("Individus non trouvés", f"Les individus suivants n'ont pas d'expressions de gène associées :\n{IndToRemove}")
            for ind in IndToRemove:
                for mut in dico_IndAndMut.keys():
                    if ind in dico_IndAndMut[mut]:
                        dico_IndAndMut[mut].remove(ind)

        dico_IndAndMut = self.filterDicoIndAndMut(dico_IndAndMut, self.minIndThreshold)

        Tableau = self.getTable(gene, dico_IndAndMut)

        NonMutatedInd = self.getNonMutatedInd(AllIndMutated)
        NonMutatedInd = self.checkIfIndExpressionAreInIndex(NonMutatedInd)

        Tableau = self.addNonMutatedIndExpression(Tableau, gene, NonMutatedInd)
        # print(Tableau)

        ##### Effectuer les tests statistiques #####
        
        if len(dico_IndAndMut) > 1:
            # ANOVA
            groups = [Tableau[Tableau["TypeMutation"] == mut]["ExpressionGene"] for mut in Tableau["TypeMutation"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.3f}"
        elif len(dico_IndAndMut) == 1: #une seule mutation : on ne compare qu'avec les non mutés
            # Mann-Whitney U test
            #mutés
            group1 = Tableau[Tableau["TypeMutation"] == list(dico_IndAndMut.keys())[0]]["ExpressionGene"]
            #non mutés
            group2 = Tableau[Tableau["TypeMutation"] == "NonMut"]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney U : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.3f}"
        else:
            messagebox.showwarning("Pas assez de types de mutations", "Il n'y a pas assez de types de mutations pour effectuer un test statistique")
            # print(Tableau)

        if p_val < 0.05:
            significance = "La différence est significative"
        else:
            significance = "La différence n'est pas significative"

        self.stats_label.config(text=stats_res_test+"\n"+significance)

        ##### Créer le graphe #####
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="TypeMutation", y="ExpressionGene", data=Tableau, ax=self.ax)
        sns.stripplot(x="TypeMutation", y="ExpressionGene", data=Tableau, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction des types de mutations")
        self.ax.set_ylabel("Expression du gène")

        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get():
            premiereLigne = f"#Boxplots d'expressions du gène {gene} en fonction des types de mutations\n"
            fileName = f"MutOnly_{gene}"
            self.CreateFileRes(Tableau, stats_res_test, premiereLigne, fileName)
        
        return
    

    ##### Générer les plots avec les features selon les mutations du gène #####

    def update_mutations(self, event):
        gene = self.gene_var.get()
        self.dico_IndAndMut, self.dico_mut = self.getTypesOfMutationsAndInd(gene)
        self.mut_var.set(self.format_mutations(self.dico_mut)[0])
        self.mut_dropdown['values'] = self.format_mutations(self.dico_mut)


    def format_mutations(self, dico_mut):
        formatted_mutations = ["Toutes les mutations"]
        for mut, details in dico_mut.items():
            count = len(self.dico_IndAndMut[mut])
            if count >= self.minIndThreshold:
                formatted_mutations.append(f"{mut}: {details[2]}>{details[1]} ({count} patients)")
        return formatted_mutations
    

    def CreateFileRes(self, Tableau, stats_res_test, premiereLigne, fileName):

        output_brut_file = f"Documents/ScriptsPrincipaux/4_Expressions/DossierRes/{fileName}.txt"
        if not(os.path.exists("Documents/ScriptsPrincipaux/4_Expressions/DossierRes/")):
            os.mkdir("Documents/ScriptsPrincipaux/4_Expressions/DossierRes/")

        with open(output_brut_file, "w") as f:
            f.write(premiereLigne)
            f.write(f"#Résultat du test statistique {stats_res_test}\n")
            f.write("#Tableau des résultats :\n")
            f.write(Tableau.to_string(index=False))

        output_graph_file = f"Documents/ScriptsPrincipaux/4_Expressions/DossierRes/{fileName}.png"
        self.fig.savefig(output_graph_file)

        print(f"Les résultats bruts et le graphique ont été sauvegardés respectivement dans le fichier {output_brut_file} et dans le fichier {output_graph_file}")

        return
    

    def generate_plot_mutations_features(self):
        Mutation = self.mut_var.get()
        if Mutation == "Toutes les mutations":
            messagebox.showwarning("Attention", "Veuillez sélectionner une seule mutation")
        Mutation = Mutation.split(":")[0]

        PatientsRelatedToMut = self.dico_IndAndMut[Mutation]
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrer un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrer une feature valide", f"La feature \"{feature}\" n'est pas valide")

        featureValues = self.data_beat_aml[feature].dropna().unique().tolist()
        featureValues = [val for val in featureValues if val not in [np.nan, "Unknown", "unknown", "UNKNOWN"]]

        ##### Avoir les individus dont on peut extraire les features #####
        SamplesAndFeatures = self.getSamplesAndFeatures(feature, featureValues)

        ##### Prendre l'intersection des deux listes #####
        inter_ind_pd = pd.DataFrame({"Sample": "", "Feature": "", "ExpressionGene": 0}, index=[0])

        for _, row in SamplesAndFeatures.iterrows():
            Id = row["Sample"]
            featureVal = row["Feature"]
            if Id in PatientsRelatedToMut and Id + "R" in self.data_expressions.columns:
                Expression = self.data_expressions.loc[self.data_expressions['display_label'] == gene, Id + "R"].values[0]
                inter_ind_pd = inter_ind_pd._append({"Sample": Id, "Feature": featureVal, "ExpressionGene": Expression}, ignore_index=True)

        inter_ind_pd = inter_ind_pd.drop(0).reset_index(drop=True)
        # print(inter_ind_pd)

        ##### Effectuer les tests statistiques #####
        if len(inter_ind_pd["Feature"].unique()) > 2:
            # ANOVA
            groups = [inter_ind_pd[inter_ind_pd["Feature"] == featureVal]["ExpressionGene"] for featureVal in inter_ind_pd["Feature"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.3f}"
        elif len(inter_ind_pd["Feature"].unique()) == 2:
            # Mann-Whitney U test
            group1 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[0]]["ExpressionGene"]
            group2 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[1]]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney U : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.3f}"
        else:
            messagebox.showwarning("Pas assez de features", "Il n'y a pas assez de valeurs de features pour effectuer un test statistique")

        if p_val < 0.05:
            significance = "La différence est significative"
        else:
            significance = "La différence n'est pas significative"

        self.stats_label.config(text=stats_res_test+"\n"+significance)

        ##### Créer le graphe #####
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, ax=self.ax)
        sns.stripplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction de la feature {feature} pour la mutation {Mutation}")
        self.ax.set_ylabel("Expression du gène")
        self.ax.set_xticklabels(self.ax.get_xticklabels(), rotation=45, horizontalalignment='right')
        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get():
            premiereLigne = f"#Boxplots d'expressions de la mutation {Mutation} ({self.dico_mut[Mutation][1]}>{self.dico_mut[Mutation][2]} aux positions {self.dico_mut[Mutation][0]}) du gène {gene}, en fonction de la feature {feature}\n"
            fileName = f"MutAndFeat_{gene}_{Mutation}_{feature}"
            self.CreateFileRes(inter_ind_pd, stats_res_test, premiereLigne, fileName)

        return



GUI(4)