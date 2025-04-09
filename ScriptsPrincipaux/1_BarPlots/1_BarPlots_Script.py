import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Charger les données
data_beat_aml = pd.read_csv("Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv", comment="#")
index_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_index.tsv"
with open(index_file) as f:
    index = f.readlines()[0]
index_list = index.split("\t")
index_list = [ind[:6] for ind in index_list]

usable_cat = ["consensus_sex", "reportedRace", "reportedEthnicity", "CEBPA_Biallelic", "consensusAMLFusions", "isRelapse", "isDenovo", "isTransformed", "specificDxAtAcquisition_MDSMPN", "nonAML_MDSMPN_specificDxAtAcquisition", "cumulativeChemo", "priorMalignancyRadiationTx", "priorMDS", "priorMDSMoreThanTwoMths", "priorMDSMPN", "priorMDSMPNMoreThanTwoMths", "priorMPN", "priorMPNMoreThanTwoMths", "ELN2017", "diseaseStageAtSpecimenCollection", "specimenType", "totalDrug", "cumulativeTreatmentRegimenCount", "cumulativeTreatmentStageCount", "responseToInductionTx", "typeInductionTx", "mostRecentTreatmentType", "vitalStatus", "causeOfDeath", "FLT3-ITD", "NPM1", "RUNX1", "ASXL1", "TP53"]

# Genes = ['ABL1', 'ADA', 'ANKRD26', 'ASXL1', 'ASXL2', 'ATM', 'ATRX', 'BCL6', 'BCOR', 'BCORL1', 'BCR', 'BIRC3', 'BLM', 'BRAF', 'BRCA1', 'BRCA2', 'CALR', 'CARD11', 'CBL', 'CBLB', 'CDKN2A', 'CEBPA', 'CHEK2', 'CREBBP', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CUX1', 'DAXX', 'DDX41', 'DNM2', 'DNMT1', 'DNMT3A', 'EED', 'EP300', 'ETNK1', 'ETV6', 'EZH2', 'FAS', 'FBXW7', 'FLRT2', 'FLT3', 'GATA1', 'GATA2', 'GNAS', 'HNRNPK', 'HRAS', 'IDH1', 'IDH2', 'IKZF1', 'IKZF3', 'IL7R', 'JAK1', 'JAK2', 'JAK3', 'KAT6A', 'KDM6A', 'KDR', 'KIT', 'KLHDC8B', 'KLHL6', 'KMT2A', 'KMT2C', 'KRAS', 'LRRC4', 'LUC7L2', 'MAP2K1', 'MLH1', 'MPL', 'MSH2', 'MSH6', 'MYC', 'MYD88', 'NBN', 'NF1', 'NOTCH1', 'NPAT', 'NPM1', 'NRAS', 'NSD1', 'NTRK3', 'P2RY2', 'PAX5', 'PDGFRA', 'PHF6', 'PML', 'PMS2', 'PRF1', 'PRPF40B', 'PRPF8', 'PTPN11', 'RAD21', 'RB1', 'RELN', 'RUNX1', 'SETBP1', 'SF1', 'SF3A1', 'SF3B1', 'SH2B3', 'SH2D1A', 'SMARCB1', 'SMC1A', 'SMC3', 'SRP72', 'SRSF2', 'STAG2', 'STAT3', 'STXBP2', 'SUZ12', 'TAL1', 'TERT', 'TET2', 'TNFRSF13B', 'TP53', 'TPMT', 'TUBA3C', 'U2AF1', 'U2AF2', 'WAS', 'WRN', 'WT1', 'XPO1', 'ZRSR2']
Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] #prevenant de Leucegene : on prend les plus abondantes pour travailler sur un nombre limité de gènes (Sample count >= 10). Seul KMT2D n'est pas dans la liste des 140 gènes
Genes.sort()


class GUI:

    def __init__(self, data_beat_aml=data_beat_aml, index_list=index_list, usable_cat=usable_cat, Genes=Genes):
        self.data_beat_aml = data_beat_aml
        self.index_list = index_list
        self.usable_cat = usable_cat
        self.Genes = Genes

        self.window = tk.Tk()  # Création d'une fenêtre
        self.window.geometry("1100x675")
        self.window.config(bg='#191830')
        self.window.title("Analyse de l'effet de mutation sur différentes features")

        # self.menu = tk.Menu(self.window)  # Création d'une barre de menu
        # self.filemenu = tk.Menu(self.menu, tearoff=0)  # Création d'un menu déroulant
        # self.menu.add_cascade(label="Options", menu=self.filemenu)  # Ajouter le menu déroulant à la barre de menu
        # self.menu.config(bg='#161533', fg='#b7d9ec')
        # self.window.config(menu=self.menu)  # Ajouter la barre de menu à la fenêtre

        # Variables pour les listes déroulantes
        self.gene_var = tk.StringVar(value=self.Genes[0])
        self.feature_var = tk.StringVar(value=self.usable_cat[0])

        # Listes déroulantes pour le gène et la feature
        self.gene_label = tk.Label(self.window, text="Sélectionnez un gène:")
        self.gene_label.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.gene_label.place(relx=0.03, rely=0.05)
        self.gene_dropdown = ttk.Combobox(self.window, textvariable=self.gene_var, values=self.Genes)
        self.gene_dropdown.place(relx=0.03, rely=0.10)

        self.feature_label = tk.Label(self.window, text="Sélectionnez une feature:")
        self.feature_label.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.feature_label.place(relx=0.03, rely=0.18)
        self.feature_dropdown = ttk.Combobox(self.window, textvariable=self.feature_var, values=self.usable_cat)
        self.feature_dropdown.place(relx=0.03, rely=0.23)

        self.DistributionVar = tk.BooleanVar()
        self.distribution_button = tk.Checkbutton(self.window, text="Montrer les distributions", variable=self.DistributionVar)
        self.distribution_button.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.distribution_button.place(relx=0.03, rely=0.30)

        self.AbondanceVar = tk.BooleanVar()
        self.abondance_button = tk.Checkbutton(self.window, text="Montrer l'abondance des mutations", variable=self.AbondanceVar)
        self.abondance_button.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.abondance_button.place(relx=0.03, rely=0.37)

        # Bouton pour générer le graphe
        self.generate_button = tk.Button(self.window, text="Générer le graphe", command=self.generate_plot)
        self.generate_button.config(width=15, font=("DejaVu Serif", 20, "bold"), highlightbackground="#370028", bg="#191830", fg="#b7d9ec")
        self.generate_button.place(relx=0.03, rely=0.47)

        # Canvas pour afficher le graphe
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(relx=0.42, rely=0.05)

        # Choix de sauvegarder les résultats
        self.SaveBool = tk.BooleanVar()
        self.save_button = tk.Checkbutton(self.window, text="Sauvegarder les résultats", variable=self.SaveBool)
        self.save_button.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.save_button.place(relx=0.03, rely=0.60)

        # Labels pour afficher la p-value et la significativité
        self.p_value_label = tk.Label(self.window, text="")
        self.p_value_label.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.p_value_label.place(relx=0.03, rely=0.70)

        self.significance_label = tk.Label(self.window, text="")
        self.significance_label.config(font=("DejaVu Serif", 13), bg="#191830", fg="#b7d9ec")
        self.significance_label.place(relx=0.03, rely=0.75)

        self.window.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.window.mainloop()


    def on_closing(self):
        self.window.quit()
        self.window.destroy()

    
    def generate_plot(self):
        if self.AbondanceVar.get() and self.DistributionVar.get():
            messagebox.showwarning("Attention", "Vous ne pouvez pas choisir les deux options en même temps")
            return
        
        if not(self.AbondanceVar.get()) and not(self.DistributionVar.get()):
            messagebox.showwarning("Attention", "Veuillez choisir au moins une option")
            return 
        
        if self.DistributionVar.get():
            self.generate_plot_without_abundance()
            return
        if self.AbondanceVar.get():
            self.generate_plot_with_abundance()
            return


    def get_number_for_bar(self, cat_name, ind_beataml, ind_list, feature):
        number_for_bar = [0 for _ in range(len(cat_name))]
        
        for ind in ind_list:
            index_ind = ind_beataml.index(ind)
            OneCat = self.data_beat_aml.iloc[index_ind][feature]

            for i, cat in enumerate(cat_name):
                if OneCat == cat:
                    number_for_bar[i] += 1
        return number_for_bar
    

    def get_number_for_bar_spe(self, cat_name:list[str], ind_beataml, ind_list:list[str], feature:str)->list[int]:
        number_for_bar = [0, 0]
        
        for ind in ind_list:
            index_ind = ind_beataml.index(ind)
            OneCat = self.data_beat_aml.iloc[index_ind][feature]

            if pd.isna(OneCat) or OneCat == 'nan':
                number_for_bar[1] += 1
            else:
                number_for_bar[0] += 1

        return number_for_bar


    def Chi2Test(self, number_for_plot, number_for_plot_nonMut):
        table = np.array([number_for_plot, number_for_plot_nonMut])
        print(table)
        table = 100 * table / np.sum(table)

        _, p, _, expected = stats.chi2_contingency(table)
        residus = (table - expected) / np.sqrt(expected)

        return p, residus


    def rearrangeZeros(self, number_for_plot, number_for_plot_nonMut):
        """S'il y a une colonne remplie de 0, on la supprime, car cela pose problème pour le test du chi2 (en plus de ne pas donner d'informations supplémentaires)"""

        L_ind = [] #liste des indices dont la colonne est remplie de 0

        for i in range(len(number_for_plot)):
            if number_for_plot[i] == 0 and number_for_plot_nonMut[i] == 0:
                L_ind.append(i)

        return L_ind
    

    def CreateFileResFeat(self, ind_beataml, gene, feature, gene_cat, number_for_plot, number_for_plot_nonMut, p, ind_geneMut, ind_geneNonMut):
        output_file = f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/Resultats_Barplot.txt"

        if not os.path.exists("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes"):
            os.mkdir("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes")

        with open(output_file, 'w') as f:
            f.write(f"Résultats de barplot pour {gene} et {feature}\n")
            f.write(f"p-value = {p:.5f} selon le test du χ²\n")
            f.write(f"Résultats bruts : \n")
            Res = pd.DataFrame({
                f"{feature}": gene_cat,
                f"{gene} muté": number_for_plot,
                f"{gene} non muté": number_for_plot_nonMut
            })
            f.write(Res.to_string())
            f.write("\n")

            f.write("Liste des ID Sample par catégorie :\n")
            f.write(f"{gene} muté :\n")
            for i, cat in enumerate(gene_cat):
                if feature in ["TP53", "RUNX1", "ASXL1"]:
                    if cat == "Positive":
                        ids = [ind for ind in ind_geneMut if not pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                    else:
                        ids = [ind for ind in ind_geneMut if pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                else:
                    ids = [ind for ind in ind_geneMut if self.data_beat_aml.loc[ind_beataml.index(ind), feature] == cat]
                f.write(f"{cat} : {', '.join(ids)}\n")

            f.write(f"{gene} non muté :\n")
            for i, cat in enumerate(gene_cat):
                if feature in ["TP53", "RUNX1", "ASXL1"]:
                    if cat == "Positive":
                        ids = [ind for ind in ind_geneNonMut if not pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                    else:
                        ids = [ind for ind in ind_geneNonMut if pd.isna(self.data_beat_aml.loc[ind_beataml.index(ind), feature])]
                else:
                    ids = [ind for ind in ind_geneNonMut if self.data_beat_aml.loc[ind_beataml.index(ind), feature] == cat]
                f.write(f"{cat} : {', '.join(ids)}\n")

        self.fig.savefig(f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/Barplot_{gene}_{feature}.png")

        print(f"Les résultats bruts et la figure ont été sauvegardés dans les fichiers {output_file} et Barplot_{gene}_{feature}.png")

        return
    

    def plot_graph_without_abundance(self, cat_name, number_for_plot, number_for_plot_nonMut):
        bar_width = 0.35

        r1 = np.arange(len(cat_name))
        r2 = [x + bar_width for x in r1]

        gene = self.gene_var.get()
        feature = self.feature_var.get()

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.bar(r1, number_for_plot, color="blue", width=bar_width, label=f"{gene} mutated")
        self.ax.bar(r2, number_for_plot_nonMut, color="green", width=bar_width, label=f"{gene} non mutated")

        self.ax.set_xticks([r + bar_width / 2 for r in range(len(cat_name))])
        self.ax.set_xticklabels(cat_name, rotation=90)
        self.ax.set_ylabel("Nombre d'échantillons")
        self.ax.legend()
        self.ax.set_title(f"Distribution de {feature} selon la mutation (ou non) de {gene}")

        self.fig.tight_layout()
        self.canvas.draw()
        return


    def generate_plot_without_abundance(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrez un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrez une feature valide", f"La feature \"{feature}\" n'est pas valide")

        data_gene_mut = pd.read_csv(f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv", sep=",")
        ind_geneMut = data_gene_mut["sampleID"].values
        ind_geneMut = [name[:6] for name in ind_geneMut]

        ind_geneNonMut = []
        samples = self.data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample"]]
        ind_beataml = []

        for _, sample in samples.iterrows():
            if pd.isna(sample["dbgap_dnaseq_sample"]):
                sample_id = sample["dbgap_rnaseq_sample"][:6]
            else:
                sample_id = sample["dbgap_dnaseq_sample"][:6]
            ind_beataml.append(sample_id)

        for ind in ind_beataml:
            if ind not in ind_geneMut:
                ind_geneNonMut.append(ind)

        ind_geneNonMut_new = []
        for ind in ind_geneNonMut:
            if ind in self.index_list:
                ind_geneNonMut_new.append(ind)

        ind_geneNonMut = ind_geneNonMut_new


        if feature in ["TP53", "RUNX1", "ASXL1"]:
            gene_cat = ["Positive", "Negative Or Nan"]
            geneMutAndCat = self.get_number_for_bar_spe(gene_cat, ind_beataml, ind_geneMut, feature)
            geneNonMutAndCat = self.get_number_for_bar_spe(gene_cat, ind_beataml, ind_geneNonMut, feature)
        else:
            gene_cat = np.unique(list((self.data_beat_aml[feature])))
            gene_cat = [cat for cat in gene_cat if not (pd.isna(cat)) or cat != 'nan']
            geneMutAndCat = self.get_number_for_bar(gene_cat, ind_beataml, ind_geneMut, feature)
            geneNonMutAndCat = self.get_number_for_bar(gene_cat, ind_beataml, ind_geneNonMut, feature)

        L_ind = self.rearrangeZeros(geneMutAndCat, geneNonMutAndCat)
        CatSupprimees = [gene_cat[i] for i in L_ind]
        gene_cat = [gene_cat[i] for i in range(len(gene_cat)) if i not in L_ind] #on enlève les catégories dont la colonne est remplie de 0
        geneMutAndCat = [geneMutAndCat[i] for i in range(len(geneMutAndCat)) if i not in L_ind]
        geneNonMutAndCat = [geneNonMutAndCat[i] for i in range(len(geneNonMutAndCat)) if i not in L_ind]

        self.plot_graph_without_abundance(gene_cat, geneMutAndCat, geneNonMutAndCat)

        if len(CatSupprimees) > 0:
            messagebox.showwarning("Attention", f"Les catégories suivantes ont été supprimées car elles étaient vides: {', '.join(CatSupprimees)}")

        p, residus = self.Chi2Test(geneMutAndCat, geneNonMutAndCat)
        self.p_value_label.config(text=f"Résultat du test \u03C72 : p-value = {p:.5f}")

        alpha = 0.05
        if p < alpha:
            self.significance_label.config(text="La différence est significative")
        else:
            self.significance_label.config(text="La différence n'est pas significative")

        if self.SaveBool.get():
            self.CreateFileResFeat(ind_beataml, gene, feature, gene_cat, geneMutAndCat, geneNonMutAndCat, p, ind_geneMut, ind_geneNonMut)

        return
    

    def NormalizationByHousekeepingGene(self, ind_list, data_abundance):

        HousekeepingGenes = ["TBP", "ABL1", "PPIA"]
        # files = ["Documents/ScriptsPrincipaux/GeneMenage/query_results_TBP.tsv", "Documents/ScriptsPrincipaux/GeneMenage/query_results_ABL1.tsv", "Documents/ScriptsPrincipaux/GeneMenage/query_results_PPIA.tsv"]

        files = [f"Documents/ScriptsPrincipaux/GeneMenage/query_results_{gene}.tsv" for gene in HousekeepingGenes]

        NormalizedExpression = pd.DataFrame({
            'ID Sample': ind_list,
            'Abundance': data_abundance,
        })
        
        for i, file in enumerate(files):
            ExpressionsKmers = pd.read_csv(file, sep="\t")
            ExpressionsKmers = ExpressionsKmers.drop(columns=["seq_name"])
            ExpressionsKmers.columns = [col[:6] for col in ExpressionsKmers.columns]
            ExpressionsKmers = ExpressionsKmers.mean()
            ExpressionsKmers = ExpressionsKmers.loc[ind_list]
            ExpressionsKmers = ExpressionsKmers.values

            NormalizedExpression[f"HousekeepingGene_{HousekeepingGenes[i]}"] = ExpressionsKmers

        NormalizedExpression["HousekeepingMean"] = [0 for _ in range(NormalizedExpression.shape[0])]
        for index, row in NormalizedExpression.iterrows():
            NormalizedExpression.at[index, "HousekeepingMean"] = np.mean([row["HousekeepingGene_" + gene] for gene in HousekeepingGenes]) #moyenne arithmétique des gènes de ménage

        NormalizedExpression["NormalizedExpression"] = NormalizedExpression["Abundance"] / NormalizedExpression["HousekeepingMean"]
        
        NormalizedExpression.set_index("ID Sample", inplace=True)

        return NormalizedExpression
    

    def getFeatForPlotAbundance(self, NormalizedExpression, feature, ind_beataml):
        NormalizedExpression[feature] = ""

        for ind, _ in NormalizedExpression.iterrows():
            featValue = self.data_beat_aml.loc[ind_beataml.index(ind)][feature]
            NormalizedExpression.at[ind, feature] = featValue

        return NormalizedExpression
    

    def MannWhitneyUTest(self, Df, feature):
        cat_name = np.unique(list(Df[feature]))
        if self.AbondanceVar.get():
            group1 = Df[Df[feature] == cat_name[0]]["NormalizedExpression"].values
            group2 = Df[Df[feature] == cat_name[1]]["NormalizedExpression"].values
        elif self.DistributionVar.get():
            group1 = Df[Df[feature] == cat_name[0]]["%ratio"].values
            group2 = Df[Df[feature] == cat_name[1]]["%ratio"].values

        _, p = stats.mannwhitneyu(group1, group2)

        return p
    

    def ANOVATest(self, Df, feature):
        if self.AbondanceVar.get():
            _, p = stats.f_oneway(*[Df[Df[feature] == cat]["NormalizedExpression"].values for cat in Df[feature].unique()])
        elif self.DistributionVar.get():
            _, p = stats.f_oneway(*[Df[Df[feature] == cat]["%ratio"].values for cat in Df[feature].unique()])
        return p
    

    def plot_graph_with_abundance(self, NormalizedExpressionAndFeat):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)

        NormalizedExpressionAndFeat = NormalizedExpressionAndFeat.reset_index(drop=True)

        sns.boxplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, ax=self.ax) #showfliers=False
        sns.stripplot(x=feature, y='NormalizedExpression', data=NormalizedExpressionAndFeat, color='black', size=3, jitter=True, ax=self.ax)

        self.ax.set_title(f"Distribution de l'expression normalisée selon {feature} pour le gène {gene}")
        self.ax.set_xlabel(feature)
        self.ax.set_ylabel("Expression normalisée par TBP")

        self.fig.tight_layout()
        self.canvas.draw()
        return
    

    def CreateFileResAbund(self, NormalizedExpressionAndFeat, p, gene, feature):
        output_file = f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/Resultats_Abundance.txt"

        if not os.path.exists("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes"):
            os.mkdir("Documents/ScriptsPrincipaux/1_BarPlots/DossierRes")

        if len(np.unique(NormalizedExpressionAndFeat[self.feature_var.get()])) == 2:
            Test = "Mann-Whitney U"
        elif len(np.unique(NormalizedExpressionAndFeat[self.feature_var.get()])) > 2:
            Test = "ANOVA"
        
        with open(output_file, 'w') as f:
            f.write(f"#Résultats bruts de l'abondance de mutation du gène {self.gene_var.get()} selon la feature {self.feature_var.get()}, normalisée par le gène de ménage TBP.\n")
            f.write(f"#Test d'hypothèse de différence significative entre les moyennes des groupes\n#p-value = {p:.6f} avec le test {Test}\n")
            f.write(NormalizedExpressionAndFeat.to_string(index=False))

        self.fig.savefig(f"Documents/ScriptsPrincipaux/1_BarPlots/DossierRes/Abundance_{gene}_{feature}.png")

        print(f"Les résultats bruts ont été sauvegardés dans le fichier {output_file}")

        return
    

    def generate_plot_with_abundance(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrez un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrez une feature valide", f"La feature \"{feature}\" n'est pas valide")

        data_gene_mut = pd.read_csv(f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv", sep=",")
        ind_geneMut = data_gene_mut["sampleID"].values
        ind_geneMut = [name[:6] for name in ind_geneMut]
        data_abundance = data_gene_mut["mean_count_kmer_alt"].values

        NormalizedExpression = self.NormalizationByHousekeepingGene(ind_geneMut, data_abundance)


        samples = self.data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample"]]
        ind_beataml = []

        for _, sample in samples.iterrows():
            if pd.isna(sample["dbgap_dnaseq_sample"]):
                sample_id = sample["dbgap_rnaseq_sample"][:6]
            else:
                sample_id = sample["dbgap_dnaseq_sample"][:6]
            ind_beataml.append(sample_id)

        NormalizedExpressionAndFeat = self.getFeatForPlotAbundance(NormalizedExpression, feature, ind_beataml)
        print(NormalizedExpressionAndFeat)
        NormalizedExpressionAndFeat = NormalizedExpressionAndFeat[~NormalizedExpressionAndFeat[feature].isin(["Unknown", "UNKNOWN", "unknown", "nan"])]

        self.plot_graph_with_abundance(NormalizedExpressionAndFeat)

        if len(NormalizedExpressionAndFeat[feature].unique()) == 2:
            test = "Mann-Whitney U"
            p = self.MannWhitneyUTest(NormalizedExpressionAndFeat, feature)
        elif len(NormalizedExpressionAndFeat[feature].unique()) > 2:
            test = "ANOVA"
            p = self.ANOVATest(NormalizedExpressionAndFeat, feature)

        self.p_value_label.config(text=f"Test {test} : p-value = {p:.5f}")

        if self.SaveBool.get():
            self.CreateFileResAbund(NormalizedExpressionAndFeat, p, gene, feature)

        alpha = 0.05
        if p < alpha:
            self.significance_label.config(text="La différence est significative")
        else:
            self.significance_label.config(text="La différence n'est pas significative")      

        return



GUI()