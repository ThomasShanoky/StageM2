import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from DistributionFuncs import *
from AbundanceFuncs import *
from FeaturesFuncs import *
from MutationsFuncs import *



#Chargement des données
data_beat_aml = pd.read_csv("Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv", comment="#")
data_beat_aml = data_beat_aml[data_beat_aml["diseaseStageAtSpecimenCollection"] == "Initial Diagnosis"] #on ne prend que les patients ayant un diagnostic initial

usable_cat = []
cats = list(data_beat_aml.columns)
for i, cat in enumerate(cats):
    if 1 < len(np.unique(data_beat_aml[cat].astype(str))) < 10 and i not in list(range(8)): #on ne prend pas les 7 premières features + on prend les catégories n'ayant que 10 ou moins de valeurs uniques
        usable_cat.append(cat)

index_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_index.tsv"
with open(index_file) as f:
    index = f.readlines()[0]
index_list = index.split("\t")
index_list = [ind[:6] for ind in index_list]

Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] #prevenant de Leucegene : on prend les plus abondantes pour travailler sur un nombre limité de gènes (Sample count >= 10). Seul KMT2D n'est pas dans la liste des 140 gènes
Genes.sort()

data_mutation_files = [f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv" for gene in Genes]
data_mutation = [pd.read_csv(file, sep=",", comment="#")[["sampleID", "localisation", "ref", "alt"]] for file in data_mutation_files]

expressions_beataml_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_NormalizedExpression.csv"
data_expressions_beataml = pd.read_csv(expressions_beataml_file, sep=",", comment="#")
expressions_kmers_file = "Documents/ScriptsPrincipaux/NormalizedExpressionsWithKmers.csv"
data_expressions_kmers = pd.read_csv(expressions_kmers_file, sep=",", comment="#")



class GUI:

    def __init__(self, data_beat_aml=data_beat_aml, index_list=index_list, usable_cat=usable_cat, Genes=Genes, data_mutation=data_mutation, data_expressions_beataml=data_expressions_beataml, data_expressions_kmers=data_expressions_kmers):
        self.data_beat_aml = data_beat_aml
        self.index_list = index_list
        self.usable_cat = usable_cat
        self.Genes = Genes
        self.data_mutation = data_mutation
        self.data_expressions_beataml = data_expressions_beataml
        self.data_expressions = data_expressions_beataml #Expressions par défaut
        self.data_expressions_kmers = data_expressions_kmers

        self.window = tk.Tk()
        self.window.geometry("1400x800")
        self.window.config(bg='#87CEEB')
        # self.window.title("Analyse de l'effet de mutation sur différentes features")

        # Variables pour les listes déroulantes
        self.gene_var = tk.StringVar(value=self.Genes[0])
        self.feature_var = tk.StringVar(value=self.usable_cat[0])
        self.dico_IndAndMut, self.dico_mut = getTypesOfMutationsAndInd(self.Genes, self.Genes[0], self.data_mutation)
        self.mut_var = tk.StringVar(value=format_mutations(self.dico_mut, self.dico_IndAndMut)[0])

        # Listes déroulantes pour le gène et la feature
        self.gene_label = tk.Label(self.window, text="Sélectionnez un gène:")
        self.gene_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.gene_label.place(x=30, y=30)  
        self.gene_dropdown = ttk.Combobox(self.window, textvariable=self.gene_var, values=self.Genes)
        self.gene_dropdown.place(x=30, y=70)

        self.mut_label = tk.Label(self.window, text="Sélectionnez une mutation :")
        self.mut_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.mut_label.place(x=30, y=110)
        self.mut_dropdown = ttk.Combobox(self.window, textvariable=self.mut_var, values=format_mutations(self.dico_mut, self.dico_IndAndMut), width=45)
        self.mut_dropdown.place(x=30, y=150)

        self.feature_label = tk.Label(self.window, text="Sélectionnez une feature:")
        self.feature_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feature_label.place(x=30, y=190)  
        self.feature_dropdown = ttk.Combobox(self.window, textvariable=self.feature_var, values=self.usable_cat)
        self.feature_dropdown.place(x=30, y=230)

        #Choix de la provenance des données d'expression
        self.switch_state = False
        self.switch_button_expr = tk.Button(self.window, text="Expressions de BEAT AML", command=self.toggle_switch, font=("DejaVu Serif", 13), bg="#FF6347", fg="#FFFFFF", width=20)
        self.switch_button_expr.place(x=30, y=280)


        # Les 5 cases pour les différentes fonctionnalités
        self.DistributionVar = tk.BooleanVar()
        self.distribution_button = tk.Checkbutton(self.window, text="1. Montrer les distributions", variable=self.DistributionVar)
        self.distribution_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.distribution_button.place(x=30, y=350)  

        self.AbondanceVar = tk.BooleanVar()
        self.abondance_button = tk.Checkbutton(self.window, text="2. Montrer l'abondance des mutations", variable=self.AbondanceVar)
        self.abondance_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.abondance_button.place(x=30, y=380)

        self.FeatureVar = tk.BooleanVar(value=False)
        self.feat_button = tk.Checkbutton(self.window, text="3. Expression selon les features", variable=self.FeatureVar)
        self.feat_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feat_button.place(x=30, y=410)

        # Bouton pour générer le graphe
        self.generate_button = tk.Button(self.window, text="Générer le graphe", command=self.generate_plot)
        self.generate_button.config(width=15, font=("DejaVu Serif", 20, "bold"), highlightbackground="#370028", bg="#87CEEB", fg="#000000")
        self.generate_button.place(x=30, y=600)  

        # Canvas pour afficher le graphe
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(x=500, y=30, width=850, height=700)

        # Choix de sauvegarder les résultats
        self.SaveBool = tk.BooleanVar()
        self.save_button = tk.Checkbutton(self.window, text="Sauvegarder les résultats", variable=self.SaveBool)
        self.save_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.save_button.place(x=30, y=650)  

        # Labels pour afficher la p-value et la significativité
        self.p_value_label = tk.Label(self.window, text="")
        self.p_value_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.p_value_label.place(x=30, y=710)  

        self.significance_label = tk.Label(self.window, text="")
        self.significance_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.significance_label.place(x=30, y=750)  

        self.window.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.window.mainloop()
    

    def on_closing(self):
        self.window.quit()
        self.window.destroy()
        return
    

    def toggle_switch(self):
        """Change l'état du switch et met à jour le texte et la couleur."""
        self.switch_state = not self.switch_state
        if self.switch_state:
            self.switch_button_expr.config(text="Expression kmers", bg="#32CD32")
            self.data_expressions = data_expressions_kmers
        else:
            self.switch_button_expr.config(text="Expressions de BEAT AML", bg="#FF6347")
            self.data_expressions = data_expressions_beataml

        return 

    
    def generate_plot(self):
        #Vérification : une seule case cochée
        
        if sum([self.DistributionVar.get(), self.AbondanceVar.get(), self.FeatureVar.get()]) != 1:
            messagebox.showwarning("Attention", "Veuillez choisir une seule fonctionnalité")
            return 
        
        if self.DistributionVar.get():
            self.generate_plot_without_abundance()
            return
        if self.AbondanceVar.get():
            self.generate_plot_with_abundance()
            return
        if self.FeatureVar.get():
            self.generate_plot_feature()
            return


    def generate_plot_without_abundance(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrez un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrez une feature valide", f"La feature \"{feature}\" n'est pas valide")

        ind_geneNonMut = []
        samples = self.data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample"]]
        ind_beataml = []

        for _, sample in samples.iterrows():
            if pd.isna(sample["dbgap_dnaseq_sample"]):
                sample_id = sample["dbgap_rnaseq_sample"][:6]
            else:
                sample_id = sample["dbgap_dnaseq_sample"][:6]
            ind_beataml.append(sample_id)

        data_gene_mut = data_mutation[Genes.index(gene)]
        ind_geneMut = data_gene_mut["sampleID"].values
        ind_geneMut = [name[:6] for name in ind_geneMut if name[:6] in ind_beataml and name[:6] in self.index_list]

        for ind in ind_beataml:
            if ind not in ind_geneMut and ind in ind_beataml and ind in self.index_list:
                ind_geneNonMut.append(ind)

        # ind_geneNonMut_new = []
        # for ind in ind_geneNonMut:
        #     if ind in self.index_list:
        #         ind_geneNonMut_new.append(ind)

        # ind_geneNonMut = ind_geneNonMut_new

        gene_cat = np.unique(list((self.data_beat_aml[feature])))
        gene_cat = [cat for cat in gene_cat if not (pd.isna(cat)) or cat != 'nan']
        geneMutAndCat = get_number_for_bar(self.data_beat_aml, gene_cat, ind_beataml, ind_geneMut, feature)
        geneNonMutAndCat = get_number_for_bar(self.data_beat_aml, gene_cat, ind_beataml, ind_geneNonMut, feature)

        L_ind = rearrangeZeros(geneMutAndCat, geneNonMutAndCat)
        CatSupprimees = [gene_cat[i] for i in L_ind]
        gene_cat = [gene_cat[i] for i in range(len(gene_cat)) if i not in L_ind] #on enlève les catégories dont la colonne est remplie de 0
        geneMutAndCat = [geneMutAndCat[i] for i in range(len(geneMutAndCat)) if i not in L_ind]
        geneNonMutAndCat = [geneNonMutAndCat[i] for i in range(len(geneNonMutAndCat)) if i not in L_ind]

        self.canvas, self.fig, self.ax = plot_graph_without_abundance(self.fig, self.canvas, gene_cat, geneMutAndCat, geneNonMutAndCat, gene, feature)

        if len(CatSupprimees) > 0:
            messagebox.showwarning("Attention", f"Les catégories suivantes ont été supprimées car elles étaient vides: {', '.join(CatSupprimees)}")

        p, residus = Chi2Test(geneMutAndCat, geneNonMutAndCat)
        self.p_value_label.config(text=f"Résultat du test \u03C72 : p-value = {p:.5f}")

        alpha = 0.05
        if p < alpha:
            self.significance_label.config(text="La différence est significative")
        else:
            self.significance_label.config(text="La différence n'est pas significative")

        if self.SaveBool.get():
            CreateFileResFeat(data_beat_aml, ind_beataml, gene, feature, gene_cat, geneMutAndCat, geneNonMutAndCat, p, ind_geneMut, ind_geneNonMut, self.fig)

        return
    


    def generate_plot_with_abundance(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrez un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrez une feature valide", f"La feature \"{feature}\" n'est pas valide")

        samples = self.data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample"]]
        ind_beataml = []

        for _, sample in samples.iterrows():
            if pd.isna(sample["dbgap_dnaseq_sample"]):
                sample_id = sample["dbgap_rnaseq_sample"][:6]
            else:
                sample_id = sample["dbgap_dnaseq_sample"][:6]
            ind_beataml.append(sample_id)

        data_gene_mut = pd.read_csv(f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv", sep=",")
        data_gene_mut = data_gene_mut[data_gene_mut["sampleID"].isin(ind_beataml)]
        ind_geneMut = data_gene_mut["sampleID"].values
        ind_geneMut = [name[:6] for name in ind_geneMut if name[:6] in ind_beataml and name[:6] in self.index_list]
        data_abundance = data_gene_mut["mean_count_kmer_alt"].values

        NormalizedExpression = NormalizationByTotKmers(ind_geneMut, data_abundance)

        NormalizedExpression.set_index("ID Sample", inplace=True)
        NormalizedExpressionAndFeat = getFeatForPlotAbundance(data_beat_aml, NormalizedExpression, feature, ind_beataml)
        NormalizedExpressionAndFeat = NormalizedExpressionAndFeat[~NormalizedExpressionAndFeat[feature].isin(["Unknown", "UNKNOWN", "unknown", "nan"])]

        self.canvas, self.fig, self.ax = plot_graph_with_abundance(self.canvas, self.fig, NormalizedExpressionAndFeat, gene, feature)

        if len(NormalizedExpressionAndFeat[feature].unique()) == 2:
            test = "Mann-Whitney U"
            p = MannWhitneyUTest(NormalizedExpressionAndFeat, feature)
        elif len(NormalizedExpressionAndFeat[feature].unique()) > 2:
            test = "ANOVA"
            p = ANOVATest(NormalizedExpressionAndFeat, feature)

        self.p_value_label.config(text=f"Test {test} : p-value = {p:.5f}")

        if self.SaveBool.get():
            CreateFileResAbund(self.fig, NormalizedExpressionAndFeat, p, gene, feature)

        alpha = 0.05
        if p < alpha:
            self.significance_label.config(text="La différence est significative")
        else:
            self.significance_label.config(text="La différence n'est pas significative")      

        return
    


    def generate_plot_feature(self):
        gene = self.gene_var.get()
        feature = self.feature_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrer un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrer une feature valide", f"La feature \"{feature}\" n'est pas valide")

        featureValues = self.data_beat_aml[feature].dropna().unique().tolist()
        featureValues = [val for val in featureValues if val not in [np.nan, "Unknown", "unknown", "UNKNOWN"]]

        ##### Avoir les échantillons dont on peut extraire les features #####
        SamplesAndFeatures = getSamplesAndFeatures(self.data_beat_aml, feature, featureValues)

        ##### Avoir les échantillons dont on peut extraire les expressions #####
        ind_expr = list(self.data_expressions.columns)[1:]

        ##### Prendre l'intersection des deux listes #####
        inter_ind_pd = pd.DataFrame({"Sample": "", "Feature": "", "ExpressionGene": 0}, index=[0])

        for _, row in SamplesAndFeatures.iterrows():
            Id = row["Sample"]
            featureVal = row["Feature"]
            if Id in ind_expr:
                Expression = self.data_expressions.loc[self.data_expressions['Gene'] == gene, Id].values[0]
                inter_ind_pd = inter_ind_pd._append({"Sample": Id, "Feature": featureVal, "ExpressionGene": Expression}, ignore_index=True)

        inter_ind_pd = inter_ind_pd.drop(0).reset_index(drop=True)

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

        self.significance_label.config(text=stats_res_test+"\n"+significance)

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
            CreateFileRes(self.fig, inter_ind_pd, stats_res_test, premiereLigne, fileName)

        return
    

    ##### 4. Expressions des différents types de mutations #####

    def generate_plot_mutations(self):
    
        gene = self.gene_var.get()
        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrer un gène valide", f"Le gène \"{gene}\" n'est pas valide")

        dico_IndAndMut, dico_mut = getTypesOfMutationsAndInd(self.data_mutation, self.Genes, gene)

        AllIndMutated = list(dico_IndAndMut.values())
        AllIndMutated = [item for sublist in AllIndMutated for item in sublist] #à utiliser pour avoir les individus NON mutés

        IndToRemove = checkIfPatientsAreInExpressions(self.data_expressions, dico_IndAndMut)

        if len(IndToRemove) > 0:
            messagebox.showwarning("Individus non trouvés", f"Les individus suivants n'ont pas d'expressions de gène associées :\n{IndToRemove}")
            for ind in IndToRemove:
                for mut in dico_IndAndMut.keys():
                    if ind in dico_IndAndMut[mut]:
                        dico_IndAndMut[mut].remove(ind)

        dico_IndAndMut = filterDicoIndAndMut(dico_IndAndMut)

        Tableau = getTable(self.data_expressions, gene, dico_IndAndMut)

        NonMutatedInd = getNonMutatedInd(self.data_expressions, AllIndMutated)
        NonMutatedInd = checkIfIndExpressionAreInIndex(self.index_list, NonMutatedInd)

        Tableau = addNonMutatedIndExpression(self.data_expressions, Tableau, gene, NonMutatedInd)
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



GUI()