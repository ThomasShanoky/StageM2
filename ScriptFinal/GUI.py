import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import shutil
from DistributionFuncs import *
from AbundanceFuncs import *
from FeaturesFuncs import *
from MutationsFuncs import *



directory = '/'.join(os.path.abspath(__file__).split('/')[:-1])


data_beat_aml_file = f"{directory}/BEATAML_Cliniques.csv"
data_beat_aml = pd.read_csv(data_beat_aml_file, sep=",", comment="#") #données cliniques de Beat-AML (métadonnées)
data_beat_aml = data_beat_aml[data_beat_aml["diseaseStageAtSpecimenCollection"] == "Initial Diagnosis"] #on ne prend que les patients ayant un diagnostic initial
data_beat_aml.set_index("ID Sample", inplace=True)

usable_cat = []
cats = list(data_beat_aml.columns)
for i, cat in enumerate(cats):
    if 1 < len(np.unique(data_beat_aml[cat].astype(str))) < 10 and i not in list(range(8)): #on prend les catégories ayant entre 2 et 9 valeurs uniques + on ne prend pas les 7 premières features 
        usable_cat.append(cat)

index_file = f"{directory}/BEATAML_index.tsv"
with open(index_file) as f:
    index = f.readlines()[0]
index_list = index.split("\t")
index_list = [ind[:6] for ind in index_list] #Liste des ID Sample de tous les échantillons indexés (sur Transipedia)

Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] #provenant de Leucegene : on prend les plus abondants pour travailler sur un nombre limité de gènes (Sample count >= 10). Seul KMT2D, présent dans Leucegene, n'est pas dans la liste des 140 gènes donné par Stéphane/Sandra
Genes.sort()
print(len(Genes))

data_mutation_files = [f"{directory}/MUTdata/{gene}_alt_perso.csv" for gene in Genes]
data_mutation = [pd.read_csv(file, sep=",", comment="#")[["sampleID", "localisation", "ref", "alt", "mean_count_kmer_alt"]] for file in data_mutation_files] # données de mutations par gène : ID sample, la position de l'altération, la séquence référente et la séquence altérée

expressions_beataml_file = f"{directory}/BEATAML_Expressions.csv"
data_expressions_beataml = pd.read_csv(expressions_beataml_file, sep=",", comment="#") #données d'expressions directement prélevées de Beat-AML
expressions_kmers_file = f"{directory}/ExpressionsWithKmers.csv"
data_expressions_kmers = pd.read_csv(expressions_kmers_file, sep=",", comment="#") #données d'expressions construites avec les kmers uniques aux gènes

tot_kmers_file = f"{directory}/TotalKmersPerSample.csv" #nombre total de kmers par échantillon


##### Interface graphique #####

class GUI:

    def __init__(self, data_beat_aml=data_beat_aml, index_list=index_list, usable_cat=usable_cat, Genes=Genes, data_mutation=data_mutation, data_expressions_beataml=data_expressions_beataml, data_expressions_kmers=data_expressions_kmers, tot_kmers_file=tot_kmers_file, directory=directory):
        self.data_beat_aml = data_beat_aml
        self.index_list = index_list
        self.usable_cat = usable_cat
        self.Genes = Genes
        self.data_mutation = data_mutation
        self.data_expressions_beataml = data_expressions_beataml
        self.data_expressions = data_expressions_beataml #Expressions par défaut
        self.data_expressions_kmers = data_expressions_kmers
        self.tot_kmers_file = tot_kmers_file
        self.directory = directory

        self.SaveAll = False #Booléen pour sauvegarder tous les résultats significatifs

        # Fenêtre
        self.window = tk.Tk()
        self.window.geometry("1400x800")
        self.window.config(bg='#87CEEB')
        self.window.title("Fenêtre")

        # MEnu
        self.menu = tk.Menu(self.window)
        self.window.config(menu=self.menu)
        self.options_menu = tk.Menu(self.window, tearoff=0)
        self.menu.add_cascade(label="Options", menu=self.options_menu)
        self.options_menu.add_command(label="Vider le dossier des résultats", command=self.clear_res_folder)

        # Variables pour les listes déroulantes
        self.gene_var = tk.StringVar(value=self.Genes[0])
        self.feature_var = tk.StringVar(value=self.usable_cat[0])
        self.dico_IndAndMut, self.dico_mut = getTypesOfMutationsAndInd(self.Genes, self.Genes[0], self.data_mutation)
        self.mut_var = tk.StringVar(value=format_mutations(self.dico_mut, self.dico_IndAndMut)[0])
        # dico_IndAndMut contient en clefs les types de mutation (mut1, mut2 ...) et en valeurs les échantillons portant la mutation
        # dico_mut contient toutes les informations liées aux différents types de mutations d'un gène

        # Listes déroulantes pour le gène, la mutation et la feature
        self.gene_label = tk.Label(self.window, text="Sélectionnez un gène:")
        self.gene_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.gene_label.place(x=30, y=30)
        self.gene_exponent_label = tk.Label(self.window, text="(1 à 5)", font=("DejaVu Serif", 8), bg="#87CEEB", fg="#000000")
        self.gene_exponent_label.place(x=225, y=30) 
        self.gene_dropdown = ttk.Combobox(self.window, textvariable=self.gene_var, values=self.Genes, width=11)
        self.gene_dropdown.place(x=30, y=70)
        self.gene_dropdown.bind("<<ComboboxSelected>>", self.update_mutations) #permet d'associer l'événement "Sélection d'un gène" et la fonction update_mutations

        self.mut_label = tk.Label(self.window, text="Sélectionnez une mutation :")
        self.mut_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.mut_label.place(x=30, y=110)
        self.mut_exponent_label = tk.Label(self.window, text="(5)", font=("DejaVu Serif", 8), bg="#87CEEB", fg="#000000")
        self.mut_exponent_label.place(x=277, y=110)
        self.mut_dropdown = ttk.Combobox(self.window, textvariable=self.mut_var, values=format_mutations(self.dico_mut, self.dico_IndAndMut), width=45)
        self.mut_dropdown.place(x=30, y=150)

        self.feature_label = tk.Label(self.window, text="Sélectionnez une métadonnée:")
        self.feature_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feature_label.place(x=30, y=190)
        self.feature_exponent_label = tk.Label(self.window, text="(1, 2, 3 et 5)", font=("DejaVu Serif", 8), bg="#87CEEB", fg="#000000")
        self.feature_exponent_label.place(x=300, y=190)
        self.feature_dropdown = ttk.Combobox(self.window, textvariable=self.feature_var, values=self.usable_cat, width=30)
        self.feature_dropdown.place(x=30, y=230)

        # Choix de la provenance des données d'expression
        self.switch_state = False
        self.mut_button_expr = tk.Button(self.window, text="Expressions de BEAT AML", command=self.toggle_switch, font=("DejaVu Serif", 13), bg="#FF6347", fg="#FFFFFF", width=20)
        self.mut_button_expr.place(x=30, y=260)

        # Les 5 cases cochables pour les différentes fonctionnalités
        self.DistributionVar = tk.BooleanVar(value=False)
        self.distribution_button = tk.Checkbutton(self.window, text="1. Montrer les distributions", variable=self.DistributionVar, command=lambda: self.toggle_checkbuttons(self.DistributionVar))
        self.distribution_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.distribution_button.place(x=30, y=320)

        self.AbondanceVar = tk.BooleanVar(value=False)
        self.abondance_button = tk.Checkbutton(self.window, text="2. Montrer l'abondance des mutations", variable=self.AbondanceVar, command=lambda: self.toggle_checkbuttons(self.AbondanceVar))
        self.abondance_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.abondance_button.place(x=30, y=350)

        self.FeatureVar = tk.BooleanVar(value=False)
        self.feat_button = tk.Checkbutton(self.window, text="3. Expression selon la métadonnée", variable=self.FeatureVar, command=lambda: self.toggle_checkbuttons(self.FeatureVar))
        self.feat_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feat_button.place(x=30, y=380)

        self.MutationVar  = tk.BooleanVar(value=False)
        self.mut_button = tk.Checkbutton(self.window, text="4. Expression selon les mutations", variable=self.MutationVar, command=lambda: self.toggle_checkbuttons(self.MutationVar))
        self.mut_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.mut_button.place(x=30, y=410)

        self.FeatAndMutVar = tk.BooleanVar(value=False)
        self.mut_feat_button = tk.Checkbutton(self.window, text="5. Expression de cette mutation selon\nla métadonnée", variable=self.FeatAndMutVar, command=lambda: self.toggle_checkbuttons(self.FeatAndMutVar))
        self.mut_feat_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.mut_feat_button.place(x=30, y=440)

        # Bouton pour générer le graphe
        self.generate_button = tk.Button(self.window, text="Générer le graphe", command=self.generate_plot)
        self.generate_button.config(width=15, font=("DejaVu Serif", 20, "bold"), highlightbackground="#370028", bg="#87CEEB", fg="#000000")
        self.generate_button.place(x=30, y=500)  

        # Canvas pour afficher le graphe
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(x=530, y=30, width=850, height=700)

        # Choix de sauvegarder les résultats ou non
        self.SaveBool = tk.BooleanVar()
        self.save_button = tk.Checkbutton(self.window, text="Sauvegarder le prochain résultat", variable=self.SaveBool)
        self.save_button.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.save_button.place(x=30, y=550)

        # Labels pour afficher la p-value et la significativité
        self.p_value_label = tk.Label(self.window, text="")
        self.p_value_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.p_value_label.place(x=30, y=580)

        # Ligne séparatrice
        self.ligne = tk.Label(self.window, text="_______________________________________________________")
        self.ligne.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.ligne.place(x=15, y=620)

        # Choix du seuil de la p-value (alpha)
        self.label_alpha = tk.Label(self.window, text="Sélectionnez le seuil de la p-value (\u03B1):")
        self.label_alpha.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.label_alpha.place(x=30, y=650)
        self.alpha_var = tk.StringVar(value=0.05)
        self.alpha_entry = tk.Entry(self.window, textvariable=self.alpha_var)
        self.alpha_entry.config(font=("DejaVu Serif", 13), bg="#ffffff", fg="#000000", width=10)
        self.alpha_entry.place(x=375, y=650)

        # Optionnel : choix d'une feature/gène dont on veut sauvegarder les résultats significatifs
        self.choice_label = tk.Label(self.window, text="Choix d'une feature ou un gène :")
        self.choice_label.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.choice_label.place(x=30, y=690)

        self.feat_choice = tk.Label(self.window, text="Feature :")
        self.feat_choice.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.feat_choice.place(x=30, y=725)
        self.choice_feat_save = ["Toutes"] + self.usable_cat
        self.feat_choice_var = tk.StringVar(value=self.choice_feat_save[0])
        self.feat_choice_combobox = ttk.Combobox(self.window, textvariable=self.feat_choice_var, values=self.choice_feat_save, width=15)
        self.feat_choice_combobox.place(x=120, y=730)

        self.gene_choice = tk.Label(self.window, text="Gène :")
        self.gene_choice.config(font=("DejaVu Serif", 13), bg="#87CEEB", fg="#000000")
        self.gene_choice.place(x=270, y=725)
        self.choice_gene_save = ["Tous"] + self.Genes
        self.gene_choice_var = tk.StringVar(value=self.choice_gene_save[0])
        self.gene_choice_combobox = ttk.Combobox(self.window, textvariable=self.gene_choice_var, values=self.choice_gene_save, width=10)
        self.gene_choice_combobox.place(x=335, y=730)

        # Bouton pour sauvegarder tous les résultats significatifs
        self.generate_all_button = tk.Button(self.window, text="Générer les résultats significatifs", command=self.generate_all_results)
        self.generate_all_button.config(width=26, font=("DejaVu Serif", 14, "bold"), highlightbackground="#370028", bg="#87CEEB", fg="#000000")
        self.generate_all_button.place(x=30, y=760)

        # Lancement de la boucle principale
        self.window.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.window.mainloop()
    

    def on_closing(self):
        """Fonction appelée lorsque l'utilisateur ferme la fenêtre, permettant de quitter proprement l'application"""
        self.window.quit()
        self.window.destroy()
        return
    
    def toggle_checkbuttons(self, selected_var):
        """Permet de décocher les autres cases lorsqu'on en coche une"""
        checkbuttons = [self.DistributionVar, self.AbondanceVar, self.FeatureVar, self.MutationVar, self.FeatAndMutVar]
        for var in checkbuttons:
            if var != selected_var:
                var.set(False)
        return
    

    def clear_res_folder(self):
        """Vide le contenu du dossier de résultats DossierRes"""
        res_folder = f"{self.directory}/DossierRes"
        if os.path.exists(res_folder):
            all_sub_folders = os.listdir(res_folder)
            for fold in all_sub_folders:
                fold_path = f"{res_folder}/{fold}"
                shutil.rmtree(fold_path)
            messagebox.showinfo("Dossier vidé", "Le dossier de résultats a été vidé")
        return
    

    def update_mutations(self, event):
        """Met à jour la liste des mutations en fonction du gène sélectionné"""
        gene = self.gene_var.get()
        self.dico_IndAndMut, self.dico_mut = getTypesOfMutationsAndInd(self.Genes, gene, self.data_mutation)
        self.mut_var.set(format_mutations(self.dico_mut, self.dico_IndAndMut)[0])
        self.mut_dropdown['values'] = format_mutations(self.dico_mut, self.dico_IndAndMut)
        return
    

    def toggle_switch(self):
        """Change l'état du bouton de choix d'origine des expressions géniques et met à jour le texte et la couleur du bouton"""

        self.switch_state = not self.switch_state
        if self.switch_state:
            self.mut_button_expr.config(text="Expression kmers", bg="#32CD32")
            self.data_expressions = data_expressions_kmers
        else:
            self.mut_button_expr.config(text="Expressions de BEAT AML", bg="#FF6347")
            self.data_expressions = data_expressions_beataml

        return 

    
    def generate_plot(self):
        """Génère le graphe et le résultat statistique en fonction de la fonctionnalité sélectionnée"""

        if sum([self.DistributionVar.get(), self.AbondanceVar.get(), self.FeatureVar.get(), self.MutationVar.get(), self.FeatAndMutVar.get()]) != 1:
            messagebox.showwarning("Attention", "Veuillez choisir une fonctionnalité")
            return 
        
        try:
            self.alpha = float(self.alpha_var.get())
        except ValueError:
            messagebox.showwarning("Attention", "Veuillez entrer un nombre")
            return
        
        if self.alpha <= 0 or self.alpha >= 1:
            messagebox.showwarning("Attention", "Veuillez entrer un nombre entre 0 et 1 exclus")
            return
        
        gene = self.gene_var.get()
        feature = self.feature_var.get()
        Mutation = self.mut_var.get()

        if gene not in self.Genes:
            messagebox.showwarning("Veuillez entrez un gène valide", f"Le gène \"{gene}\" n'est pas valide")
        if feature not in self.usable_cat:
            messagebox.showwarning("Veuillez entrez une feature valide", f"La feature \"{feature}\" n'est pas valide")
        
        if self.DistributionVar.get():
            self.generate_plot_without_abundance(gene, feature)
            return
        if self.AbondanceVar.get():
            self.generate_plot_with_abundance(gene, feature)
            return
        if self.FeatureVar.get():
            self.generate_plot_feature(gene, feature, self.data_expressions)
            return
        if self.MutationVar.get():
            self.generate_plot_mutations(gene, self.data_expressions)
            return
        if self.FeatAndMutVar.get():
            if Mutation == "Toutes les mutations":
                messagebox.showwarning("Attention", "Veuillez sélectionner une seule mutation")
                return
            format_mutations(self.dico_mut, self.dico_IndAndMut)
            if Mutation not in self.mut_dropdown['values']:
                messagebox.showwarning("Attention", "Veuillez sélectionner une mutation valide")
            Mutation = Mutation.split(":")[0]
            self.generate_plot_feat_and_mut(gene, Mutation, feature, self.data_expressions)
            return


    ##### 1. Distribution des valeurs de la catégorie #####

    def generate_plot_without_abundance(self, gene, feature):

        ind_geneNonMut = []
        ind_beataml = self.data_beat_aml.index.tolist()

        data_gene_mut = data_mutation[Genes.index(gene)] #données de mutations du gène sélectionné (contenant les ID samples, la position de l'altération, la séquence référente et la séquence altérée)
        ind_geneMut = data_gene_mut["sampleID"].values
        ind_geneMut = [name[:6] for name in ind_geneMut if name[:6] in ind_beataml and name[:6] in self.index_list] #échantillons porteurs de la mutation

        for ind in ind_beataml:
            if ind not in ind_geneMut and ind in ind_beataml and ind in self.index_list:
                ind_geneNonMut.append(ind) #échantillons non porteurs de la mutation

        gene_cat = np.unique(list((self.data_beat_aml[feature])))
        gene_cat = [cat for cat in gene_cat if not (pd.isna(cat)) or cat != 'nan'] #Valeurs possibles de la feature sélectionnée (on enlève les NaN)
        geneMutAndCat = get_number_for_bar(self.data_beat_aml, gene_cat, ind_beataml, ind_geneMut, feature)
        geneNonMutAndCat = get_number_for_bar(self.data_beat_aml, gene_cat, ind_beataml, ind_geneNonMut, feature)

        L_ind = rearrangeZeros(geneMutAndCat, geneNonMutAndCat) #indices des colonnes (=catégories) dont la valeur est 0 pour les deux groupes (muté et non muté)
        CatSupprimees = [gene_cat[i] for i in L_ind]
        gene_cat = [gene_cat[i] for i in range(len(gene_cat)) if i not in L_ind] #on enlève les catégories dont la colonne est remplie de 0
        geneMutAndCat = [geneMutAndCat[i] for i in range(len(geneMutAndCat)) if i not in L_ind]
        geneNonMutAndCat = [geneNonMutAndCat[i] for i in range(len(geneNonMutAndCat)) if i not in L_ind]

        p = Chi2Test(geneMutAndCat, geneNonMutAndCat)
        self.p_value_label.config(text=f"Test \u03C72 : p-value = {p:.5f}")

        self.canvas, self.fig, self.ax = plot_graph_without_abundance(self.fig, self.canvas, gene_cat, geneMutAndCat, geneNonMutAndCat, gene, feature, p)

        if len(CatSupprimees) > 0 and not(self.SaveAll): #on n'affiche pas le message si on a demandé de sauvegarder tous les résultats significatifs
            messagebox.showwarning("Attention", f"Les catégories suivantes ont été supprimées car elles étaient vides: {', '.join(CatSupprimees)}")

        if p < self.alpha:
            SaveSignificant = True
        else:
            SaveSignificant = False

        if self.SaveBool.get() or (self.SaveAll and SaveSignificant):
            CreateFileResFeat(data_beat_aml, gene, feature, gene_cat, geneMutAndCat, geneNonMutAndCat, p, ind_geneMut, ind_geneNonMut, self.fig, self.SaveAll, self.directory)

        return
    

    ##### 2. Abondance des mutations #####

    def generate_plot_with_abundance(self, gene, feature):

        ind_beataml = self.data_beat_aml.index.tolist()

        # data_gene_mut = pd.read_csv(f"{self.directory}/newMUTdata/{gene}_alt_perso.csv", sep=",")
        data_gene_mut = data_mutation[self.Genes.index(gene)]
        data_gene_mut = data_gene_mut[data_gene_mut["sampleID"].isin(ind_beataml)] #filtration des échantillons non présent dans les échantillons Beat-AML déjà filtrés au début (qui sont donc qu'au diagnostic initial)
        ind_geneMut = data_gene_mut["sampleID"].values #échantillons porteurs de la mutation
        ind_geneMut = [name[:6] for name in ind_geneMut if name[:6] in ind_beataml and name[:6] in self.index_list]
        data_abundance = data_gene_mut["mean_count_kmer_alt"].values #abondance des mutations

        NormalizedExpression = NormalizationByTotKmers(ind_geneMut, data_abundance, self.tot_kmers_file)

        NormalizedExpression.set_index("ID Sample", inplace=True)
        NormalizedExpressionAndFeat = getFeatForPlotAbundance(data_beat_aml, NormalizedExpression, feature)
        NormalizedExpressionAndFeat = NormalizedExpressionAndFeat[~NormalizedExpressionAndFeat[feature].isin(["Unknown", "UNKNOWN", "unknown", "nan"])] # ~ = négation, on veut donc ici ENLEVER les valeurs de catégories de type "Unknown"

        if len(np.unique(list(NormalizedExpressionAndFeat[feature]))) == 1 and not(self.SaveAll): 
            messagebox.showwarning("Attention", "Il n'y a pas assez de valeurs de features")
            return 
        if len(np.unique(list(NormalizedExpressionAndFeat[feature]))) == 1 and self.SaveAll:
            return

        if len(NormalizedExpressionAndFeat[feature].unique()) == 2: #Comparaison de 2 moyennes
            test = "Mann-Whitney"
            p = MannWhitneyUTest(NormalizedExpressionAndFeat, feature)
        elif len(NormalizedExpressionAndFeat[feature].unique()) > 2: # Comparaison de plusieurs moyennes
            test = "ANOVA"
            p = ANOVATest(NormalizedExpressionAndFeat, feature)
        else:
            test = "Pas assez de features"
            p = 1

        self.canvas, self.fig, self.ax = plot_graph_with_abundance(self.canvas, self.fig, NormalizedExpressionAndFeat, gene, feature, p)
            
        self.p_value_label.config(text=f"Test {test} : p-value = {p:.5f}")

        if p < self.alpha:
            SaveSignificant = True
        else:
            SaveSignificant = False

        if self.SaveBool.get() or (self.SaveAll and SaveSignificant):
            CreateFileResAbund(self.fig, NormalizedExpressionAndFeat, test, p, gene, feature, self.SaveAll, self.directory)   

        return
    

    ##### 3. Expression selon les features #####

    def generate_plot_feature(self, gene, feature, expressions):

        featureValues = self.data_beat_aml[feature].dropna().unique().tolist()
        featureValues = [val for val in featureValues if val not in [np.nan, "Unknown", "unknown", "UNKNOWN"]] #valeurs uniques de la feature sélectionnée (on enlève les NaN et les Unknown)

        SamplesAndFeatures = getSamplesAndFeatures(self.data_beat_aml, feature, featureValues) #échantillons et features associées 

        ind_expr = list(expressions.columns)[1:] #échantillons qui ont une expression de gène

        inter_ind_pd = pd.DataFrame({"Sample": "", "Feature": "", "ExpressionGene": 0}, index=[0]) 

        for _, row in SamplesAndFeatures.iterrows():
            Id = row["Sample"]
            featureVal = row["Feature"]
            if Id in ind_expr:
                Expression = expressions.loc[expressions['Gene'] == gene, Id].values[0] #expression du gène pour l'échantillon
                inter_ind_pd = inter_ind_pd._append({"Sample": Id, "Feature": featureVal, "ExpressionGene": Expression}, ignore_index=True)

        inter_ind_pd = inter_ind_pd.drop(0).reset_index(drop=True) #suppression de la première ligne vide

        if len(inter_ind_pd["Feature"].unique()) > 2:
            #ANOVA
            groups = [inter_ind_pd[inter_ind_pd["Feature"] == featureVal]["ExpressionGene"] for featureVal in inter_ind_pd["Feature"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif len(inter_ind_pd["Feature"].unique()) == 2:
            #Mann-Whitney U 
            group1 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[0]]["ExpressionGene"]
            group2 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[1]]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif not(self.SaveAll):
            messagebox.showwarning("Pas assez de features", "Il n'y a pas assez de valeurs de features pour effectuer un test statistique")
            stats_res_test = ""
        else:
            p_val = 1
            stats_res_test = ""

        if p_val < self.alpha:
            SaveSignificant = True
        else:
            SaveSignificant = False

        self.p_value_label.config(text=stats_res_test)

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, ax=self.ax)
        sns.stripplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction de la feature {feature}")
        self.ax.set_ylabel("Expression du gène")
        self.ax.tick_params(axis='x', rotation=45)

        if p_val < self.alpha:
            if p_val < 0.0001:
                stars = "****"
            elif p_val < 0.001:
                stars = "***"
            elif p_val < 0.01:
                stars = "**"
            else:
                stars = "*"

            x1, x2 = 0, 1
            y, h, col = inter_ind_pd["ExpressionGene"].max() + 1, 1, "black"
            self.ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
            self.ax.text((x1 + x2) * 0.5, y + h, stars, ha="center", va="bottom", color=col)

        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get() or (self.SaveAll and SaveSignificant):
            premiereLigne = f"#Boxplots d'expressions du gène {gene} en fonction de la feature {feature}\n"
            fileName = f"3_FeatureOnly_{gene}_{feature}"
            CreateFileRes(self.directory, self.fig, inter_ind_pd, stats_res_test, premiereLigne, fileName, self.SaveAll)

        return
    

    ##### 4. Expressions des différents types de mutations #####

    def generate_plot_mutations(self, gene, expressions):

        dico_IndAndMut, _ = getTypesOfMutationsAndInd(self.Genes, gene, self.data_mutation)

        AllIndMutated = list(dico_IndAndMut.values())
        AllIndMutated = [item for sublist in AllIndMutated for item in sublist] #à utiliser pour avoir les échantillons NON mutés

        IndToRemove = checkIfPatientsAreInExpressions(expressions, dico_IndAndMut)

        if len(IndToRemove) > 0:
            for ind in IndToRemove:
                for mut in dico_IndAndMut.keys():
                    if ind in dico_IndAndMut[mut]:
                        dico_IndAndMut[mut].remove(ind)
            if not(self.SaveAll):
                messagebox.showwarning("Echantillons non trouvés", f"Les échantillons suivants n'ont pas d'expressions de gène associées :\n{IndToRemove}")

        dico_IndAndMut = filterDicoIndAndMut(dico_IndAndMut)

        Tableau = getTable(expressions, gene, dico_IndAndMut)

        NonMutatedInd = getNonMutatedInd(expressions, AllIndMutated)
        NonMutatedInd = checkIfIndExpressionAreInIndex(self.index_list, NonMutatedInd)

        Tableau = addNonMutatedIndExpression(expressions, Tableau, gene, NonMutatedInd)
        
        if len(dico_IndAndMut) > 1:
            # ANOVA
            groups = [Tableau[Tableau["TypeMutation"] == mut]["ExpressionGene"] for mut in Tableau["TypeMutation"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif len(dico_IndAndMut) == 1: #une seule mutation : on ne compare qu'avec les non mutés
            # Mann-Whitney U test
            #mutés
            group1 = Tableau[Tableau["TypeMutation"] == list(dico_IndAndMut.keys())[0]]["ExpressionGene"]
            #non mutés
            group2 = Tableau[Tableau["TypeMutation"] == "NonMut"]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif not(self.SaveAll):
            messagebox.showwarning("Pas assez de types de mutations", "Il n'y a pas assez de types de mutations pour effectuer un test statistique")
            stats_res_test = ""
        else:
            p_val = 1
            stats_res_test = ""

        if p_val < self.alpha:
            SaveSignificant = True
        else:
            SaveSignificant = False

        self.p_value_label.config(text=f"Test {stats_res_test}")

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="TypeMutation", y="ExpressionGene", data=Tableau, ax=self.ax)
        sns.stripplot(x="TypeMutation", y="ExpressionGene", data=Tableau, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction des types de mutations")
        self.ax.set_ylabel("Expression du gène")
        if p_val < self.alpha:
            if p_val < 0.0001:
                stars = "****"
            elif p_val < 0.001:
                stars = "***"
            elif p_val < 0.01:
                stars = "**"
            else:
                stars = "*"

            x1, x2 = 0, len(Tableau["TypeMutation"].unique()) - 1
            y, h, col = Tableau["ExpressionGene"].max() + 1, 1, "black"
            self.ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
            self.ax.text((x1 + x2) * 0.5, y + h, stars, ha="center", va="bottom", color=col)

        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get() or (self.SaveAll and SaveSignificant):
            premiereLigne = f"#Boxplots d'expressions du gène {gene} en fonction des types de mutations\n"
            fileName = f"4_MutOnly_{gene}"
            CreateFileRes(self.directory, self.fig, Tableau, stats_res_test, premiereLigne, fileName, self.SaveAll)
        
        return
    

    ##### 5. Expression de cette mutation selon la feature #####

    def generate_plot_feat_and_mut(self, gene, Mutation, feature, expressions):

        PatientsRelatedToMut = self.dico_IndAndMut[Mutation]

        featureValues = self.data_beat_aml[feature].dropna().unique().tolist()
        featureValues = [val for val in featureValues if val not in [np.nan, "Unknown", "unknown", "UNKNOWN"]]

        SamplesAndFeatures = getSamplesAndFeatures(self.data_beat_aml, feature, featureValues)

        inter_ind_pd = pd.DataFrame({"Sample": "", "Feature": "", "ExpressionGene": 0}, index=[0])

        for _, row in SamplesAndFeatures.iterrows():
            Id = row["Sample"]
            featureVal = row["Feature"]
            if Id in PatientsRelatedToMut and Id in expressions.columns:
                Expression = expressions.loc[expressions['Gene'] == gene, Id].values[0]
                inter_ind_pd = inter_ind_pd._append({"Sample": Id, "Feature": featureVal, "ExpressionGene": Expression}, ignore_index=True)

        inter_ind_pd = inter_ind_pd.drop(0).reset_index(drop=True)

        if len(inter_ind_pd["Feature"].unique()) == 1 and not(self.SaveAll):
            messagebox.showwarning("Pas assez de features", "Il n'y a pas assez de valeurs de features pour effectuer un test statistique")
            return
        if len(inter_ind_pd["Feature"].unique()) == 1 and self.SaveAll:
            return

        if len(inter_ind_pd["Feature"].unique()) > 2:
            # ANOVA
            groups = [inter_ind_pd[inter_ind_pd["Feature"] == featureVal]["ExpressionGene"] for featureVal in inter_ind_pd["Feature"].unique()]
            anova_result = stats.f_oneway(*groups)
            p_val = anova_result.pvalue
            stats_res_test = f"ANOVA : F = {anova_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif len(inter_ind_pd["Feature"].unique()) == 2:
            # Mann-Whitney U test
            group1 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[0]]["ExpressionGene"]
            group2 = inter_ind_pd[inter_ind_pd["Feature"] == list(inter_ind_pd["Feature"].unique())[1]]["ExpressionGene"]
            mannwhitney_result = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            p_val = mannwhitney_result.pvalue
            stats_res_test = f"Mann-Whitney : U = {mannwhitney_result.statistic:.3f}, p-value = {p_val:.5f}"
        elif not(self.SaveAll):
            messagebox.showwarning("Pas assez de features", "Il n'y a pas assez de valeurs de features pour effectuer un test statistique")
            stats_res_test = ""
        else:
            p_val = 1
            stats_res_test = ""

        if p_val < self.alpha:
            SaveSignificant = True
        else:
            SaveSignificant = False

        self.p_value_label.config(text=f"Test {stats_res_test}")

        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        sns.boxplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, ax=self.ax)
        sns.stripplot(x="Feature", y="ExpressionGene", data=inter_ind_pd, color='black', alpha=0.5, jitter=True, ax=self.ax)
        self.ax.set_title(f"Expression du gène {gene} en fonction de la feature {feature} pour la mutation {Mutation}")
        self.ax.set_ylabel("Expression du gène")
        self.ax.tick_params(axis='x', rotation=45)
        if p_val < self.alpha:
            if p_val < 0.0001:
                stars = "****"
            elif p_val < 0.001:
                stars = "***"
            elif p_val < 0.01:
                stars = "**"
            else:
                stars = "*"

            x1, x2 = 0, len(inter_ind_pd["Feature"].unique()) - 1
            y, h, col = inter_ind_pd["ExpressionGene"].max() + 1, 1, "black"
            self.ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
            self.ax.text((x1 + x2) * 0.5, y + h, stars, ha="center", va="bottom", color=col)
        self.fig.tight_layout()
        self.canvas.draw()

        if self.SaveBool.get() or (self.SaveAll and SaveSignificant):
            premiereLigne = f"#Boxplots d'expressions de la mutation {Mutation} ({self.dico_mut[Mutation][1]}>{self.dico_mut[Mutation][2]} aux positions {self.dico_mut[Mutation][0]}) du gène {gene}, en fonction de la feature {feature}\n"
            fileName = f"5_MutAndFeat_{gene}_{Mutation}_{feature}"
            CreateFileRes(self.directory, self.fig, inter_ind_pd, stats_res_test, premiereLigne, fileName, self.SaveAll)

        return
    

    ##### Génération de tous les résultats significatifs #####

    def generate_all_results(self):

        self.SaveAll = True

        try:
            self.alpha = float(self.alpha_var.get())
        except ValueError:
            messagebox.showwarning("Attention", "Veuillez entrer un nombre")
            return
        
        if self.alpha <= 0 or self.alpha >= 1:
            messagebox.showwarning("Attention", "Veuillez entrer un nombre entre 0 et 1 exclus")
            return
        
        if self.feat_choice_var.get() == "Toutes":
            list_features = self.usable_cat
        else:
            list_features = [self.feat_choice_var.get()]

        if self.gene_choice_var.get() == "Tous":
            list_genes = self.Genes
        else:
            list_genes = [self.gene_choice_var.get()]
        

        progress_window = tk.Toplevel(self.window)
        progress_window.title("Progression")
        progress_window.geometry("400x100")
        progress_window.resizable(False, False)

        progress_label = tk.Label(progress_window, text="Génération des résultats significatifs...")
        progress_label.place(x=10, y=10)

        progress_bar = ttk.Progressbar(progress_window, orient="horizontal", length=300, mode="determinate")
        progress_bar.place(x=10, y=40)

        tot_tasks = 2*len(list_genes)*len(list_features)
        progress_bar["maximum"] = tot_tasks
        progress_bar["value"] = 0

        for gene in list_genes:
            for expressions in [data_expressions_kmers, data_expressions_beataml]:
                self.generate_plot_mutations(gene, expressions)
                for feature in list_features:
                    progress_bar["value"] += 1
                    progress_window.update_idletasks()
                    self.generate_plot_feature(gene, feature, expressions)
                    self.generate_plot_with_abundance(gene, feature)
                    self.generate_plot_without_abundance(gene, feature)
                    for Mutation in format_mutations(self.dico_mut, self.dico_IndAndMut):
                        if Mutation != "Toutes les mutations":
                            Mutation = Mutation.split(":")[0]
                            self.generate_plot_feat_and_mut(gene, Mutation, feature, expressions)
        progress_window.destroy()
        messagebox.showinfo("Résultats générés", "Tous les résultats significatifs ont été générés et sauvegardés dans le dossier DossierRes")

        self.SaveAll = False

        return



GUI()