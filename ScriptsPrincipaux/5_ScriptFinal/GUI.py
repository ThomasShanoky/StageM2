import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from DistributionFuncs import *



#Chargement des données
data_beat_aml = pd.read_csv("Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_Cliniques.csv", comment="#")

usable_cat = []
cats = list(data_beat_aml.columns)
for i, cat in enumerate(cats):
    if len(np.unique(data_beat_aml[cat].astype(str))) < 10 and i not in list(range(8)): #on ne prend pas les 7 premières features + on prend les catégories n'ayant que 10 ou moins de valeurs uniques
        usable_cat.append(cat)

index_file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_index.tsv"
with open(index_file) as f:
    index = f.readlines()[0]
index_list = index.split("\t")
index_list = [ind[:6] for ind in index_list]


# class GUI: