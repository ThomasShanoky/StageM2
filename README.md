# Projet de stage : Développement d’un outil d’analyse de données issues de RNAseq en vue d’identifier des mutations associées à des paramètres cliniques de la leucémie myéloïde aigüe

## Description

Ce projet de stage a abouti au développement d'un outil possédant une interface graphique, permettant à l'utilisateur de sélectionner plusieurs paramètres pour executer différents types d'analyses sur des données (Beat-AML par défaut). 


## Structure du projet

    ScriptFinal/
        |
        ├── README.md
        ├── GUI.py
        ├── DistributionFuncs.py
        ├── AbundanceFuncs.py
        ├── FeaturesFuncs.py
        ├── MutationsFuncs.py 
        ├── BEATAML_Cliniques2.csv
        ├── BEATAML_Expressions.csv
        ├── ExpressionsWithKmers.csv
        ├── MUTdata/ #données de mutations 
        |       ├── GENE_alt_perso.csv *
        |       ├── ...
        |
        └── DossierRes/
                ├── N_ResFONCTIONNALITE_GENE_FEATURE/
                |     ├── Resultat_Bruts.txt
                |     ├── Resultat_plot.png
                |     └── Resultat_Table.tsv
                └── ...

*GUI.py* est le script principal lançant l'interface graphique. Les autres fichiers python (en *.py*) contiennent les fonctions nécessaires à tous les différents types d'analyses. 

## Technologies utilisées

Langage : Python
Panckages associés : 


### Cloner le dépôt

t


## Auteurs