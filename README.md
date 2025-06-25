# Projet de stage : Développement d’un outil d’analyse de données issues de RNAseq en vue d’identifier des mutations associées à des paramètres cliniques de la leucémie myéloïde aigüe

## Description

Ce projet de stage a abouti au développement d'un outil possédant une interface graphique, permettant à l'utilisateur de sélectionner plusieurs paramètres pour executer différents types d'analyses sur des données. Cet outil interactif d'exploration et de visualisation de données de patients est basé sur l'utilisation de k-mers (index des RNA-seq des cohortes de patients) et les métadonnées cliniques. Il offre une alternative modulable où les données sont fournies par l'utilisateur (Beat-AML par défaut). L'objectif est d'explorer ces données en menant plusieurs analyses (affichage graphique et tests statistiques) afin d'extraire des informations pertinentes pour lier les métadonnées et les événements de variation génétique.

L'outil propose six fonctionnalités analytiques, incluant des tests statistiques automatisés pour chacun, afin d'étudier les corrélations entre mutations géniques, expressions géniques et métadonnées. Appliqué à la cohorte Beat-AML (utilisée comme preuve de concept), il a permis de retrouver des associations déjà connues et d'identifier de nouvelles hypothèses biologiques. Cet outil est à visée exploratoire pour différents jeux de données et peut être étendu à d'autres pathologies autres que l'AML. Les données utilisables sont issues d'expériences de séquençage ARN (RNA-seq), permettant d'extraire des informations de mutations géniques (voire des variants géniques, incluant les variants d'épissage alternatif) ainsi que des métadonnées cliniques et phénotypiques des individus.


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