# Projet de stage : Développement d’un outil d’analyse de données issues de RNAseq en vue d’identifier des mutations associées à des paramètres cliniques de la leucémie myéloïde aigüe

## Description du projet

Ce projet de stage a abouti au développement d'un outil possédant une interface graphique, permettant à l'utilisateur de sélectionner plusieurs paramètres pour executer différents types d'analyses sur des données. Cet outil interactif d'exploration et de visualisation de données de patients est basé sur l'utilisation de k-mers (index des RNA-seq des cohortes de patients) et les métadonnées cliniques. Il offre une alternative modulable où les données sont fournies par l'utilisateur (Beat-AML par défaut). L'objectif est d'explorer ces données en menant plusieurs analyses (affichage graphique et tests statistiques) afin d'extraire des informations pertinentes pour lier les métadonnées et les événements de variation génétique.

L'outil propose six fonctionnalités analytiques, incluant des tests statistiques automatisés pour chacun, afin d'étudier les corrélations entre mutations géniques, expressions géniques et métadonnées. Appliqué à la cohorte Beat-AML (utilisée comme preuve de concept), il a permis de retrouver des associations déjà connues et d'identifier de nouvelles hypothèses biologiques. Cet outil est à visée exploratoire pour différents jeux de données et peut être étendu à d'autres pathologies autres que l'AML. Les données utilisables sont issues d'expériences de séquençage ARN (RNA-seq), permettant d'extraire des informations de mutations géniques (voire des variants géniques, incluant les variants d'épissage alternatif) ainsi que des métadonnées cliniques et phénotypiques des individus.



### Technologies utilisées

Langage : python (version 3.10.12)
Panckages associés : 
- numpy (version 2.1.3)
- pandas (version 2.2.3)
- tkinter (version 8.6)
- matplotlib (version 3.10.0)
- seaborn (version 0.13.2)
- scipy (version 1.15.2)
- shutil
- os



### Structure du projet

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



### Explication des différentes fonctionnalités

Il y a 6 fonctionnalités executable via l'interface graphique.

1. **Comparaison mutation / métadonnée**  
   L’utilisateur choisit un gène et une métadonnée. Les échantillons sont répartis en deux groupes : mutés et non mutés pour ce gène. Pour chaque valeur de la métadonnée, une barre est affichée dans un diagramme en barre. Une table de contingence est construite (avec la présence de mutation en ligne, et la métadonnée en colonne) afin de visualiser la distribution et effectuer un test du $\chi^2$ pour vérifier s'il existe une association significative entre la présence de mutation sur un gène et la distribution des catégories d'une métadonnée.

2. **Boxplots d’abondance de mutations selon une métadonnée**  
   Pour un gène donné, l’outil prend les échantillons mutés, extrait leur abondance de mutation et les normalise par le nombre total de k-mers. Ensuite, il classe ces abondances selon les valeurs de la métadonnée sélectionnée. Des boxplots (ou boite à moustache) sont générés et affichés pour chaque valeur possible, et une p-valeur est calculée pour vérifier la significativité de la différence entre groupes. 

3. **Boxplots d’expression génique selon une métadonnée**  
   Fonctionnalité similaire à (2), mais ici ce sont les valeurs d’expression génique qui sont utilisées (données par l'utilisateur ou basées sur les k-mers). Une table est construite avec les identifiants d’échantillons, les valeurs de la métadonnée, et les expressions géniques. Le script génère des boxplots comparant les expressions selon la métadonnée. Un test statistique est également effectué. 

4. **Effet des variants d’un gène sur son expression**  
   Après avoir sélectionné un gène, l’outil identifie les différents variants observés (définis par leur position et séquences altérée et de référence) chez les patients. Chaque variant est associé à un groupe d’échantillons. Les patients sans valeur d’expression sont exclus. Les variants représentés par moins de 4 individus sont filtrés. Une table est créée avec les expressions par variant (et les non-mutés sont ajoutés). L’expression du gène est comparée entre variants via des boxplots et un test statistique est executé.

5. **Effet d’un variant précis sur une métadonnée**  
   L’utilisateur choisit un variant spécifique d’un gène. L’outil récupère les échantillons concernés par ce variant et associe leurs valeurs de métadonnée. Il trace ensuite un boxplot pour visualiser l’effet du variant sur la métadonnée par son expression. Les variants doivent être représentés par au moins 4 individus. Aucun graphe n’est tracé si une seule valeur de métadonnée est présente. Ce test est exploratoire, car souvent peu d’échantillons.

6. **Corrélation entre expressions de deux gènes**  
   L’utilisateur sélectionne deux gènes. Le script récupère les expressions correspondantes pour chaque échantillon et construit une table commune. Une corrélation de Pearson est calculée, et le nuage de points est affiché avec les expressions des deux gènes sur les axes. Aucun test statistique n’est réalisé, seule la valeur du coefficient $r$ est donnée. 


L'utilisateur peut faire le choix, au lieu de voir les résultats un par un, de sauvegarder un ensemble de résultats significatifs. Il peut choisir le seuil de significativité (par défaut $p=0.05$) pour les fonctionnalités 1 à 5, et un seuil de coefficient de corrélation ($|r|=0.5$) pour la dernière fonction. \
De plus, s'il a un gène et/ou une métadonnée d'intérêt, il peut choisir de sauvegarder un ensemble de résultats significatifs en lien avec ce gène et/ou cette métadonnée. 



### Guide d'utilisation

Par défaut, l'outil utilise les données de Beat-AML : les métadonnées, les expressions et les mutations. L'outil est déjà prêt à l'emploi pour ces données pour le tester. \
Si l'utilisateur veut executer les différentes fonctionnalités, il doit remplacer les données de Beat-AML par les siennes. 


- Données de mutation : 
L'utilisateur doit posséder un fichier nommé `data_mutations.vcf` avec les colonnes suivantes au minimum : `CHROM` indiquant le chromosome de la mutation, `POS` donnant la position génomique du variant, `ID` qui est l'identifiant d'échantillon portant le variant, `GENE` qui est le gène où est présent le variant, et enfin `REF` et `ALT` qui sont les séquences de référence et altérée. Le fichier peut bénéficier d'une autre colonne, `ABUNDANCE`, qui contient les valeurs d'abondance en k-mers du variant pour l'échantillon. Cette colonne est optionnel pour l'execution de l'outil mais indispensable si l'utilisateur souhaite générer les résultats de la deuxième fonctionnalité. \
Il est à noté que la colonne `ID` n'est pas une colonne d'identifiant pour les variants mais bien pour les **échantillons**, nous permettant de lier les échantillons aux variants qu'ils portent. 

- Nombre total de k-mers par échantillons :
Si l'utilisateur les possède, il peut fournir dans un fichier `TotalKmersPerSample.csv` le nombre total de k-mers présent dans chaque échantillon. La première colonne doit être les identifiants d'échantillon avec comme header `ID` et la deuxième colonne doit être les nombres totaux de k-mers avec en header `TotalKmers`. Ce fichier n'est pas nécessaire à l'execution de l'outil principal (interface graphique), il permet de normaliser les données d'abondances précédemments

- Données d'expression en RPKM :
L'utilisateur peut fournir des données d'expression déjà prêtes s'il en possède. Le nom du fichier doit être `ExpressionsDonnees.csv`, et il doit contenir : 
  - En colonne, les gènes. Le nom du header de la colonne doit obligatoirement être `Gene`
  - En ligne, les identifiants des échantillons. Ces derniers doivent faire la même taille de caractères (Par exemple pour Beat-AML, les identifiants sont composés de 6 caractères).

- Données d'expression via les k-mers : 
L'utilisateur doit avoir les données d'expression par k-mers des gènes qu'il veut étudier. S'il ne les a pas, il doit posséder un environnement dans lequel un ou des index de ses données RNA-seq existe (obtenable par exemple, avec REINDEER), il peut ensuite suivre ces étapes : 
  1. Créer un fichier `Genes.tsv` contenant l'ensemble de ses gènes.
  2. Récupérer les contigs spécifiques de ces gènes (par exemple, à l'aide de l'outil Kmerator qui lui donnera l'ensemble des contigs spécifiques aux gènes).
  3. Récupérer les abondances de ces contigs dans l'index dans un dossier nommé `GenesKmers`.
  4. Si la requête s'est effectuée sur plusieurs index, il faut regrouper les comptages en un seul fichier, nommé `CountKmers.tsv`.
  5. Executer le script présent dans le dossier `Elmer/ExpressionKmers` avec le fichier `CountKmers.tsv` présent dans ce même dossier, avec `Genes.tsv`, le dossier `GenesKmers` et enfin un fichier `SampleID.csv` qui contient tous les identifiants d'échantillons avec comme header `Sample`.
  6. L'utilisateur peut copier/coller ou déplacer le fichier final `ExpressionsWithKmers.csv` dans le dossier `ScriptFinal/`.

  Attention, si les expressions d'un gène sont indisponibles dans un fichier d'expression, il doit être supprimé du fichier VCF.

- Métadonnées : 
Un fichier sous le nom de `Metadata.csv` contenant les métadonnées phénotypiques et cliniques des patients doit être présent. Il doit obligatoirement comporter la colonne `ID Sample` (avec ce même nom d'header) contenant les identifiants d'échantillons. Il est à noté que les échantillons doivent être uniques aux patients (chaque patient doit posséder exactement un échantillon RNA-seq) afin d'assurer l'indépendance des échantillons lors des tests statistiques effectués en aval. 


Une fois ces fichiers présents dans le dossier `ScriptFinal/`, l'interface graphique est prête à être affichée en executant le script `GUI.py`. Ce dernier utilise les fonctions écrites dans les autres fichiers python : `DistributionFuncs.py`, `AbundanceFuncs.py`, `FeaturesFuncs.py` et `MutationsFuncs.py`.

Le dossier `DossierRes` contient tous les résultats sauvegardés, de manière automatique via le bouton "Générer les résultats significatifs" ou manuellement par l'utilisateur via le bouton "Sauvegarder le prochain résultat".



## Cloner le dépôt

git clone https://github.com/ThomasShanoky/StageM2.git



## Auteurs

Thomas LOUVET\
thomaslouvet7@gmail.com

Stage encadré par Chloé BESSIERE, Sandra DAILHAU et Stéphane PYRONNET. 

