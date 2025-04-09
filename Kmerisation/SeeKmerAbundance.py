# Tracer un graphique où en abscisse on a les kmers (dans l'ordre d'apparition dans le fichier fasta) et en ordonnée le nombre d'occurences de chaque kmer pour chaque individus

# Bases d'index utilisées : Leucegene non stranded (203 experiments), Leucegene stranded (251 experiments)
#Count Method : Normalized

import matplotlib.pyplot as plt
import pandas as pd
import gc



def read_table(file:str):
    """Lis le fichier de résultat issu de Transipedia et le retourne sous forme de tableau"""

    return pd.read_csv(file, sep="\t")


def read_kmers_abundance_per_ind(tableau:pd.DataFrame):
    """A partir d'un tableau résultat issu d'une requête sur Transipedia, on extrait, pour chaque individu, le nombre d'occurences de chaque kmer.
    La sortie est un dictionnaire où chaque clé est le nom de l'individu et la valeur est une liste de nombres"""

    EchNames = list(tableau.columns)[1:] #on enlève le 'seq_name'
    
    KmersAbundance = {}

    for ech in EchNames:
        KmersAbundance[ech] = tableau[ech].values

    return KmersAbundance


def print_kmers_abundance(KmersAbundance:dict, grapheName:str):
    """Affiche un graphique où chaque tracé est un individu avec son nombre d'occurences pour chaque kmer"""

    for ech in KmersAbundance:
        plt.plot(KmersAbundance[ech], label=ech)
    isoName = grapheName.split("_")[1]
    plt.title("Abondance des K-mers pour chaque individu de la baseLeucegene\n(stranded et non-stranded) pour l'ARNm " + isoName)
    plt.ylabel("Occurences des K-mers")
    plt.xlabel("K-mers")
    plt.savefig(f"{grapheName}.png")
    # plt.legend()
    # plt.show()

    return


if __name__ == "__main__":

    Transipedia_files = ["query_results.tsv"]
    Transipedia_files = ["TransipediaResults/" + file for file in Transipedia_files]

    graphes_name = ["FLT3_transcrit_abundance"]

    for i, file in enumerate(Transipedia_files):
        Table = read_table(file)
        KmersPerInd = read_kmers_abundance_per_ind(Table)
        print_kmers_abundance(KmersPerInd, graphes_name[i])


#drastique baisse d'occurence des kmer (formant des rectangles) => mutation sur cet individu ?