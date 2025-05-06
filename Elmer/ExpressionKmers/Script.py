import os
import pandas as pd
from pprint import pprint


# 1. Création d'un fichier Genes.tsv 
# 2. Kmerator pour générer les kmers de tous nos gènes
# 3. Rdeer sur les fichiers contigs.fa (sur les deux index de Beat-AML)
# 4. Merge les deux fichiers en un seul : CountKmers.tsv
# 5. Normalisation par le longueur des contigs + par le nombre total de kmers par échantillon

all_genes = "Scripts/ExpressionKmers/Genes.tsv"
contigs_sequence_file = "Scripts/ExpressionKmers/GenesKmers/contigs.fa"
contigs_count_file = "Scripts/ExpressionKmers/CountKmers.tsv"
kmers_per_sample = "Scripts/ExpressionKmers/TotalKmersPerSample.csv"

k = 31 #longueurs des kmers


def get_list_of_genes(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        genes = lines[0].split('\t')[:-1]
    return genes

Genes = get_list_of_genes(all_genes)
# print(Genes)


def contigs_lengths(file):
    length_dico = {}

    with open(file, 'r') as f:
        for l, line in enumerate(f.readlines()):
            if l%2 == 0: #ligne paire => en-tête de la séquence
                seq_name = line[1:][:-1][:49] #rdeer coupe le nom à 50 caractères
            else:
                seq = line[:-1]
                length_dico[seq_name] = len(seq)

    return length_dico

lengths_dico = contigs_lengths(contigs_sequence_file)
# pprint(lengths_dico)


def get_sample_list(file):

    with open(file, 'r') as f:
        lines = f.readlines()
        samples = lines[0].split("\t")
        samples = samples[1:] #on enlève "seq_name"
        samples = [samp[:6] for samp in samples] #on garde que les 6 premiers caractères

    return samples

SamplesList = get_sample_list(contigs_count_file)


def get_total_kmers_per_sample(file):
    kmers_per_sample = {}

    with open(file, 'r') as f:
        for l, line in enumerate(f.readlines()):
            if l != 0: #on enlève l'en-tête
                line_list = line.split(',')
                sample = line_list[0]
                val = float(line_list[1][:-1])
                kmers_per_sample[sample] = val

    return kmers_per_sample

TotKmers = get_total_kmers_per_sample(kmers_per_sample)
# pprint(TotKmers)



contigsCount = pd.read_csv(contigs_count_file, sep="\t")
contigsCount.set_index("seq_na", inplace=True)
contigsCount.index.name = None


##### Normalisation par le nombre de kmers dans chaque contig #####

for index, row in contigsCount.iterrows():
    den = lengths_dico[index]-k+1 #nombre de kmers dans un contig
    for sample in SamplesList:
        contigsCount.at[index, sample] = contigsCount.at[index, sample] / den # 1ère normalisation


##### Comptage moyen des contigs

new_contigsCount = pd.DataFrame(0, index=Genes, columns=contigsCount.columns)

genesComptabilises = []
# nb_contigs = 0
for index, row in contigsCount.iterrows():
    gene = index.split(':')[0]
    if gene not in genesComptabilises:
        if len(genesComptabilises) != 0:
            new_contigsCount.loc[genesComptabilises[-1], :] = new_contigsCount.loc[genesComptabilises[-1], :] / nb_contigs

        nb_contigs = 1
        new_contigsCount.loc[gene, :] = row
        genesComptabilises.append(gene)
    else:
        new_contigsCount.loc[gene, :] += row
        nb_contigs += 1
        
# print(new_contigsCount)

##### Normalisation par le nombre total de kmers dans chaque échantillon #####

# for sample in new_contigsCount:
#     expressions = new_contigsCount[sample]
#     for index, value in expressions.items():
#         new_contigsCount.at[index, sample] = (10**9) * new_contigsCount.at[index, sample] / TotKmers[sample] # 2ème normalisation

##### Enlever les ID samples qui ne sont pas indexés #####

def get_indexed_samples(file):
    indexed = []
    with open(file, 'r') as f:
        for l, line in enumerate(f.readlines()):
            if l != 0:
                line_list = line.split(',')
                sample = line_list[0]
                print(sample)
                indexed.append(sample)
    return indexed

indexed_samples = get_indexed_samples(kmers_per_sample)

for sample in new_contigsCount:
    if sample not in indexed_samples:
        new_contigsCount.drop(sample, axis=1, inplace=True)


print(new_contigsCount)
new_contigsCount.to_csv("Scripts/ExpressionKmers/ExpressionsWithKmers.csv", sep=",", index=True, header=True)

