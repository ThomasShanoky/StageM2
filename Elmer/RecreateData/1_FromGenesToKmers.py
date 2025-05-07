import os
import gzip
import pandas as pd
from pprint import pprint
import shutil



Genes = ["NPM1", "FLT3", "DNMT3A", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"]
Genes.sort()

directory = '/'.join(os.path.abspath(__file__).split('/')[:-1])

annotation_ref_file = f"{directory}/Homo_sapiens.GRCh37.87.gtf.gz"
genome_ref_file = f"{directory}/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"


for gene in Genes:
    with gzip.open(annotation_ref_file, 'rt') as f:
        for line in f:
            if line[0] != '#':
                gene_name = line.split('"')[5]
                if gene_name == gene:
                    chrm = line.split()[0]
                    start_gene = int(line.split()[3])
                    end_gene = int(line.split()[4])

    WriteBool = False
    chr_seq = ""

    with gzip.open(genome_ref_file, 'rt') as f:
        for i, line in enumerate(f):
            if WriteBool and line[0] == ">": #on arrête dès qu'on est sur un autre chromosome
                WriteBool = False
            if WriteBool: #on écrit la séquence du chromosome
                chr_seq += line.strip()
            if line[0:1+len(chrm)] == ">"+chrm : #si on est au chromosome sur lequel est notre gène
                WriteBool = True

    Gene_seq_ref = chr_seq[start_gene-1:end_gene]

    Variants_file = f"{directory}/Full_AML_GeneVar/Full_AML_{gene}_variants.csv"
    Variants = pd.read_csv(Variants_file)[["chr", "start", "end", "ref", "alt", "in_cosmic"]]

    k = 31
    kmer_var = []
    kmer_non_mutated = []
    m = 0 
    for index, row in Variants.iterrows():
        m += 1 
        start = row["start"] - start_gene
        end = row["end"] - start_gene
        ref = row["ref"]
        alt = row["alt"]
        
        Gene_seq_mutated_list = list(Gene_seq_ref)

        if Gene_seq_mutated_list[start] != ref[:1]: #on ne prend que le premier nucléotide de la référence pour les délétions
            raise ValueError(f"Erreur : La référence est {ref[:1]} mais on a {Gene_seq_mutated_list[start]} à la position {start}.")
        
        Gene_seq_mutated_list[start:end+1] = list(alt)
        Gene_seq_mutated = "".join(Gene_seq_mutated_list)

        set_kmer = []
        set_kmer_non_mutated = []
        if len(alt) >= len(ref): #pour les insertions/substitution
            list_pos_mut = [i for i in range(start, end+len(alt))] 
        else:
            list_pos_mut = [i for i in range(start, start+len(alt))]

        for n in range(0, len(Gene_seq_mutated)-k+1): #on parcourt les positions possibles du gène
            start_kmer = n
            end_kmer = n+k

            t = 0
            for pos in list_pos_mut: #on vérifie que toutes les positions de la mutation sont bien comprises entre start et end
                if start_kmer <= pos < end_kmer:
                    t += 1

            if t == len(list_pos_mut):
                set_kmer.append(Gene_seq_mutated[start_kmer:end_kmer])
                set_kmer_non_mutated.append(Gene_seq_ref[start_kmer:end_kmer])

        kmer_var.append(set_kmer)
        kmer_non_mutated.append(set_kmer_non_mutated)


    output_files = [f"{directory}/KmerToSearch/kmer_to_search_{gene}_var{i}.fasta" for i in range(len(kmer_var))]
    for i, file in enumerate(output_files):
        with open(file, 'w') as f:
            for l, kmer in enumerate(kmer_var[i]):
                f.write(f">{gene}_var{i}_kmer{l}\n")
                f.write(f"{kmer}\n")

    for file in output_files:
        shutil.copy(file, "Elmer/scratch/users/tlouvet/Scripts/RecreateData/KmerToSearch/")
        # directement sauvegardés dans le dossier sur le serveur Elmer, pour executer kmerator et rdeer
        # Avant, j'executais ce script localement et je copiais les fichiers dans le dossier du serveur Elmer, car Elmer était lent à executer ce script (localement, ce script prend 20-30min)


    output_files = [f"{directory}/KmerToSearch/kmer_to_search_{gene}_ref{i}.fasta" for i in range(len(kmer_non_mutated))]
    for i, file in enumerate(output_files):
        with open(file, 'w') as f:
            for l, kmer in enumerate(kmer_non_mutated[i]):
                f.write(f">{gene}_ref{i}_kmer{l}\n")
                f.write(f"{kmer}\n")

    for file in output_files:
        shutil.copy(file, "Elmer/scratch/users/tlouvet/Scripts/RecreateData/KmerToSearchRef/")
        # directement sauvegardés dans le dossier sur le serveur Elmer, pour executer kmerator et rdeer


