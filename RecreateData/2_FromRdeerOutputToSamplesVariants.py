import pandas as pd
import numpy as np
import os 


# Seuil minimal de la somme des counts des kmers : en moyenne, chaque kmer doit être présent au moins 10 fois (compris)
seuil = 10 

Genes = ["NPM1", "FLT3", "DNMT3A", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"]
Genes.sort()


for gene in Genes:

    Gene_variants_file = f"Documents/RecreateData/Full_AML_GeneVar/Full_AML_{gene}_variants.csv"
    Gene_variants = pd.read_csv(Gene_variants_file)[["chr", "start", "end", "ref", "alt", "in_cosmic", "variant_classification", "gene", "total_freq"]]

    TableMutants = pd.DataFrame(columns=["sampleID", "gene", "ref", "alt", "localisation", "in_cosmic", "variant_classification", "sum_count_kmer_alt", "nb_kmers_alt", "mean_count_kmer_alt", "sum_count_kmer_ref", "nb_kmers_ref", "mean_count_kmer_ref", "ratio"])


    i = 0
    ResTheorEtPrat = []
    for index, row in Gene_variants.iterrows():
        ref = row["ref"]
        alt = row["alt"]
        in_cosmic = row["in_cosmic"]
        variant_classification = row["variant_classification"]
        gene = row["gene"]
        total_freq = row["total_freq"]
        localisation = "chr" + str(row["chr"]) + ":" + str(row["start"]) + "-" + str(row["end"])

        if in_cosmic == 1: #on ne prend que les variants dont on est sûr qu'ils sont oncogéniques

            file_alt1 = f"Documents/RecreateData/RdeerOutput/kmer_to_search_{gene}_var{i}_output_rdeer_part1_alt.tsv"
            file_alt2 = f"Documents/RecreateData/RdeerOutput/kmer_to_search_{gene}_var{i}_output_rdeer_part2_alt.tsv"
            if os.path.exists(file_alt1) and os.path.exists(file_alt2):
                df_alt1 = pd.read_csv(file_alt1, sep="\t")
                df_alt2 = pd.read_csv(file_alt2, sep="\t")

                df_alt1.set_index(df_alt1.columns[0], inplace=True)
                df_alt2.set_index(df_alt2.columns[0], inplace=True)

                df_alt = pd.concat([df_alt1, df_alt2], axis=1)

                nb_kmers_alt = df_alt.shape[0] - 2

                file_ref1 = f"Documents/RecreateData/RdeerRefOutput/kmer_to_search_{gene}_ref{i}_output_rdeer_part1_ref.tsv"
                file_ref2 = f"Documents/RecreateData/RdeerRefOutput/kmer_to_search_{gene}_ref{i}_output_rdeer_part2_ref.tsv"

                if os.path.exists(file_ref1) and os.path.exists(file_ref2):
                    df_ref1 = pd.read_csv(file_ref1, sep="\t")
                    df_ref2 = pd.read_csv(file_ref2, sep="\t")

                    df_ref1.set_index(df_ref1.columns[0], inplace=True)
                    df_ref2.set_index(df_ref2.columns[0], inplace=True)

                    df_ref = pd.concat([df_ref1, df_ref2], axis=1)

                    nb_kmers_ref = df_ref.shape[0] - 2

                n = 0
                for column in df_alt.columns:
                    sum_count_kmer_alt = sum(df_alt[column])
                    if sum_count_kmer_alt >= seuil*nb_kmers_alt:
                        n += 1
                        if os.path.exists(file_ref1) and os.path.exists(file_ref2):
                            sum_count_kmer_ref = sum(df_ref[column])
                            if sum_count_kmer_ref == 0 or nb_kmers_ref == 0:
                                sum_count_kmer_ref = np.nan
                                nb_kmers_ref = np.nan
                                mean_count_kmer_ref = np.nan
                                ratio = np.nan
                            else:
                                mean_count_kmer_ref = sum_count_kmer_ref / nb_kmers_ref
                                ratio = 100*(sum_count_kmer_alt/nb_kmers_alt)/((sum_count_kmer_alt/nb_kmers_alt)+mean_count_kmer_ref)
                        else:
                            sum_count_kmer_ref = np.nan
                            nb_kmers_ref = np.nan
                            mean_count_kmer_ref = np.nan
                            ratio = np.nan
                        TableMutants.loc[len(TableMutants)] = [column, gene, ref, alt, localisation, in_cosmic, variant_classification, sum_count_kmer_alt, nb_kmers_alt, sum_count_kmer_alt/nb_kmers_alt, sum_count_kmer_ref, nb_kmers_ref, mean_count_kmer_ref, ratio]
                ResTheorEtPrat.append((942*total_freq, n))

        i += 1

    TableMutants["sampleID"] = TableMutants["sampleID"].str[:6]
    TableMutants.to_csv(f"Documents/ScriptsPrincipaux/newMUTdata/{gene}_alt_perso.csv", index=False)
    # Sauvegardés dans le dossier des scripts principaux 

    print(f"{gene}")
    print(ResTheorEtPrat)