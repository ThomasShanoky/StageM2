#kmerator -f kmer_to_search.fasta -G 0
#si kmers mutés => -G 0 | si kmers non mutés => -G 1 (une seule fois présent)

# Script à modifier selon qu'on requête les kmers non mutés ou les kmers mutés
import os


directory = "Scripts/RecreateData"


# Kmers de référence
folder_kmers = f"{directory}/KmerToSearchRef"
files = os.listdir(folder_kmers)

for file in files:
    cmd = f"kmerator -f {folder_kmers}/{file} -G 1 -T 100 -o {directory}/KmerToSearchRefOutput/{file.split('.')[0]}_output"
    os.system(cmd)

    # os.remove(f"{folder_kmers}/{file}")


# Kmers mutés
folder_kmers = f"{directory}/KmerToSearch"
files = os.listdir(folder_kmers)

for file in files:
    cmd = f"kmerator -f {folder_kmers}/{file} -G 0 -T 0 -o {directory}/KmerToSearchOutput/{file.split('.')[0]}_output"
    os.system(cmd)

    # os.remove(f"{folder_kmers}/{file}")