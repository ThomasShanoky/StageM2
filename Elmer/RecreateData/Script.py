import gzip
import time


start = time.time()

annotation_ref_file = "Scripts/RecreateData/Homo_sapiens.GRCh37.87.gtf.gz"
with gzip.open(annotation_ref_file, 'rt') as f:
    for line in f:
        if line[0] != '#':
            gene_name = line.split('"')[5]
            if gene_name == "NPM1":
                chrm = line.split()[0]
                start_gene = int(line.split()[3])
                end_gene = int(line.split()[4])

# print(chrm)
# print(start_gene)
# print(end_gene)


genome_ref_file = "Scripts/RecreateData/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
WriteBool = False
chr_seq = ""

with gzip.open(genome_ref_file, 'rt') as f:
    for line in f:
        if WriteBool and line[0] == ">": #on arrête dès qu'on est sur un autre chromosome
            WriteBool = False
        if WriteBool: #on écrit la séquence du chromosome
            chr_seq += line.strip()
        if line[0:2] == ">"+chrm : #si on est au chromosome sur lequel est notre gène
            # print(line.strip())
            WriteBool = True

print("OK")

NPM1_seq = chr_seq[start_gene-1:end_gene]
# print(len(NPM1_seq))


var_pos = 170818335
# var_pos = var_pos
var_pos_in_gene = var_pos - start_gene #dans start_gene, on prend deja en compte l'indexation 1-based 


k = 15 #pour faire des kmer de 31
kmer_alt = NPM1_seq[var_pos_in_gene-k:var_pos_in_gene] + NPM1_seq[var_pos_in_gene] + NPM1_seq[var_pos_in_gene+1:var_pos_in_gene+k+1]

kmer_to_search_file = "Scripts/RecreateData/kmer_to_search.fasta"
with open(kmer_to_search_file, 'w') as f:
    f.write(">kmer_alt\n")
    f.write(kmer_alt)


end = time.time()

print(f"Temps : {end-start}s")