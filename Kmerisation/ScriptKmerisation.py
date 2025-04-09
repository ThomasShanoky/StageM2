#Idée générale : à partir de fichiers fasta de gènes, on veut un fichier qui contient ces gènes en k-mers (avec k = 31)



##################################
# Fonctions                      #
##################################

def read_fasta(FastaFile:str):
    """A partir d'un fichier fasta, on extrait les séquences et les noms des gènes"""

    SequencesList = []
    NamesList = [] #un fichier fasta contient plusieurs séquences (allèles ?)

    with open(FastaFile, 'r') as f:
        for line in f.readlines():
            if line[0] == '>': #on a un nouvel identifiant
                Name = line[1:].split(' (')[1].split(')')[0] 
                # print(Name)
                NamesList.append(Name)
                SequencesList.append("")
            else: #on continue la séquence
                SequencesList[-1] += line.strip() #strip => enlever le '\n'

    return NamesList, SequencesList


def kmerize(Sequence:str, k:int):
    """A partir d'une séquence, on extrait les k-mers"""
    kmers = []
    for i in range(len(Sequence) - k + 1): #ne pas dépasser la taille de la séquence
        kmers.append(Sequence[i:i+k])
    return kmers


def write_kmers(Names:list[str], Sequences:list[str], k:int, output_file:str):
    """On écrit les k-mers dans un fichier fasta"""

    with open(output_file, 'w') as newf:
        for index_name, seq in enumerate(Sequences):
            kmersList = kmerize(seq, k)
            for index_kmer, kmer in enumerate(kmersList):
                id_to_write = f">{Names[index_name]} kmer{index_kmer}\n"
                id_to_write = id_to_write.replace(" ", "_") #avoir un identifiant sans espace (nécessaire pour Transipedia)
                newf.write(id_to_write)
                newf.write(f"{kmer}\n")

    return



##################################
# Main : pour plusieurs fichiers #
##################################

def main():

    ### Pour plusieurs fichiers fasta ###

    fasta_files = ["NPM1_transcript_iso1.fasta", "NPM1_transcript_iso1bis.fasta", "NPM1_transcript_iso2.fasta", "NPM1_transcript_iso3.fasta", "NPM1_transcript_iso4.fasta", "NPM1_transcript_iso5.fasta", "NPM1_transcript_iso6.fasta", "NPM1_transcript_iso7.fasta"]

    output_files = ["FastaKmers/" + file_name[:-6] + "_kmers.fna" for file_name in fasta_files]

    fasta_files = ["FastaTranscriptFiles/" + file for file in fasta_files]

    k = 31
    
    
    for i in range(len(fasta_files)):
        Names, Sequences = read_fasta(fasta_files[i])
        write_kmers(Names, Sequences, k, output_files[i])

    return 



if __name__ == "__main__":
    main()