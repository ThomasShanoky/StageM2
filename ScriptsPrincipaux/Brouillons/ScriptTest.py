

tot_kmers_file_to_remake = "Documents/ScriptsPrincipaux/reindeer_matrix_eqc_info_beataml-474.txt"
tot_kmers_output_file = "Documents/ScriptsPrincipaux/TotalKmersPerSample.csv"


fr = open(tot_kmers_file_to_remake, 'r')
fw = open(tot_kmers_output_file, 'w')

fw.write("Sample,TotalKmers\n")


for line in fr.readlines():
    if line[0] != '#':
        sample = line.split(':')[0][:6]
        tot_kmers = line.split(':')[1].split()[0]

        fw.write(f"{sample},{tot_kmers}\n")