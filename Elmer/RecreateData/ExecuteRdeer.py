#les deux index à utiliser : 
# BEAT-AML-tumor-474
# BEAT-AML-part2-224
# Script à modifier selon qu'on requête les kmers non mutés ou les kmers mutés

# rdeer query INDEX -q /chemin/vers/kmer_to_search_NPM1_varX_output/kmers.fa

import os


directory = "Scripts/RecreateData"

# Kmers de référence
dirs = os.listdir(f"{directory}/KmerToSearchRefOutput")

for dir in dirs:
    if os.path.exists(f"{directory}/KmerToSearchRefOutput/{dir}/kmers.fa"): 
        cmd1 = f"rdeer query BEAT-AML-tumor-474 -q {directory}/KmerToSearchRefOutput/{dir}/kmers.fa -o {directory}/RdeerRefOutput/{dir}_rdeer_part1_ref.tsv -p 12300"
        # print(cmd1)
        cmd2 = f"rdeer query BEAT-AML-part2-224 -q {directory}/KmerToSearchRefOutput/{dir}/kmers.fa -o {directory}/RdeerRefOutput/{dir}_rdeer_part2_ref.tsv -p 12300"
        os.system(cmd1)
        os.system(cmd2)


# Kmers mutés
dirs = os.listdir(f"{directory}/KmerToSearchOutput")

for dir in dirs:
    if os.path.exists(f"{directory}/KmerToSearchOutput/{dir}/kmers.fa"): 
        cmd1 = f"rdeer query BEAT-AML-tumor-474 -q {directory}/KmerToSearchOutput/{dir}/kmers.fa -o {directory}/RdeerOutput/{dir}_rdeer_part1_alt.tsv -p 12300"
        # print(cmd1)
        cmd2 = f"rdeer query BEAT-AML-part2-224 -q {directory}/KmerToSearchOutput/{dir}/kmers.fa -o {directory}/RdeerOutput/{dir}_rdeer_part2_alt.tsv -p 12300"
        os.system(cmd1)
        os.system(cmd2)