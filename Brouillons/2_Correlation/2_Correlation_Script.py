import os



variants_files_dir = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAMLvariants"
variants_files = os.listdir(variants_files_dir)

file = variants_files[0]
print(file)

with open(f"{variants_files_dir}/{file}", "r") as f:
    for l, line in enumerate(f.readlines()):
        if l != 0:
            line_list = line.split(',')

            chrm = line_list[0]
            pos_start = line_list[1]
            pos_end = line_list[2]

            print(f"chr{chrm}:{pos_start}-{pos_end}")