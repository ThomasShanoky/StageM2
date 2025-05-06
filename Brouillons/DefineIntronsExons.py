import gzip



genome_ref_file = "Documents/RecreateData/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
annotations_file = "Documents/RecreateData/Homo_sapiens.GRCh37.87.gtf.gz"

geneToSearch = "ATF4"


l = 0
with gzip.open(annotations_file, 'rt') as f:
    for line in f.readlines():
        columns = line.strip().split('\t')
        if line[0] != '#' and f'gene_name "{geneToSearch}"' in line and columns[2] == "exon":
            l += 1
            print(line)
print(l)
