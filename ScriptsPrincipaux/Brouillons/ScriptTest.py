import pandas as pd



file = "Documents/ScriptsPrincipaux/BEATAMLdata/BEATAML_NormalizedExpression.tsv"
output_file = "Documents/ScriptsPrincipaux/BEATAMLdata/TestExpressions.csv"

Genes = ["NPM1", "DNMT3A", "FLT3", "TET2", "NRAS", "TP53", "RUNX1", "IDH2", "ASXL1", "WT1", "KRAS", "IDH1", "PTPN11", "SRSF2", "CEBPA", "KIT", "NF1", "STAG2", "GATA2", "EZH2", "BCOR", "JAK2", "SMC1A", "RAD21", "SF3B1", "CBL"] 

with open(file, 'r') as f:
    for l, line in enumerate(f.readlines()):
        line_list = line.strip().split("\t")
        if l == 0:
            header = line_list[4:]
            header = ["Gene"] + [sample[:6] for sample in header]
            expressions_df = pd.DataFrame(columns=header)
        else:
            gene = line_list[1]
            if gene in Genes:
                expressions = line_list[4:]
                expressions_df.loc[gene] = [gene] + expressions

# print(expressions_df)
expressions_df.set_index("Gene", inplace=True)
expressions_df.to_csv(output_file, sep=",", index=True)