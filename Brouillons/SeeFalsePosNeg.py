import numpy as np
import pandas as pd

beat_aml_file = "BEATAMLdata/BEATAML_Cliniques.csv"
mutation_file = "MUTdata/FLT3_alt.tsv"
data_beat_aml = pd.read_csv(beat_aml_file, sep=",")
data_mutation = pd.read_csv(mutation_file, sep="\t", comment="#")[["ID_sample", "type_alt", "seq_ref", "seq_alt", "kmer_ref", "kmer_alt"]]
data_mutation["ID_sample"] = [Id[:6] for Id in data_mutation["ID_sample"]]

FauxPos = ["BA2025", "BA2134", "BA2248", "BA2446", "BA2488", "BA2994", "BA3109", "BA3169", "BA3240", "BA3452", 
 "BA2031", "BA2106", "BA2123", "BA2226", "BA2226", "BA2236", "BA2499", "BA2502", "BA2512", "BA2670", 
 "BA2906", "BA2906", "BA2939", "BA3084", "BA3126", "BA3141", "BA3423", "BA2104", "BA2928", "BA3070", 
 "BA3071", "BA3072", "BA3369", "BA3381", "BA3381", "BA2054", "BA2455", "BA3362", "BA2066", "BA2386", 
 "BA2865", "BA3110", "BA3178", "BA2082", "BA2162", "BA3038", "BA2087", "BA2301", "BA2094", "BA2427", 
 "BA2427", "BA2456", "BA2561", "BA2631", "BA2631", "BA2879", "BA2879", "BA2930", "BA3447", "BA3447", 
 "BA2116", "BA2715", "BA2088", "BA2195", "BA2441", "BA2492", "BA2514", "BA2514", "BA2723", "BA2817", 
 "BA2866", "BA3202", "BA3318", "BA2789", "BA2792", "BA3222", "BA3353", "BA3374", "BA3403", "BA2101", 
 "BA2233", "BA2244", "BA2363", "BA2477", "BA3089", "BA2288", "BA2306", "BA2325", "BA2172", "BA2375", 
 "BA2905", "BA2489", "BA2524", "BA2344", "BA2577", "BA2748", "BA2592", "BA2714", "BA3138", "BA3229", 
 "BA3388", "BA2772", "BA3400", "BA3001", "BA3334", "BA3317", "BA3319", "BA3415", "BA3422", "BA3449"]

FauxNeg = []

samples = data_beat_aml[["dbgap_dnaseq_sample", "dbgap_rnaseq_sample"]]
ind_beataml = []

for _, sample in samples.iterrows():
    if pd.isna(sample["dbgap_dnaseq_sample"]):
        sample_id = sample["dbgap_rnaseq_sample"][:6]
    else:
        sample_id = sample["dbgap_dnaseq_sample"][:6]
    ind_beataml.append(sample_id)

NewCols = list(data_beat_aml.columns)
NewCols.insert(3, "kmer_alt")
NewCols.insert(3, "kmer_ref")
NewCols.insert(3, "seq_alt")
NewCols.insert(3, "seq_ref")
NewCols.insert(3, "type_alt")
NewCols.insert(3, "Classification")

# print(NewCols)

Res = pd.DataFrame(columns=NewCols)

for i, ind in enumerate(ind_beataml):
    if ind in FauxPos:
        ClassValue = "FauxPositif"
        type_altValue = data_mutation[data_mutation["ID_sample"] == ind]["type_alt"].values[0]
        seq_refValue = data_mutation[data_mutation["ID_sample"] == ind]["seq_ref"].values[0]
        seq_altValue = data_mutation[data_mutation["ID_sample"] == ind]["seq_alt"].values[0]
        kmer_refValue = data_mutation[data_mutation["ID_sample"] == ind]["kmer_ref"].values[0]
        kmer_altValue = data_mutation[data_mutation["ID_sample"] == ind]["kmer_alt"].values[0]
    elif ind in FauxNeg:
        ClassValue = "FauxNegatif"
        type_altValue = None
        seq_refValue = None
        seq_altValue = None
        kmer_refValue = None
        kmer_altValue = None

    if ind in FauxNeg or ind in FauxPos:
        # print(ind)
        row = data_beat_aml.iloc[i].copy()
        row["Classification"] = ClassValue
        row["type_alt"] = type_altValue
        row["seq_ref"] = seq_refValue
        row["seq_alt"] = seq_altValue
        row["kmer_ref"] = kmer_refValue
        row["kmer_alt"] = kmer_altValue
        Res = Res._append(row, ignore_index=True)

Res.to_csv("FLT3_MetaDataPourFauxPos-FauxNeg.csv", sep=",", index=False)