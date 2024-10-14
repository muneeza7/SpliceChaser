import os
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
import math

# Load libraries (assuming they are already installed)
from pandas import read_csv, read_table
from string import punctuation

# Set working directory
os.chdir("./")

results_dir = "./results/"
bams_dir = "./bams/tmp_bams/"
output_dir = "./bams/tmp_bams/output/"

bname = [os.path.basename(filename.replace("_result_final.csv", "")) for filename in os.listdir(results_dir) if filename.endswith("final.csv")]
print(bname)

for i in range(len(bname)):
    raw_dat_1 = pd.read_csv(results_dir + "/" + bname[i] + "_result_final.csv", header=None)
    raw_dat_1 = raw_dat_1.iloc[1:].reset_index(drop=True)
    dat = pd.read_table(bams_dir + "/" + bname[i] + ".star.b37.bam.all_jxs.tsv", header=None, dtype=str)
    
    id1 = bname[i]
    print("processing ...", id1)

    if len(raw_dat_1) <= 2:
        raw_dat_1.to_csv(f"{output_dir}/{id1}_Shannon_score.csv", index=False)
    else:
        for junc in range(len(raw_dat_1)):
            dat0 = dat[(dat[2] == raw_dat_1.iloc[junc, 1]) & (dat[3] == raw_dat_1.iloc[junc, 2]) & (dat[6] == "1")]

            # Use .loc to avoid SettingWithCopyWarning
            dat0 = dat0.copy()  # Ensure we're working on a copy

            dat0.loc[:, "n_count"] = dat0[5].str.count("N")
            dat0.loc[:, "first_m"] = dat0[5].str.extract(r"([0-9]+)M").astype(int)
            dat0.loc[:, "last_m"] = dat0[5].str.extractall(r"([0-9]+)M").groupby(level=0).last().astype(int)

            dat1 = dat0.groupby([1, 2, 3, "last_m", 6]).size().reset_index(name="n")
            dat1["frac"] = dat1["n"] / dat1["n"].sum()
            dat1["total"] = dat1["n"].sum()

            dat2 = dat0.groupby([1, 2, 3, "first_m", 6]).size().reset_index(name="m")
            dat2["frac2"] = dat2["m"] / dat2["m"].sum()
            dat2["total2"] = dat2["m"].sum()

            H_end = round(-(dat1["frac"] * dat1["frac"].apply(lambda x: math.log(x))).sum(), 2)
            H_start = round(-(dat2["frac2"] * dat2["frac2"].apply(lambda x: math.log(x))).sum(), 2)

            start = int(raw_dat_1.loc[junc, 1])
            stop = int(raw_dat_1.loc[junc, 2])

            # Update raw_dat_1 columns
            raw_dat_1.loc[junc, "H_start"] = H_start
            raw_dat_1.loc[junc, "fraction_start"] = "_".join(round(dat2["frac2"], 3).astype(str))
            raw_dat_1.loc[junc, "first_match"] = "_".join(dat2["first_m"].astype(str))
            raw_dat_1.loc[junc, "n_groups_start"] = "_".join(dat2["m"].astype(str))

            raw_dat_1.loc[junc, "H_end"] = H_end
            raw_dat_1.loc[junc, "fraction_end"] = "_".join(round(dat1["frac"], 3).astype(str))
            raw_dat_1.loc[junc, "last_match"] = "_".join(dat1["last_m"].astype(str))
            raw_dat_1.loc[junc, "n_groups_end"] = "_".join(dat1["n"].astype(str))

            raw_dat_1.loc[junc, "del_size"] = abs(start - stop) + 1
            raw_dat_1.loc[junc, "diff_H_start_end"] = abs(H_start - H_end)

        final = raw_dat_1.rename(
            columns={0: "chr", 1: "start", 2: "stop", 3: "startAnno", 4: "stopAnno", 5: "intron_motif", 6: "unique_count",
                     7: "multimapped_count", 8: "max_overhang", 9: "geneStrand", 10: "transcriptName", 11: "totalEx", 
                     12: "geneName", 13: "readStrand", 14: "prop_uniq_multi", 15: "first5", 16: "last5", 
                     17: "first3", 18: "last3", 19: "event5", 20: "event3", 21: "location5", 22: "location3", 
                     23: "strandDirection", 24: "event", 25: "nExSkip", 26: "exon", 27: "withinEx5", 28: "withinEx3", 
                     29: "nGene", 30: "nTranscripts", 31: "distance_from_last5", 32: "distance_from_first3", 
                     33: "Gene", 34: "GTExFraction", 35: "counts", 36: "JunctionRatio", 37: "JRPM", 38:"t_len", 39: "JPKM", 40: "delEx", 
                     41: "normalJunction", 42: "backSpliceInfo", 43: "sampleName"})

        final.to_csv(f"{output_dir}/{id1}_Shannon_score.csv", index=False)

    # print("processing completed ...", id1)
