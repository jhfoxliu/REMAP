import sys
import argparse
import pandas as pd
import numpy as np

priorities = {"snRNA":0, "snoRNA":1, "protein_coding":2, "lncRNA": 3}

def handle_line(line):
    line = line.strip().split("\t")
    chr, pos_0, pos_1, _, count, strand, _, tx_start, tx_end, gene, enst, _, bedtools_dist = line
    BIO = dict_biotypes.get(enst)
    if BIO is None:
        bio1 = "NONE"
        bio2 = "NONE"
    else:
        bio1, bio2 = BIO
    pri = priorities.get(bio2)
    if pri is None:
        pri = 9999
    if strand == "+":
        real_dist = int(pos_0) - int(tx_start)
    else:
        real_dist = int(tx_end) - 1 - int(pos_0)
    bedtools_dist = int(bedtools_dist)
    return chr, pos_0, pos_1, count, strand, gene, enst, bio1, bio2, pri, bedtools_dist, real_dist # , real_dist

def process_line(row):
    global output
    chr, pos_0, pos_1, info, count, strand = row
    key = (chr, pos_0, pos_1, strand)
    all_exons = dict_all_exons.get(key)
    if all_exons is None:
        print("\t".join(row), file=output)
    else:
        dict_temp = {"enst": {}, "gene": {}, "prior": {}, "bedtools_dist": {}, "real_dist": {}, "N": {}, "biotype": {}, "biotype2": {}}
        N = 0
        for gene, enst, bio1, bio2, pri, bedtools_dist, real_dist in all_exons:
            dict_temp["bedtools_dist"][N] = int(bedtools_dist)
            dict_temp["real_dist"][N] = int(real_dist)
            dict_temp["enst"][N] = enst
            dict_temp["gene"][N] = gene
            dict_temp["prior"][N] = pri
            dict_temp["biotype"][N] = bio1
            dict_temp["biotype2"][N] = bio2
            dict_temp["N"][N] = N
            N += 1
        df_temp = pd.DataFrame(dict_temp)
        df_temp["abs_dist"] = np.abs(df_temp["real_dist"])
        df_temp["abs_bedtools_dist"] = np.abs(df_temp["bedtools_dist"])
        df_temp = df_temp.sort_values(by=["abs_bedtools_dist", "prior", "abs_dist"])
        # print(df_temp)
        if len(df_temp["gene"].unique()) == 1:
            GENE = df_temp.iloc[0]["gene"]
            biotype = df_temp.iloc[0]["biotype"]
            biotype2 = df_temp.iloc[0]["biotype2"]
            if df_temp.iloc[0]["bedtools_dist"] == 0:
                dist = 0
            else:
                dist = df_temp.iloc[0]["real_dist"]
        else:
            if df_temp["bedtools_dist"].min() == 0:
                df_temp = df_temp[df_temp["bedtools_dist"] == 0]
                if len(df_temp["gene"].unique()) == 1:
                    GENE = df_temp.iloc[0]["gene"]
                    biotype = df_temp.iloc[0]["biotype"]
                    biotype2 = df_temp.iloc[0]["biotype2"]
                    dist = 0
                else:
                    # lowest_prior = df_temp["prior"].min()
                    # closest_bedtools = df_temp["abs_bedtools_dist"].min()
                    closest_real = df_temp["abs_dist"].min()
                    if len(df_temp["prior"].unique()) == 1 or len(df_temp["abs_bedtools_dist"].unique()) == 1:
                        GENE = df_temp.iloc[0]["gene"]
                        biotype = df_temp.iloc[0]["biotype"]
                        biotype2 = df_temp.iloc[0]["biotype2"]
                        dist = 0
                    else:
                        df_temp2 = df_temp[df_temp["abs_dist"] <= closest_real + 10]
                        df_temp2 = df_temp2.sort_values(by=["prior", "abs_dist"])
                        GENE = df_temp.iloc[0]["gene"]
                        biotype = df_temp.iloc[0]["biotype"]
                        biotype2 = df_temp.iloc[0]["biotype2"]
                        dist = 0
            else:
                closest_real = df_temp["abs_dist"].min()
                if len(df_temp["prior"].unique()) == 1 or len(df_temp["abs_bedtools_dist"].unique()) == 1:
                    GENE = df_temp.iloc[0]["gene"]
                    biotype = df_temp.iloc[0]["biotype"]
                    biotype2 = df_temp.iloc[0]["biotype2"]
                    dist = df_temp.iloc[0]["real_dist"]
                else:
                    df_temp2 = df_temp[df_temp["abs_dist"] <= closest_real + 10]
                    df_temp2 = df_temp2.sort_values(by=["prior", "abs_dist"])
                    GENE = df_temp.iloc[0]["gene"]
                    biotype = df_temp.iloc[0]["biotype"]
                    biotype2 = df_temp.iloc[0]["biotype2"]
                    dist = df_temp.iloc[0]["real_dist"]
        if dist > options.thresh_max or dist < options.thresh_min:
            GENE = "*{}_{}".format(GENE, dist)
            INFO = "UNKNOWN|UNKNOWN|UNKNOWN|UNKNOWN|NA" # .format(GENE, biotype, biotype2, dist)
        else:
            INFO = "{}|{}|{}|UPSTREAM|{}".format(GENE, biotype, biotype2, dist)
        print("{}\t{}\t{}\t{}\t{}\t{}".format(chr, pos_0, pos_1, INFO, count, strand), file=output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #Required
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i" ,dest="input",required=True,help="Input BED file, collapsed")
    group_required.add_argument("-I" ,dest="input2",required=True,help="Input BED file, closest")
    group_required.add_argument("-a", dest="anntation",required=True,help="Annotation")
    group_required.add_argument("-o", dest="output",required=True,help="Output name")
    group_required.add_argument("--thresh-min", dest="thresh_min",required=False, type=int, default=-2000,help="Distance limit, downstream")
    group_required.add_argument("--thresh-max", dest="thresh_max",required=False, type=int, default=1,help="Distance limit, upstream")
    options = parser.parse_args()

    dict_biotypes = {}
    with open(options.anntation, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t")
            enst = line[1]
            GENE_BIOTYPE, TRANS_BIOTYPE = line[0].split(",")

            dict_biotypes[enst] = (GENE_BIOTYPE, TRANS_BIOTYPE)

    dict_all_exons = {}
    with open(options.input2, "r") as input:
        for line in input.readlines():
            line = line.strip().split()
            chr, pos_0, pos_1, _, count, strand, _, tx_start, tx_end, gene, enst, _, bedtools_dist = line
            key = (chr, pos_0, pos_1, strand)
            if key not in dict_all_exons:
                dict_all_exons[key] = []
            BIO = dict_biotypes.get(enst)
            if BIO is None:
                bio1 = "NONE"
                bio2 = "NONE"
            else:
                bio1, bio2 = BIO
            pri = priorities.get(bio2)
            if pri is None:
                pri = 9999
            if strand == "+":
                real_dist = int(pos_0) - int(tx_start)
            else:
                real_dist = int(tx_end) - 1 - int(pos_0)
            dict_all_exons[key].append([gene, enst, bio1, bio2, pri, bedtools_dist, real_dist])

    with open(options.input, "r") as input, open(options.output, "w") as output:
        for line in input.readlines():
            row = line.strip().split("\t")
            if row[3].startswith("*") == False:
                print("\t".join(row), file=output)
            else:
                process_line(row)

