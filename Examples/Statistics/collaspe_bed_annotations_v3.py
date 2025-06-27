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

def select_annotation(annot):
    # If only one gene found, use the closest one.
    # We consider absolute distance first, if we found two annotations within 10 bp, then consider priority.
    if len(annot) == 1:
        return annot
    else:
        dict_temp = {"enst": {}, "gene": {}, "prior": {}, "bedtools_dist": {}, "real_dist": {}, "N": {}, "biotype": {}}
        N = 0
        for chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist in annot:
            dict_temp["bedtools_dist"][N] = bedtools_dist
            dict_temp["real_dist"][N] = real_dist
            dict_temp["enst"][N] = enst
            dict_temp["gene"][N] = gene
            dict_temp["prior"][N] = pri
            dict_temp["biotype"][N] = biotype1
            dict_temp["N"][N] = N
            N += 1
        df_temp = pd.DataFrame(dict_temp)
        df_temp["abs_dist"] = np.abs(df_temp["real_dist"])
        df_temp["abs_bedtools_dist"] = np.abs(df_temp["bedtools_dist"])
        df_temp = df_temp.sort_values(by=["abs_bedtools_dist", "prior", "abs_dist"])
        # print(df_temp)
        if len(df_temp["gene"].unique()) == 1:
            return annot[df_temp.iloc[0]["N"]]
        else:
            # lowest_prior = df_temp["prior"].min()
            # closest_bedtools = df_temp["abs_bedtools_dist"].min()
            closest_real = df_temp["abs_dist"].min()

            if len(df_temp["prior"].unique()) == 1 or len(df_temp["abs_bedtools_dist"].unique()) == 1:
                return annot[df_temp.iloc[0]["N"]]
            else:
                df_temp2 = df_temp[df_temp["abs_dist"] <= closest_real + 5]
                df_temp2 = df_temp2.sort_values(by=["prior", "abs_dist"])
                return annot[df_temp2.iloc[0]["N"]]



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #Required
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i" ,dest="input",required=True,help="Input BED file")
    group_required.add_argument("-a", dest="anntation",required=True,help="Annotation")
    group_required.add_argument("-o", dest="output",required=True,help="Output name")
    group_required.add_argument("--thresh-min", dest="thresh_min",required=False, type=int, default=-100,help="Distance limit, upstream to known TSS")
    group_required.add_argument("--thresh-max", dest="thresh_max",required=False, type=int, default=100,help="Distance limit, upstream to known TSS")
    options = parser.parse_args()

    dict_biotypes = {}
    with open(options.anntation, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t")
            enst = line[1]
            GENE_BIOTYPE, TRANS_BIOTYPE = line[0].split(",")

            dict_biotypes[enst] = (GENE_BIOTYPE, TRANS_BIOTYPE)

    with open(options.input, "r") as input, open(options.output, "w") as output:
        last_key = None
        temp = []
        line = input.readline()
        while(line):
            chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist = handle_line(line)
            if real_dist < options.thresh_min or real_dist > options.thresh_max:
                gene = "*" + gene + "_" + str(real_dist)
            key =(chr, pos_0, pos_1, count, strand)
            if last_key is None:
                last_key = key
                temp.append((chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist))
                line = input.readline()
            else:
                if key == last_key:
                    temp.append((chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist))
                    line = input.readline()
                elif key != last_key:
                    if len(temp) == 1:
                        info = "{}|{}|{}|TSS|{}".format(temp[0][5], temp[0][7], temp[0][8], temp[0][-1])
                        print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)

                        # clear buffer
                        temp = []
                        last_key = key
                        temp.append((chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist))
                        line = input.readline()
                    else:  # non-unique record
                        best_one = select_annotation(temp)
                        info = "{}|{}|{}|TSS|{}".format(best_one[5], best_one[7], best_one[8], best_one[-1])
                        print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)

                        # clear buffer
                        temp = []
                        last_key = key
                        temp.append((chr, pos_0, pos_1, count, strand, gene, enst,  biotype1, biotype2, pri, bedtools_dist, real_dist))
                        line = input.readline()
       
        # the last one
        key =(chr, pos_0, pos_1, count, strand)
        if key == last_key:  # non-unique
            temp.append((chr, pos_0, pos_1, count, strand, gene, enst,  biotype1, biotype2, pri, bedtools_dist, real_dist))
            best_one = select_annotation(temp)
            info = "{}|{}|{}|TSS|{}".format(best_one[5], best_one[7], best_one[8], best_one[-1])
            print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)
        elif key != last_key:
            if len(temp) == 1:
                info = "{}|{}|{}|TSS|{}".format(best_one[5], best_one[7], best_one[8], best_one[-1])
                print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)

                temp = []
                last_key = key
                temp.append((chr, pos_0, pos_1, count, strand, gene, enst, biotype1, biotype2, pri, bedtools_dist, real_dist))
            else:
                best_one = select_annotation(temp)
                info = "{}|{}|{}|TSS|{}".format(best_one[5], best_one[7], best_one[8], best_one[-1])
                print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)
                
                temp = []
                last_key = key
                temp.append((chr, pos_0, pos_1, count, strand, gene, enst,  biotype1, biotype2, pri, bedtools_dist, real_dist))
            
            info = "{}|{}|{}|TSS|{}".format(temp[0][5], temp[0][7], temp[0][8], temp[0][-1])
            print("{}\t{}\t{}\t{}\t{}\t{}".format(last_key[0], last_key[1], last_key[2], info, last_key[3], last_key[4]), file=output)
