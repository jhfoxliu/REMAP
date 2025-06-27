import pandas as pd
import sys, os, pysam, argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement

output = {"A": {}, "G": {}, "C":{}, "T": {}}

def add_to_dict(read, pos_0, KEY, target_dict, strand):
    if strand == "+":
        for query_pos, ref_pos in read.get_aligned_pairs(with_seq=False):
            if query_pos is not None and ref_pos is not None:
                if ref_pos == pos_0:
                    query_base = read.query_sequence[query_pos]
                    if query_base != "N":
                        target_dict[query_base][KEY] += 1
                    break
    elif strand == "-":
        for query_pos, ref_pos in read.get_aligned_pairs(with_seq=False):
            if query_pos is not None and ref_pos is not None:
                if ref_pos == pos_0:
                    query_base = read.query_sequence[query_pos]
                    query_base = reverse_complement(query_base)
                    if query_base != "N":
                        target_dict[query_base][KEY] += 1
                    break

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    # group_required.add_argument("-i","--input_bed",dest="input_bed",required=True,help="input BED file")
    group_required.add_argument("-b","--input_bam",dest="input_bam",required=True,help="input BAM file")
    group_required.add_argument("-o","--output_file",dest="output_file",required=True,help="output file")
    group_required.add_argument("-c","--coverage",dest="coverage",required=False, type=int, default=10,help="input read count")
    group_required.add_argument("-r","--ref",dest="reference",required=True,help="reference fasta")
    options = parser.parse_args()

    dfs = []

    for fn in os.listdir("./"):
        if fn.endswith(".insert.csv") or fn.endswith(".insert0.csv"):
            try:
                subdf = pd.read_csv(fn, index_col=[0,1,2,3,4], header=0)
            except IndexError:
                continue
            dfs.append(subdf)

    dict_sites = {}
    df_sites = pd.concat(dfs)
    with open("temp.bed", "w") as fn_out:
        for idx in df_sites.index:
            dict_sites[(idx[0], idx[1], idx[2])] = 1
        for keys in dict_sites.keys():
            fn_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(keys[0], keys[1]-1, keys[1], 999, 999, keys[2]))

    BAM = options.input_bam
    reference_genome = {}
    for seq in SeqIO.parse(options.reference,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    with open("temp.bed", "r") as input_sites, pysam.AlignmentFile(BAM, "rb") as input_BAM:
        chr_names = {}
        for SQ in input_BAM.header["SQ"]:
            chr_names[SQ["SN"]] = 1

        for line in input_sites.readlines():
            line = line.strip().split("\t")
            chr, pos_0, pos_1, _, read_count, strand = line

            # if int(read_count) < options.coverage:
            #     continue

            pos_0 = int(pos_0)
            pos_1 = int(pos_1)
            ref_base = reference_genome[chr][pos_0]
            if strand == "-":
                ref_base = reverse_complement(ref_base)
            KEY = (chr, pos_1, strand, ref_base)

            output["A"][KEY] = 0
            output["G"][KEY] = 0
            output["C"][KEY] = 0
            output["T"][KEY] = 0

            chr_v2 = chr
            if "M" in chr_names:
                if chr_v2 == "MT":
                    chr_v2 = "M"
            elif "MT" in chr_names:
                if chr_v2 == "M":
                    chr_v2 = "MT"

            for read in input_BAM.fetch(chr_v2, pos_0, pos_1, multiple_iterators=True):
                if strand == "+" and read.is_reverse == False:
                    if  read.reference_start != pos_0:
                        continue
                    if read.cigartuples[0][0] == 4:
                        continue
                    else:
                        add_to_dict(read, pos_0, KEY, output, strand)
                    
                elif strand == "-" and read.is_reverse == True:
                    if read.reference_end != pos_1:
                        continue
                    if read.cigartuples[-1][0] == 4:
                        continue
                    else:
                        add_to_dict(read, pos_0, KEY, output, strand)

    df = pd.DataFrame(output)
    df2 = df.melt(value_name="Count", var_name="Base_type", ignore_index=False)
    df2["Insert"] = 0
    df2 = df2[df2["Count"]>0]
    df2.to_csv("SoftClip0.refetch.insert.csv")

    dfs_2 = []

    for fn in os.listdir("./"):
        if fn.endswith(".insert.csv"):
            try:
                subdf = pd.read_csv(fn, index_col=[0,1,2,3,4], header=0)
            except IndexError:
                continue
            dfs_2.append(subdf)

    df = pd.concat(dfs_2).reset_index()
    df.columns=["chr", "pos", "strand", "ref_base", "softclip_bases", "Count", "Softclip_size"]
    df["Count"] = df["Count"].astype(int)
    df = df.sort_values(by=["chr", "pos", "strand", "Softclip_size"])
    # df = df.reindex()
    df.to_csv(options.output_file + ".out.csv", index=False)

    df["Count sum"] =  df.groupby(by=["chr", "pos", "strand","ref_base"])[["Count"]].transform(lambda x: x.sum())
    df2 = df[df["Count sum"]>options.coverage]
    df2.to_csv(options.output_file + ".filtered_{}.csv".format(options.coverage), index=False)