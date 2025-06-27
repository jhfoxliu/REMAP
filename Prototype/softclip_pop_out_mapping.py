from Bio import SeqIO
import pandas as pd
import argparse
import pysam
import pandas as pd
from Bio.Seq import reverse_complement

output = {"A": {}, "G": {}, "C":{}, "T": {}}
fastq_out = {}

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

def match_fastq(fastq_in, fastq_out, dict_used, read1=False):
    fastq_out = fastq_out.split("/")[-1]
    with open(fastq_in, 'r') as input, open(fastq_out, "w") as output:
        line = input.readline()
        i = 0
        while (line):
            if i == 0:
                line = line.strip().split(" ")
                name = line[0][1:]
                if name in dict_used:
                    softclip_num, softclip_base = dict_used[name]
                    output.write("@{}_{}\n".format(name, softclip_base))
                line = input.readline()
                i += 1
            elif i == 1:
                if name in dict_used:
                    if read1 == False:
                        output.write("{}".format(line))
                    else:
                        softclip_num, softclip_base = dict_used[name]
                        new_seq = line[0: softclip_num] +  line[softclip_num+1: ]
                        output.write("{}".format(new_seq))
                line = input.readline()
                i += 1
            elif i == 2:
                if name in dict_used:
                    output.write("{}".format(line))
                line = input.readline()
                i +=1
            elif i == 3:
                if name in dict_used:
                    if read1 == False:
                        output.write("{}".format(line))
                    else:
                        softclip_num, softclip_base = dict_used[name]
                        new_qual = line[0: softclip_num] +  line[softclip_num+1:]
                        output.write("{}".format(new_qual))
                line = input.readline()
                i = 0
            

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i","--input_bed",dest="input_bed",required=True,help="input BED file")
    group_required.add_argument("-b","--input_bam",dest="input_bam",required=True,help="input BAM file")
    group_required.add_argument("-o","--output_file",dest="output_file",required=True,help="output file")
    group_required.add_argument("-1","--read1",dest="read1",required=True, help="input fastq read1")
    group_required.add_argument("-2","--read2",dest="read2",required=False, help="input fastq read2")
    group_required.add_argument("-c","--coverage",dest="coverage",required=False, type=int, default=10,help="input read count")
    group_required.add_argument("-r","--ref",dest="reference",required=True,help="reference fasta")
    options = parser.parse_args()

    BAM = options.input_bam
    bed_file = options.input_bed
    reference_genome = {}
    for seq in SeqIO.parse(options.reference,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    with open(bed_file, "r") as input_sites, pysam.AlignmentFile(BAM, "rb") as input_BAM:
        chr_names = {}
        for SQ in input_BAM.header["SQ"]:
            chr_names[SQ["SN"]] = 1

        for line in input_sites.readlines():
            line = line.strip().split("\t")
            chr, pos_0, pos_1, _, read_count, strand = line

            if int(read_count) < options.coverage:
                continue

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
                if strand == "+" and read.is_read1 == True and read.is_reverse == False:
                    if  read.reference_start != pos_0:
                        continue
                    if read.cigartuples[0][0] == 4:
                        pop_base = read.query_sequence[read.cigartuples[0][1]]
                        fastq_out[read.query_name] = (read.cigartuples[0][1], pop_base)
                    else:
                        add_to_dict(read, pos_0, KEY, output, strand)
                    
                elif strand == "-" and read.is_read1 == True and read.is_reverse == True:
                    if read.reference_end != pos_1:
                        continue
                    if read.cigartuples[-1][0] == 4:
                        pop_base = read.query_sequence[len(read.query_sequence) - 1 - read.cigartuples[-1][1]]
                        pop_base = reverse_complement(pop_base)
                        fastq_out[read.query_name] = (read.cigartuples[-1][1], pop_base)
                    else:
                        add_to_dict(read, pos_0, KEY, output, strand)

df = pd.DataFrame(output)
df.to_csv(options.output_file)
df2 = df.melt(value_name="Count", var_name="Base_type", ignore_index=False)
df2["Insert"] = 0
df2 = df2[df2["Count"]>0]
df2.to_csv(options.output_file + ".insert0.csv")

match_fastq(options.read1, options.read1.replace(".fastq", ".pop_1.fastq"), fastq_out, read1=True)
if options.read2:
    match_fastq(options.read2, options.read2.replace(".fastq", ".pop_1.fastq"), fastq_out)
