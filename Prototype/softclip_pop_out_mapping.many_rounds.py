from Bio import SeqIO
import pandas as pd
import argparse, os, sys
import pysam
import pandas as pd
from Bio.Seq import reverse_complement
import re

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

def get_report(BED, BAM, read1, read2=None, out_prefix=None, next_read1=None, next_read2=None):
    global current_round
    output = {"A": {}, "G": {}, "C":{}, "T": {}}
    output_insert = {}
    fastq_out = {}
    with open(BED, "r") as input_sites, pysam.AlignmentFile(BAM, "rb") as input_BAM:
        chr_names = {}
        for SQ in input_BAM.header["SQ"]:
            chr_names[SQ["SN"]] = 1

        for line in input_sites.readlines():
            line = line.strip().split("\t")
            chr, pos_0, pos_1, _, read_count, strand = line

            # disable read count threshold
            
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
                if strand == "+" and read.is_read1 == True and read.is_reverse == False:
                    if  read.reference_start != pos_0:
                        continue
                    if read.cigartuples[0][0] == 4:
                        pop_base = read.query_sequence[read.cigartuples[0][1]]
                        fastq_out[read.query_name] = (read.cigartuples[0][1], pop_base)
                    else:
                        insert_type = read.query_name.split("_")[-1]
                        if insert_type not in output_insert:
                            output_insert[insert_type] = {}
                        if KEY not in output_insert[insert_type]:
                            output_insert[insert_type][KEY] = 0
                        output_insert[insert_type][KEY] += 1
                        add_to_dict(read, pos_0, KEY, output, strand)
                    
                elif strand == "-" and read.is_read1 == True and read.is_reverse == True:
                    if read.reference_end != pos_1:
                        continue
                    if read.cigartuples[-1][0] == 4:
                        pop_base = read.query_sequence[len(read.query_sequence) - 1 - read.cigartuples[-1][1]]
                        pop_base = reverse_complement(pop_base)
                        fastq_out[read.query_name] = (read.cigartuples[-1][1], pop_base)
                    else:
                        insert_type = read.query_name.split("_")[-1]
                        if insert_type not in output_insert:
                            output_insert[insert_type] = {}
                        if KEY not in output_insert[insert_type]:
                            output_insert[insert_type][KEY] = 0
                        output_insert[insert_type][KEY] += 1
                        add_to_dict(read, pos_0, KEY, output, strand)

    df = pd.DataFrame(output)
    df.to_csv(out_prefix + ".csv")
    df = pd.DataFrame(output_insert)
    
    df2 = df.melt(value_name="Count", var_name="Base_type", ignore_index=False)
    df2["Insert"] = current_round
    df2 = df2[df2["Count"]>0]
    # df2.to_csv(options.output_file + ".insert.csv")
    df2.to_csv(out_prefix + ".insert.csv")

    next_read1 = re.sub("pop_.*.fastq", "pop_{current_round}.fastq".format(current_round=current_round+1), read1)
    next_read1 = next_read1.split("/")[-1]

    if read2:
        next_read2 = re.sub("pop_.*.fastq", "pop_{current_round}.fastq".format(current_round=current_round+1),read2)
        next_read2 = next_read2.split("/")[-1]
    else:
        next_read2 = None

    match_fastq(read1, next_read1, fastq_out, read1=True)
    if read2:
        match_fastq(read2, next_read2, fastq_out, read1=False)
    
    return next_read1, next_read2

def match_fastq(fastq_in, fastq_out, dict_used, read1=False):
    with open(fastq_in, 'r') as input, open(fastq_out, "w") as output:
        line = input.readline()
        i = 0
        # print(dict_used)
        while (line):
            if i == 0:
                line = line.strip().split(" ")
                name = line[0][1:]
                # print(name)
                if name in dict_used:
                    softclip_num, softclip_base = dict_used[name]
                    output.write("@{}{}\n".format(name, softclip_base))
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

def mapping(read1, read2=None, out_prefix=None, p=4):
    if read1 and read2:
        os.system("{hisat2} --no-temp-splicesite --no-unal --fr --rna-strandness FR --no-discordant --no-mixed -p {p} -x {index} -1 {read1} -2 {read2} -S {out_sam} ".format(
            hisat2 = options.hisat2,
            read1=read1, read2=read2, 
            index=options.index, 
            out_sam = out_prefix + ".sam", 
            p=p
        ))

        os.system("samtools view -q 60 -bS {out_sam} > {out_bam}".format(
            out_sam = out_prefix + ".sam", 
            out_bam = out_prefix + ".bam", 
        ))
        os.system("samtools sort -@ 4 -m 2G -o {out_sorted_bam} {out_bam}".format(
            out_bam = out_prefix + ".bam",
            out_sorted_bam = out_prefix + ".sorted.bam",
        ))
        os.system("samtools index {out_sorted_bam}".format(
            out_sorted_bam = out_prefix + ".sorted.bam",
        ))
        os.system("rm {out_sam} {out_bam} ".format(
            out_sam = out_prefix + ".sam", 
            out_bam = out_prefix + ".bam",
        ))
        os.system("/home/fox/Software/bin/python /home/fox/Scripts/ReCappable-seq/Align_to_cap_starts_PE_fix_softclips.all.py {out_sorted_bam}  > {bed_out}".format(
            out_sorted_bam = out_prefix + ".sorted.bam",
            bed_out = out_prefix + ".bed"
        ))

if __name__ == "__main__":
    description = """"""
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
    #Require
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-1","--read1",dest="read1",required=True, help="input fastq read1")
    group_required.add_argument("-2","--read2",dest="read2",required=False, help="input fastq read2")
    group_required.add_argument("--hisat2", dest="hisat2",required=False, default="/home/fox/Software/hisat2/2.2.1/hisat2", help="hisat2")
    group_required.add_argument("-x","--hisat2_index",dest="index",required=False, default="/home/fox/Database/hisat2/GRCh38_r104/Homo_sapiens.GRCh38.dna_sm.primary_assembly.format", help="hisat2 index")
    group_required.add_argument("-r","--reference",dest="reference_genome",required=False, default="/home/fox/Database/Genome/Human/Ensembl/GRCh38_r104/Homo_sapiens.GRCh38.dna_sm.primary_assembly.format.fa", help="reference")
    group_required.add_argument("-n", dest="trim_times",required=False, type=int, default=10, help="Trimming trails")
    options = parser.parse_args()

    reference_genome = {}
    for seq in SeqIO.parse(options.reference_genome,"fasta"):
        reference_genome[seq.id] = str(seq.seq).upper()

    N = 0
    # primer round
    print("Softclip 1 base:")

    current_round = 1
    read1 = options.read1
    read2 = options.read2
    mapping(read1, read2=read2, out_prefix="Softclip_{current_round}".format(current_round=current_round))

    read1, read2 = get_report("Softclip_{current_round}.bed".format(current_round=current_round), 
                                "Softclip_{current_round}.sorted.bam".format(current_round=current_round), 
                                read1, 
                                read2=read2, 
                                out_prefix="Softclip_{current_round}".format(current_round=current_round),
                                )

    N += 1

    print("=" * 20)
    print(" " * 20)
    
    # iterate until fastq file empty
    
    while (N < options.trim_times):
        current_round = current_round + 1
        read1_size = os.path.getsize("./{}".format(read1))
        N += 1
        if read1_size > 0:
            print("Softclip {} base:".format(current_round))
            mapping(read1, 
                    read2=read2, 
                    out_prefix="Softclip_{current_round}".format(current_round=current_round)
                    )

            read1, read2 =  get_report("Softclip_{current_round}.bed".format(current_round=current_round), 
                                        "Softclip_{current_round}.sorted.bam".format(current_round=current_round), 
                                        read1, 
                                        read2=read2, 
                                        out_prefix="Softclip_{current_round}".format(current_round=current_round),
                                        )

            print("=" * 20)
            print(" " * 20)
        else:
            break
