import sys, os
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import reverse_complement, complement
from itertools import product
import argparse
import re
import multiprocessing
from multiprocessing import Process,Pool
import signal
import time
from time import gmtime, strftime

dict_homopolymer = {}
for i in range(2, 20):
    for j in ["A", "T", "C", "G"]:
        dict_homopolymer[j*i] = 1

def handle_softclipping(read, forced_softclip=None):
    if forced_softclip is None:
        if read.is_reverse == False:
            start_pos = read.reference_start
            genomic_N3 = ref_genome[read.reference_name][start_pos: start_pos+3]
            genomic_upstream_one = ref_genome[read.reference_name][start_pos-1]
            genomic_upstream_two = ref_genome[read.reference_name][start_pos-2]
            query_sequence = read.query_sequence
            N_softclip = read.cigartuples[0][1]
            seq_softclip = query_sequence[0: N_softclip]
        else:
            start_pos = read.reference_end - 1 
            genomic_N3 = reverse_complement(ref_genome[read.reference_name][start_pos-3: start_pos+1])
            genomic_upstream_one = reverse_complement(ref_genome[read.reference_name][start_pos])
            genomic_upstream_two = reverse_complement(ref_genome[read.reference_name][start_pos+1])
            query_sequence = reverse_complement(read.query_sequence)
            N_softclip = read.cigartuples[-1][1]
            seq_softclip = query_sequence[0: N_softclip]
    else:
        if read.is_reverse == False:
            start_pos = read.reference_start + len(forced_softclip)
            genomic_N3 = ref_genome[read.reference_name][start_pos: start_pos+3]
            genomic_upstream_one = ref_genome[read.reference_name][start_pos-1]
            genomic_upstream_two = ref_genome[read.reference_name][start_pos-2]
            query_sequence = read.query_sequence
            N_softclip = len(forced_softclip)
            seq_softclip = forced_softclip
            # print(forced_softclip, N_softclip, seq_softclip, read.reference_start, start_pos)
        else:
            start_pos = read.reference_end - len(forced_softclip) - 1
            genomic_N3 = reverse_complement(ref_genome[read.reference_name][start_pos-3: start_pos+1])
            genomic_upstream_one = reverse_complement(ref_genome[read.reference_name][start_pos])
            genomic_upstream_two = reverse_complement(ref_genome[read.reference_name][start_pos+1])
            query_sequence = reverse_complement(read.query_sequence)
            N_softclip = len(forced_softclip)
            seq_softclip = forced_softclip

    # let's begin with the simplest situations
    if N_softclip == 1:
        # one base expansion
        if seq_softclip == genomic_N3[0]:
            return 1, seq_softclip, "H", start_pos  # SL, SS, ST, GP
        if seq_softclip == genomic_upstream_one:
            if read.is_reverse == False:
                return 1, query_sequence[1], "H1", start_pos-1  # SL, SS, ST, GP
            else:
                return 1, query_sequence[1], "H1", start_pos+1  # SL, SS, ST, GP
        if seq_softclip == genomic_upstream_two:
            if read.is_reverse == False:
                return -1, "D", "SK", start_pos-2  # SL, SS, ST, GP
            else:
                return -1, "D", "SK", start_pos+2  # SL, SS, ST, GP
        # otherwise 
        # test if there is a SNP:
        if options.snp:
            with pysam.VariantFile(options.snp, mode="rb", index_filename=options.snp+".csi") as BCF:
                if read.is_reverse == False:
                    for records in BCF.fetch(contig=chr, start=start_pos-1, end=start_pos):
                        return 0, seq_softclip , "S", start_pos-1  # SL, SS, ST, GP
                else:
                    for records in BCF.fetch(contig=chr, start=start_pos+1, end=start_pos+2):
                        return 0, seq_softclip , "S", start_pos+1  # SL, SS, ST, GP
        
        # Then, another possiblity is reanneal, but leave one nt overhang
        seq_softclip_first_three = seq_softclip[0: 3]
        possible_repeats = [a.start() for a in re.finditer(seq_softclip_first_three, query_sequence[:20])]
        if len(possible_repeats) > 1:
            N_softclip = possible_repeats[1]
            seq_softclip = query_sequence[0: N_softclip]
            # guess the repeated motif
            five_prime = query_sequence[0: 3]
            gapped = query_sequence[3: N_softclip]
            for i in range(3, N_softclip):
                five_prime = query_sequence[0: i]
                TSS_motif = query_sequence[N_softclip: N_softclip + i]
                if five_prime == TSS_motif:
                    gapped = query_sequence[i: N_softclip]
                else:
                    five_prime = query_sequence[0: i-1]
                    gapped = query_sequence[i-1: N_softclip]
                    break
            
            if read.is_reverse == False:
                real_mapped_pos = start_pos + possible_repeats[1] - 1
            else:
                real_mapped_pos = start_pos - possible_repeats[1] + 1

            return 1, seq_softclip , "RA_{}_{}".format(five_prime, gapped), real_mapped_pos

        return 1, seq_softclip , "UN", start_pos  # SL, SS, ST, GP

    # the softclipping bases are identical to the homopolymer
    if seq_softclip in dict_homopolymer:
        if seq_softclip[0:2] == genomic_N3[0:2]:
            return N_softclip, seq_softclip , "H", start_pos  # SL, SS, ST, GP

    # then let's remove the first base
    seq_softclip_first = seq_softclip[0]
    seq_softclip_remove_first = seq_softclip[1:]
    if read.is_reverse == False:
        genomic_upstream_one_base = ref_genome[read.reference_name][start_pos-1]
    else:
        genomic_upstream_one_base = reverse_complement(ref_genome[read.reference_name][start_pos+1])
    if (len(seq_softclip_remove_first) == 1 and seq_softclip_remove_first == genomic_N3[0] and  seq_softclip_first == genomic_upstream_one_base) or (seq_softclip_remove_first in dict_homopolymer and seq_softclip_remove_first[0:2] == genomic_N3[0:2] and seq_softclip_first == genomic_upstream_one_base):
        if read.is_reverse == False:
            return N_softclip-1, seq_softclip_remove_first, "H1" , start_pos - 1 # SL, SS, ST, GP
        else:
            return N_softclip-1, seq_softclip_remove_first, "H1" , start_pos + 1 # SL, SS, ST, GP

    # then let's remove the first two base
    seq_softclip_first_two = seq_softclip[0: 2]
    seq_softclip_remove_first_two = seq_softclip[2:]
    if read.is_reverse == False:
        genomic_upstream_two_base = ref_genome[read.reference_name][start_pos-2: start_pos]
    else:
        genomic_upstream_two_base = reverse_complement(ref_genome[read.reference_name][start_pos+1: start_pos+3])

    if (len(seq_softclip_remove_first_two) == 1 and genomic_N3[0] == seq_softclip_remove_first_two and seq_softclip_first_two == genomic_upstream_two_base) or (seq_softclip_remove_first_two in dict_homopolymer and seq_softclip_remove_first_two[0:2] == genomic_N3[0:2] and seq_softclip_first_two == genomic_upstream_two_base):
        if read.is_reverse == False:
            return N_softclip-2, seq_softclip_remove_first_two, "H2", start_pos - 2 # SL, SS, ST, GP
        else:
            return N_softclip-2, seq_softclip_remove_first_two, "H2", start_pos + 2 # SL, SS, ST, GP
        
    # then let's remove the first three base
    seq_softclip_first_three = seq_softclip[0: 3]
    seq_softclip_remove_first_three = seq_softclip[3:]
    if read.is_reverse == False:
        genomic_upstream_three_base = ref_genome[read.reference_name][start_pos-3: start_pos]
    else:
        genomic_upstream_three_base = reverse_complement(ref_genome[read.reference_name][start_pos+1: start_pos+4])
    if (len(seq_softclip_remove_first_three) == 1 and genomic_N3[0] == seq_softclip_remove_first_three and seq_softclip_first_three == genomic_upstream_three_base) or (seq_softclip_remove_first_three in dict_homopolymer and seq_softclip_remove_first_three[0:2] == genomic_N3[0:2] and seq_softclip_first_three == genomic_upstream_three_base):
        if read.is_reverse == False:
            return N_softclip-2, seq_softclip_remove_first_three, "H3", start_pos - 3 # SL, SS, ST, GP
        else:
            return N_softclip-2, seq_softclip_remove_first_three, "H3", start_pos + 3 # SL, SS, ST, GP

    # now let's handle some other situations for reannealing
    # note: this one only consider that the next repeat is well mapped!
    seq_softclip_first_three = seq_softclip[0: 3]
    possible_repeats = [a.start() for a in re.finditer(seq_softclip_first_three, query_sequence[:20])]
    if len(possible_repeats) > 1:
        N_softclip = possible_repeats[1]
        seq_softclip = query_sequence[0: N_softclip]
        # guess the repeated motif
        five_prime = query_sequence[0: 3]
        gapped = query_sequence[3: N_softclip]
        for i in range(3, N_softclip):
            five_prime = query_sequence[0: i]
            TSS_motif = query_sequence[N_softclip: N_softclip + i]
            if five_prime == TSS_motif:
                gapped = query_sequence[i: N_softclip]
            else:
                five_prime = query_sequence[0: i-1]
                gapped = query_sequence[i-1: N_softclip]
                break

        return N_softclip, seq_softclip , "RA_{}_{}".format(five_prime, gapped), start_pos

    return N_softclip, seq_softclip, "UN", start_pos
    

def check_5prime_mismatch(read, max_len=10): # by ChatGPT
    md_tag = read.get_tag("MD")
    md_parts = re.findall(r'(\d+)|([A-Z])|\^', md_tag)  # Extract matches, mismatches, and deletions
    query_pos = 0
    last_mutation_pos = None
    if read.is_reverse == False:
        for part in md_parts:
            if part[0]:  # Number: matched bases
                query_pos += int(part[0])
            elif part[1]:  # Letter: mismatch
                if query_pos < max_len:  # Check if within the first 10 nt
                    last_mutation_pos = query_pos + 1  # Convert to 1-based position
                query_pos += 1  # Move to the next query position

            if query_pos >= max_len:  # Stop processing after 10 nucleotides
                break
    else:
        for part in md_parts[::-1]:
            if part[0]:  # Number: matched bases
                query_pos += int(part[0])
            elif part[1]:  # Letter: mismatch
                if query_pos < max_len:  # Check if within the first 10 nt
                    last_mutation_pos = query_pos + 1  # Convert to 1-based position
                query_pos += 1  # Move to the next query position

            if query_pos >= max_len:  # Stop processing after 10 nucleotides
                break
        if last_mutation_pos is not None and read.is_reverse:
            read_length = read.query_length
            last_mutation_pos = read_length - last_mutation_pos + 1
        
    return last_mutation_pos # only return the position as the read

def correct_softclipping_in_bam_file(chr, bam_out):
    with pysam.AlignmentFile(options.bam, "rb") as BAM_IN, \
         pysam.AlignmentFile(bam_out, "wb", template=BAM_IN) as BAM_OUT:
        #  pysam.AlignmentFile(bam_out+".mismatch.bam", "wb", template=BAM_IN) as BAM_OUT_mut :
        for read in BAM_IN.fetch(contig=chr):
            if options.single_end == False and read.is_read2 == True:
                continue
            if read.is_secondary == True:
                continue

            if read.is_reverse == False and read.cigartuples[0][0] != 4: # perfect 5' end
                last_five_prime_mistach = check_5prime_mismatch(read)
                if last_five_prime_mistach is None:
                    read.set_tag("SL", 0, value_type="i")  # Soft Clipping Length
                    read.set_tag("SS", "", value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("ST", "N", value_type="Z")  # Soft Clipping Type
                    read.set_tag("GP", read.reference_start, value_type="i") # Genomic position, 0-based
                    BAM_OUT.write(read)
                else:
                    forced_softclip = read.query_sequence[0: last_five_prime_mistach]
                    SL, SS, ST, GP = handle_softclipping(read, forced_softclip=forced_softclip)
                    read.set_tag("SL", SL, value_type="i")  # Soft Clipping Length
                    read.set_tag("SS", SS, value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("ST", ST, value_type="Z")  # Soft Clipping Type
                    read.set_tag("GP", GP, value_type="i") # Genomic position, 0-based
                    BAM_OUT.write(read)

            elif read.is_reverse == True and read.cigartuples[-1][0] != 4:
                last_five_prime_mistach = check_5prime_mismatch(read)
                if not last_five_prime_mistach:
                    read.set_tag("SL", 0, value_type="i")  # Soft Clipping Length
                    read.set_tag("SS", "", value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("ST", "N", value_type="Z")  # Soft Clipping Type
                    read.set_tag("GP", read.reference_end-1, value_type="i") # Genomic position, 0-based
                    BAM_OUT.write(read)
                else:
                    forced_softclip = reverse_complement(read.query_sequence[last_five_prime_mistach-1:])
                    SL, SS, ST, GP = handle_softclipping(read, forced_softclip=forced_softclip)
                    read.set_tag("SL", SL, value_type="i")  # Soft Clipping Length
                    read.set_tag("SS", SS, value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("ST", ST, value_type="Z")  # Soft Clipping Type
                    read.set_tag("GP", GP, value_type="i") # Genomic position, 0-based
                    BAM_OUT.write(read)

            else:
                '''
                ST code:
                N: Non-softclipping
                S: SNP
                H: Homopolymer immediately
                H1: Homopolymer with first base overhang
                H2: Homopolymer with first tow bases overhang
                UN: Unknown
                RA_{seq1}_{seq2}: Ranneal. {seq1} is the repeat motif, seq2 is the gapped sequence
                SK: Skip
                '''
                try:
                    SL, SS, ST, GP = handle_softclipping(read)
                    read.set_tag("SL", SL, value_type="i")  # Soft Clipping Length
                    read.set_tag("SS", SS, value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("ST", ST, value_type="Z")  # Soft Clipping Sequence
                    read.set_tag("GP", GP, value_type="i") # Genomic position, 0-based
                    BAM_OUT.write(read)
                # read.set_tag("FP", read.query_sequence[0:10], value_type="Z") # Genomic position, 0-based
                except IndexError:
                    pass
            

def signal_handler(sig,frame):
	pool.terminate()
	sys.exit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="amend_extrabases_annotations",fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
    #Required
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-b",dest="bam",required=True,help="Input BAM file")
    group_required.add_argument("-f",dest="fasta",required=True,help="Fasta file")
    group_required.add_argument("-s",dest="snp",required=False,help="SNP annotation, bcf with index")
    group_required.add_argument("--output","-o",dest="output",required=True,help="Output")
    group_required.add_argument("--bam-out",dest="bam_out",required=False, default=None,help="Output bam file (read1 only)")
    group_required.add_argument("--skip-correct",dest="skip_correct",required=False, default=False, action="store_true",help="Skip correcting BAM file, start with file of --bam-out")
    group_required.add_argument("-p",dest="processors",required=False, default=4, type=int,help="number of processors")
    group_required.add_argument("--single-end",dest="single_end",required=False, default=False, action="store_true",help="Single end")
    group_other = parser.add_argument_group("Other")
    group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
    options = parser.parse_args()
    
    print("[{t}] Loading FASTA...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    ref_genome = {}
    for seq in SeqIO.parse(options.fasta, "fasta"):
        ref_genome[seq.id] = str(seq.seq)
    print("[{t}] FASTA loaded.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    if options.bam_out is None:
        options.bam_out = options.bam.replace(".bam", ".softclip_corrected.bam")

    if options.skip_correct == False:
        BAM_reads = 0
        with pysam.AlignmentFile(options.bam, "rb") as BAM_IN:
            chr_names = {}
            for SQ in BAM_IN.header["SQ"]:
                chr_reads = BAM_IN.count(SQ["SN"])
                BAM_reads += chr_reads
                if chr_reads > 0:
                    chr_names[SQ["SN"]] = 1
        print("[{t}] Reads in BAM file (read1 + read2 + seconadary + etc.): {b}".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime()), b=BAM_reads))
  
        pid = os.getpid()


        # test
        # all_temp_bam = []
        # for chr in chr_names.keys():
        #     temp_bam = "temp_{}_{}.bam".format(chr, pid)
        #     all_temp_bam.append(temp_bam)
        #     correct_softclipping_in_bam_file(chr, temp_bam)
        # sys.exit()

        signal.signal(signal.SIGINT,signal_handler)
        # call_all_bases_from_a_read("X", BAM, main_pid)
        pool = multiprocessing.Pool(options.processors)

        print("[{t}] Analyzing BAM file with {p} processors.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime()), p=options.processors))
        try:
            all_temp_bam = []
            for chr in chr_names.keys():
                # correct_softclipping_in_bam_file(chr, temp_bam)
                temp_bam = "temp_{}_{}.bam".format(chr, pid)
                all_temp_bam.append(temp_bam)
                pool.apply_async(correct_softclipping_in_bam_file,args=(chr, temp_bam,))
            pool.close()
            pool.join()
                
            pysam.merge(options.bam_out+"_temp.bam", *all_temp_bam)
            pysam.sort("-o", options.bam_out, "-@ 4", "-m 4G", options.bam_out+"_temp.bam")
            for i in all_temp_bam:
                os.remove(i)
            os.remove(options.bam_out+"_temp.bam")
            pysam.index(options.bam_out)
            print("[{t}] BAM file with corrected 5' annotation: {b}".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime()), b=options.bam_out))

        finally:
            pool.terminate()
    
    print("[{t}] Reading corrected BAM file...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    with pysam.AlignmentFile(options.bam_out, "rb") as BAM_IN:
        data = {"Count": {}}
        for read in BAM_IN:
            chr = read.reference_name
            if read.is_reverse == False:
                strand = "+"
            else:
                strand = "-"
            pos = read.get_tag("GP") # 0-based
            ref_base = ref_genome[chr][pos]
            if strand == "-":
                ref_base = reverse_complement(ref_base)
            softclip_size = read.get_tag("SL")
            softclip_bases = read.get_tag("SS")
            event =  read.get_tag("ST")
            # FP =  read.get_tag("FP")
            if event == "N":
                COMMENT = "M"
                softclip_bases = ref_base
            elif event == "S":
                COMMENT = "+1 SNP"
            elif event == "H":
                COMMENT = "EXPANDED"
            elif event == "H1":
                COMMENT = "+2 EXPANDED"
            elif event == "H2":
                COMMENT = "+3 EXPANDED"
            elif event == "SK":
                COMMENT = "SKIPPED"
            elif event == "UN":
                COMMENT = "UNKNOWN"
            elif event.startswith("RA"):
                COMMENT = "Reanneal_{}".format(event.replace("RA_", ""))
                

            key = (chr, pos+1, strand, ref_base, softclip_bases, softclip_size, COMMENT)
            # key = (chr, pos+1, strand, ref_base, softclip_bases, softclip_size, COMMENT, FP)

            if key not in data["Count"]:
                data["Count"][key] = 0
            data["Count"][key] += 1
   
    df_out = pd.DataFrame(data)
    if df_out.shape[0] >0:
        # df_out.index.names = ["chr", "pos", "strand", "ref_base","softclip_bases", "Softclip_size", "COMMENT", "Five prime"]
        df_out.index.names = ["chr", "pos", "strand", "ref_base","softclip_bases", "Softclip_size", "COMMENT"]
        df_out = df_out.reset_index()
        # df_out = df_out[["chr", "pos", "strand", "ref_base","softclip_bases", "Count", "Softclip_size", "COMMENT", "Five prime"]]
        df_out = df_out[["chr", "pos", "strand", "ref_base","softclip_bases", "Count", "Softclip_size", "COMMENT"]]
        df_out.to_csv(options.output)
        
        print("[{t}] Total reads (Read1 only. Expected 1/2 as paired-end read counts): {T}".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime()), T=df_out["Count"].sum()))

    else:
        print("No output!!!")


