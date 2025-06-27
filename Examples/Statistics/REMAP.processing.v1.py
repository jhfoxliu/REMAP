import multiprocessing
from multiprocessing import Process,Pool
import time
import signal
import sys,os
from Bio import SeqIO
from Bio.Seq import reverse_complement
import pysam
import psutil
from itertools import chain
import scipy.stats
import numpy as np
import pandas as pd
import argparse
from time import gmtime, strftime
pid = os.getpid()

def process_file(fn, name, level, cov):
    global pid
    allowed = {}
    for i in range(1, 12):
        allowed["A" * i] = 1
        allowed["T" * i] = 1
        allowed["C" * i] = 1
        allowed["G" * i] = 1

    subdf = pd.read_csv(fn, low_memory=False)
    # EXPANDED = subdf[(subdf["COMMENT"]=="EXPANDED") | (subdf["COMMENT"]=="+2 EXPANDED") | (subdf["COMMENT"]=="+3 EXPANDED")]
    subdf = subdf[(subdf["COMMENT"]=="+1 EXPANDED") | (subdf["COMMENT"]=="+2 EXPANDED") | (subdf["COMMENT"]=="+3 EXPANDED") | (subdf["COMMENT"]=="M")]
    # subdf = subdf[subdf["softclip_bases"].isin(allowed) == True]
    # exclude = subdf[(subdf["softclip_bases"] != subdf["ref_base"]) & (subdf["Softclip_size"] == 0)]
    # subdf = subdf.loc[subdf.index.difference(exclude.index)]
    subdf_sc = pd.pivot_table(data=subdf, index=["chr", "pos", "strand"], columns=["Softclip_size"], values="Count", aggfunc="sum").fillna(0)
    subdf_count_sum = subdf_sc.sum(axis=1)
    subdf_sc = subdf_sc.loc[(subdf_count_sum[subdf_count_sum>= cov]).index]
    subdf_sc_weighted = subdf_sc*subdf_sc.columns
    subdf_sc_level = pd.DataFrame(subdf_sc_weighted.sum(axis=1)) 
    subdf_sc_level.columns = [name]
    subdf_sc_level["sum"] = subdf_count_sum.loc[subdf_sc_level.index]
    subdf_sc_level[name] = subdf_sc_level[name] / subdf_sc_level["sum"]
    subdf_sc_level = subdf_sc_level[[name]]
    subdf_sc_level = subdf_sc_level[subdf_sc_level[name]>=level]
    subdf_sc_level.to_csv("temp_by_single_sample_{}_{}.csv".format(pid, name))

    subdf_frac_exp = subdf[(subdf["COMMENT"]=="+1 EXPANDED") | (subdf["COMMENT"]=="+2 EXPANDED") | (subdf["COMMENT"]=="+3 EXPANDED")]
    subdf_frac_non = subdf[(subdf["COMMENT"]=="M")]

    subdf_frac_exp = subdf_frac_exp.groupby(by=["chr", "pos", "strand"])[["Count"]].sum()
    subdf_frac_non = subdf_frac_non.groupby(by=["chr", "pos", "strand"])[["Count"]].sum()
    subdf_frac = pd.concat([subdf_frac_exp, subdf_frac_non], axis=1)
    subdf_frac.columns = ["exp", "non_exp"]
    subdf_frac = subdf_frac.fillna(0)
    subdf_frac["ratio"] = subdf_frac["exp"] / (subdf_frac["exp"] + subdf_frac["non_exp"])
    subdf_frac = subdf_frac[["ratio"]]
    subdf_frac.columns = [name]
    subdf_frac.to_csv("temp_by_single_sample_fraction_{}_{}.csv".format(pid, name))


def process_file_all(fn, name, site_list, cov):
    global pid
    allowed = {}
    for i in range(1, 20):
        allowed["A" * i] = 1
        allowed["T" * i] = 1
        allowed["C" * i] = 1
        allowed["G" * i] = 1

    subdf = pd.read_csv(fn, low_memory=False)
    subdf = subdf[subdf["softclip_bases"].isin(allowed) == True]
    exclude = subdf[(subdf["softclip_bases"] != subdf["ref_base"]) & (subdf["Softclip_size"] == 0)]
    subdf = subdf.loc[subdf.index.difference(exclude.index)]
    subdf_sc = pd.pivot_table(data=subdf, index=["chr", "pos", "strand"], columns=["Softclip_size"], values="Count", aggfunc="sum").fillna(0)
    subdf_count_sum = subdf_sc.sum(axis=1)
    subdf_sc = subdf_sc.loc[(subdf_count_sum[subdf_count_sum>= cov]).index]
    subdf_sc = subdf_sc.loc[subdf_sc.index.intersection(site_list)]
    subdf_sc_weighted = subdf_sc*subdf_sc.columns
    subdf_sc_level = pd.DataFrame(subdf_sc_weighted.sum(axis=1))
    subdf_sc_level.columns = [name]
    subdf_sc_level["sum"] = subdf_count_sum.loc[subdf_sc_level.index]
    subdf_sc_level[name] = subdf_sc_level[name] / subdf_sc_level["sum"]
    subdf_sc_level = subdf_sc_level[[name]]
    subdf_sc_level.to_csv("temp_by_all_sites_{}_{}.csv".format(pid, name))
    subdf_sc_cov = pd.DataFrame(subdf_count_sum).loc[subdf_sc_level.index]
    subdf_sc_cov.columns = [name]
    subdf_sc_cov.to_csv("temp_coverage_by_all_sites_{}_{}.csv".format(pid, name))

if __name__ == "__main__":
    #Parser
    parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',formatter_class=argparse.RawTextHelpFormatter)
    group_required = parser.add_argument_group("Required")
    group_required.add_argument("-i",dest="input",help="Input list")
    group_required.add_argument("-o","--output",dest="output",help="Output prefix",default="[prefix]_[cov_cutoff]_[level_cutoff].csv")
    group_required.add_argument("-P",dest="process",type=int,default=4,help="Process number, default=4")
    group_required.add_argument("-c","--coverage",dest="coverage", type=int, default=50, help="Coverage cutoff for site filtering, default=50")
    group_required.add_argument("--c2","--coverage2",dest="coverage2", type=int, default=50, help="Coverage cutoff for concaternating, default=50")
    group_required.add_argument("-l","--level",dest="minimum_level", type=float, default=0.1, help="Cutoff for extra base levels, default=0.1")
    group_required.add_argument("--bedtools",dest="bedtools", default="/home/fox/Software/bedtools/2.29.1/bedtools", help="Bedtools directory")
    group_required.add_argument("-f","--fasta",dest="fasta",help="Fasta file containing the names", default="/home/fox/Database/Genome/Human/Gencode/v45/GRCh38.primary_assembly.genome.nochr.fa")
    group_required.add_argument("-t","--tss",dest="tss_db",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.tss.bed",help="")
    group_required.add_argument("-e","--exons",dest="exons_db",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.tx.bed",help="")
    group_required.add_argument("-a","--annot",dest="annot",required=False, default="/home/fox/Database/Genome/Human/Gencode/v45/gencode.v45.primary_assembly.annotation.nochr.anno",help="")

    options = parser.parse_args()
    print("[{t}] Analysis begins.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    bedtools = options.bedtools

    sample_key_values = []
    column_order = []
    with open(options.input, "r") as input:
        for line in input.readlines():
            line = line.strip().split("\t")
            sample_name, fn = line
            sample_key_values.append((sample_name, fn))
            column_order.append(sample_name)

    def signal_handler(sig,frame):
        pool.terminate()
        sys.exit()
    
    signal.signal(signal.SIGINT,signal_handler)

    print("[{t}] Looking for extra bases...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    pool = multiprocessing.Pool(options.process)
    
    try:
        for key, value in sample_key_values:
            pool.apply_async(process_file, args=(value, key, options.minimum_level, options.coverage, ))
        pool.close()
        pool.join()
        dfs = []
        dfs_frac = []
        for fn in os.listdir("./"):
            if fn.startswith("temp_by_single_sample_{}".format(pid)):
                dfs.append(pd.read_csv(fn, index_col=[0, 1, 2], header=0, low_memory=False))
            elif fn.startswith("temp_by_single_sample_fraction_{}".format(pid)):
                dfs_frac.append(pd.read_csv(fn, index_col=[0, 1, 2], header=0, low_memory=False))
        df_merged = pd.concat(dfs, axis=1)
    finally:
        pool.terminate()
        for fn in os.listdir("./"):
            if fn.startswith("temp_by_single_sample_{}".format(pid)):
                os.remove(fn)
            elif fn.startswith("temp_by_single_sample_fraction_{}".format(pid)):
                os.remove(fn)

    print("[{t}] Merging tables...".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    pool = multiprocessing.Pool(options.process)
    
    try:
        for key, value in sample_key_values:
            pool.apply_async(process_file_all, args=(value, key, df_merged.index, options.coverage, ))
        pool.close()
        pool.join()
        dfs = []
        for fn in os.listdir("./"):
            if fn.startswith("temp_by_all_sites_{}".format(pid)):
                dfs.append(pd.read_csv(fn, index_col=[0, 1, 2], header=0, low_memory=False))
        df_merged = pd.concat(dfs, axis=1)
    finally:
        pool.terminate()
        for fn in os.listdir("./"):
            if fn.startswith("temp_by_all_sites_{}".format(pid)):
                os.remove(fn)
        dfs_cov = []
        for fn in os.listdir("./"):
            if fn.startswith("temp_coverage_by_all_sites_{}".format(pid)):
                dfs_cov.append(pd.read_csv(fn, index_col=[0, 1, 2], header=0, low_memory=False))
                os.remove(fn)
        df_cov_merged = pd.concat(dfs_cov, axis=1)
    
    df_merged = df_merged[column_order]
    df_cov_merged = df_cov_merged[column_order]
    
    print("[{t}] Data fetching ends. Start annotation.".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # annotate sequence
    def get_flanking_pysam(ref_genome,chr,pos,strand,flanking=None,up=0,down=0):
        pos_0 = pos - 1
        if flanking is not None:
            up = flanking
            down = flanking
        try:
            if strand == "+":
                seq = ref_genome.fetch(reference=chr, start=pos_0-up, end=pos_0+1+down) # [start, end) 0-based
                if seq and len(seq) == 1+up+down:
                    return seq.upper()
            elif strand == "-":
                seq = ref_genome.fetch(reference=chr, start=pos_0-down, end=pos_0+1+up)
                seq = reverse_complement(seq)
                if seq and len(seq) == 1+up+down:
                    return seq.upper()
        except ValueError:
            pass
        except IndexError:
            pass
        except KeyError:
            if options.nona == True:
                pass

    
    with pysam.FastaFile(options.fasta) as ref_genome:
        df_merged["N20"] = df_merged.apply(lambda x: get_flanking_pysam(ref_genome, x.name[0], int(x.name[1]), x.name[2], up=0, down=10), axis=1)
        df_merged["N10"] = df_merged["N20"].str[0:10]
        df_merged["N4"] = df_merged["N20"].str[0:4]

    # annotate gene
    with open(str(pid) + ".bed", "w") as bed_out:
        for idx, row in df_merged.iterrows():
            bed_out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(idx[0], idx[1] - 1, idx[1], ".", ".", idx[2]))

    # Bedtools closest
    cmd = "{bedtools} sort -i {bed_in} > {bed_out}".format(bedtools = bedtools, bed_in = str(pid) + ".bed", bed_out = str(pid) + ".sorted.bed")
    print(cmd)
    os.system(cmd) 
    cmd = "{bedtools} closest -s -D b -a {bed_out} -b {database} > {closest_tss}".format(bedtools=bedtools, bed_out = str(pid) + ".sorted.bed", closest_tss=str(pid) + ".closest_tss.bed", name=str(pid), database=options.tss_db)
    os.system(cmd)
    print(cmd)
    cmd = "{bedtools} closest -S -D b -a {bed_out} -b {database} > {closest_tss}".format(bedtools=bedtools, bed_out = str(pid) + ".sorted.bed", closest_tss= str(pid) + ".closest_tss.upstream.bed", name=str(pid), database=options.tss_db)
    os.system(cmd)
    print(cmd)
    cmd = "{bedtools} closest -s -D b -a {bed_out} -b {database} > {closest_exons}".format(bedtools=bedtools, bed_out = str(pid) + ".sorted.bed", closest_exons= str(pid) + ".closest_exons.bed", name=str(pid), database=options.exons_db)
    os.system(cmd)
    print(cmd)

    # Merge annotations
    os.system("python3 collaspe_bed_annotations_v3.py --thresh-min -1000 -i {closest} -a {annot} -o {collasped}".format(closest= str(pid) + ".closest_tss.bed", annot=options.annot, collasped = str(pid) + ".collasped.1.bed"))

    os.system("python3 collaspe_bed_annotations_fix_other_exons_v2.py -a {annot} -i {collasped} -I {closest_exons} -o {collasped_2} ".format(collasped = str(pid) + ".collasped.1.bed", closest_exons= str(pid) + ".closest_exons.bed", annot=options.annot, collasped_2= str(pid) + ".collasped.2.bed"))

    os.system("python3 collaspe_bed_annotations_fix_upstream_TSS.py  --thresh-min -2000 --thresh-max 2000 -a {annot} -i {collasped} -I {closest_exons} -o {collasped_final} ".format(collasped = str(pid) + ".collasped.2.bed", closest_exons= str(pid) + ".closest_tss.upstream.bed", annot=options.annot, collasped_final= str(pid) + ".collasped_final.bed"))


    data_annot_gene = {}
    data_annot_biotype = {}
    data_annot_method = {}
    data_annot_dist = {}

    with open(str(pid) + ".collasped_final.bed", "r") as annotated_bed:
        for line in annotated_bed.readlines():
            line = line.strip().split("\t")
            key = (line[0], int(line[2]), line[5])
            annot = line[3]
            gene, gene_biotype, trans_biotype, method, dist = line[3].split("|")

            data_annot_gene[key] = gene
            data_annot_biotype[key] = gene_biotype
            data_annot_method[key] = method
            data_annot_dist[key] = dist

    os.system("rm {}".format(str(pid) + "*.bed"))

    df_merged["Gene"] = df_merged.apply(lambda x: data_annot_gene.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df_merged["Gene_biotype"] = df_merged.apply(lambda x: data_annot_biotype.get((x.name[0], x.name[1], x.name[2])), axis=1)
    df_merged["Annotate_to_closest"] = df_merged.apply(lambda x: data_annot_method.get((x.name[0], x.name[1], x.name[2])), axis=1)

    dfs_frac_merged = pd.concat(dfs_frac, axis=1)


    df_merged.to_csv("{}_extrabases_cov{}_level{}.csv".format(options.output, options.coverage2, options.minimum_level))
    df_cov_merged.to_csv("{}_coverages_cov{}_level{}.csv".format(options.output, options.coverage2, options.minimum_level))
    dfs_frac_merged.loc[df_cov_merged.index.intersection(df_merged.index)].to_csv("{}_fractions.csv".format(options.output, options.coverage2, options.minimum_level))
    print("[{t}] Done!".format(t=strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
