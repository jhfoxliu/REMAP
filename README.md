# RNA 5' Expansion Mapping and Annotation Pipeline (REMAP)

It is belived that in eukaryotic cells, RNA Pol II synthesizes the RNA strand faithfully based on the template, and therefore the RNA 5' end should perfectly mathces the reference genome. However, I found that it is not true, that in many cases, there are extra nucleotides in the RNA 5' ends. I called this phenomenon "RNA 5' expansion", which is highly likely due to a dehybridization-reanneal-reinitation of transcription initiation. 

![Graphical Abstract_v3 - 1 - 1200x1200](https://github.com/user-attachments/assets/b497b6a9-4c26-4e8d-addb-d2494e4090fb)

Currently, all aligner cannot well handle RNA 5' expansion. For some cases, the aligner will return "soft clipping" (unaligned bases) at the 5' end of the reads; for some other instances, the aligner will report no soft clipping but 5' mismatches; for some extreme cases, when the expansion is too long, the aligner might generate an artifical splicing, where the expanded bases are mapped far away upstream of the TSS. 

To dissect RNA 5' expansion, I desgined `REMAP` to amend the 5' mapping results from `Hisat2`. This pipeline contains three steps:

(1) The script will go over the `BAM` file, and correct the mistakenly assigned 5' ends. This step will generate a amended `BAM` file.

(2) The script will go over the amended `BAM` file, and return a `CSV` file describing the amended results.

(3) I also provide a tool to find and compare 5' expansions among different samples.

## Tested dependencies

| Package   | Version |
| :---------| :-------|
| Pysam     | 0.19.1  |
| Biopython | 1.79    |
| Pandas    | 2.2.2   |

All based on `Python 3.9.14`.

## Installation

Installation is not required. Please clone the scripts in this repo.

## Usage

### 1. Mapping

Map the adapter-trimmed reads with aligners. I designed the pipeline based on `Hisat2` (v2.2.1). Please make sure that softclipping is enabled (default in `Hisat2`). Please convert the `SAM` output into sorted and indexed `BAM`. 

**Suggestions**

* Including SNPs in `Hisat2` indexes can improve mapping accuracy. Please refer to [`Hisat2 manual`](https://daehwankimlab.github.io/hisat2/manual/) about how to add SNPs. I suggest using UCSC common SNPs (MAF > 1%) rather than all SNPs, which might lead to fasle negatives.

* Normally I include both SNPs and splicing in the indexes. It can take up to 200 Gb RAM during index prepration and ~8 Gb RAM during mapping.

**Warnings**

* Please do not use `--no-temp-splicesite` option in `Hisat2`. Temporary splicing is very important for alignment when the expansion is extremely long!!!

* Unexpected errors might occur while using other aligners such as `STAR`, because the strategy in handling soft clipping might vary.

**Example mapping cmd (single-end)**:

`hisat2 -x {hisat2_index} --no-discordant --no-mixed -p {threads} -U read.fastq | samtools view -bS -@ {threads} -F 4 > hisat2.bam`

**Example mapping cmd (paired-end)**:

 `hisat2 -x {hisat2_index} --no-discordant --no-mixed --fr --rna-strandness FR -p {threads} -1 read1.fastq -2 read2.fastq | samtools view -bS -@ {threads} -F 4 > hisat2.bam`


### 2. Amend 5' alignments

Run the script `REMAP.v1.py` with the `BAM` file. 

Options:
    -b              The input `BAM` file.
    -f              The `FASTA` file of the reference genome.
    -p              The number of processors. Each processor will handle one chromosome.
    -o              The name of output file in `CSV` format.
    -s              `BCF` file. If provided, the program will mask the sites with known SNPs in the TSS.
    --bam-out       The name of output `BAM` file. `{input_bam}.corrected.bam` by defult.
    --skip-correct  Skip the BAM file generation. Will continue with the the amended `BAM` file.
    -p              Number of processors to use. One processor for one chromosome.
    --single-end    Please use this option if the reads are single-end.

The `REMAP.v1.py` will first go over the `BAM` input. For each read (read1 if paired), the script will amend the 5' end mapping result if softclipping occurs. Once amendment is done, the script will generate a corrected `BAM` file (determined by `--bam-out`, `{input_bam}.corrected.bam` by defult). The output `BAM` file only contains the uniquely mapped read1 with the following tags:

| Tag|  Explanation                       | Data type |
| :--| :----------------------------------| :-------  |
| SL |  The length of softclipping.       | INT       |
| SS |  The sequence of expended sequence.| String    |
| ST |  The type of expansion (see below).| String    |
| PG |  Genomic position (0-based).       | INT       |

The types of expansion (in the ST tag):

| BAM file value     | CSV file value (COMMENTS column) | Explanation  |
| :------------------| :--------------------------------| :----------  |
|   M                |  M                               | Non-expanded |
|   S                |  SNP                             | SNP (only if `-s` is used) | 
|   H                |  +1 EXPANDED                     | Homopolymer from +1. e.g., A+1AAA | 
|   H1               |  +2 EXPANDED                     | Homopolymer from +2. e.g., C+1TTT |
|   H2               |  +3 EXPANDED                     | Homopolymer from +3. e.g., T+1CAAAA |
|   H3               |  +4 EXPANDED                     | Homopolymer from +4. e.g., G+1T1CAAAA |
|   UN               |  UNKNOWN                         | Unknown. When it is unknown, the `SS` tag will return the sequence of the softclipping. |
|   RA_{seq1}_{seq2} |  Reanneal_{seq1}_{seq2}          | Ranneal event. {seq1} is the repeat motif, seq2 is the gapped sequence |
|   SK               |  SKIPPED                         | Skip (deletion) of the homopolymer |


The script will then analyze the generated BAM file to perform statistics on the 5' end mapping.

The script will automatically reports the number of reads.

### 3. Compare different samples

Once you get multiple `CSV` files, you can use the script `REMAP.processing.v1.py` to compare RNA 5' expansion among samples.


## Mechanism

<img width="2255" height="570" alt="Examples" src="https://github.com/user-attachments/assets/66d0be44-98f3-4500-aa6c-b2c1f55c38ef" />

## Prototype

The prototype relies on removing the last base from the mapped 5' ends until the trimmed reads are fully mapped. This strategy is very slow and no longer used.

## Bug report

Please contact `Fox` (`jhfoxliu@gmail.com` or `jil4026@med.cornell.edu`) if you encounter any bug. Appriciate if you can post in `Issues`, so that the other users will know where is the bug and whether it get fixed.

## License

MIT.

## Citation

The paper is under submission.
