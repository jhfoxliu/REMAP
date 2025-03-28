# A pipeline for amending of RNA 5' end mapping data

## Dependencies

python3
pysam == 


## Installation

Installation is not required. Please clone the scripts in this repo.

## Usage

1. Map the reads with aligners such as `hisat2` and `STAR`. I designed the pipeline with `hisat2`. Please make sure that softclipping is enabled. By defult, both `hisat2` and `STAR` enables softclipping by default. Please convert the `SAM` output into sorted and indexed `BAM`. 

2. Run the script `amend_bam.py` with the `BAM` file. 

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

The `amend_bam.py` will first go over the `BAM` input. For each read (read1 if paired), the script will amend the 5' end mapping result if softclipping occurs. Please refer to `amendment_logic.pdf` to see the decision tree. Once amendment is done, the script will generate a corrected `BAM` file (determined by `--bam-out`, `{input_bam}.corrected.bam` by defult). The output `BAM` file only contains the uniquely mapped read1 with the following tags:

Tags:

    SL  The length of softclipping. INT. 
    SS  The sequence of expended sequence. String. 
    ST  The type of expansion (see below). String. 
    PG  Genomic position (0-based). INT. 

The types of expansion (in the ST tag):

    N                   Non-expanded
    S                   SNP (only if `-s` is used)
    H                   Homopolymer from +1. e.g., A+1AAA
    H1                  Homopolymer from +2. e.g., C+1TTT
    H2                  Homopolymer from +3. e.g., T+1CAAAA
    UN                  Unknown. When it is unknown, the `SS` tag will return the sequence of the softclipping.
    RA_{seq1}_{seq2}    Ranneal event. {seq1} is the repeat motif, seq2 is the gapped sequence
    SK                  Skip (deletion) of the homopolymer


The script will then analyze the generated BAM file to perform statistics on the 5' end mapping.

The script will automatically reports the number of reads.

3. Once you get multiple `CSV` files, you can use the script `xxxx.py` to compare RNA 5' expansion among samples.


## Examples

![Fig SX algorithm](https://github.com/user-attachments/assets/1895b3a6-fb74-4e60-be6d-a6fa09bce335)

## Prototype

The algorithm in the `amend_bam.py` script is optimized based on the knowledge 

## Bug report

Please contact `Fox` (`jhfoxliu@gmail.com` or `jil4026@med.cornell.edu`) if you encounter any bug. Appriciate if you can post in `Issues`, so that the other users will know where is the bug and whether it get fixed.

## License

MIT.

## Citation

Paper is under submission.
