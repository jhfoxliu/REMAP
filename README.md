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

    N   Non-expanded
    S   SNP
    H: Homopolymer immediately
    H1: Homopolymer with first base overhang
    H2: Homopolymer with first tow bases overhang
    UN: Unknown
    RA_{seq1}_{seq2}: Ranneal. {seq1} is the repeat motif, seq2 is the gapped sequence
    SK: Skip


The script will then analyze the generated BAM file to perform statistics on the 5' end mapping.

3. Once you get multiple `CSV` files, you can use the script `xxxx.py` to compare RNA 5' expansion among samples.


## Examples

Please check the `Examples` folder for details:

## Prototype

The algorithm in the `amend_bam.py` script is optimized based on the knowledge 


## License

MIT.

## Citation

Paper is under submission.