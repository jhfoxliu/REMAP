import sys
import pysam

with pysam.AlignmentFile(sys.argv[1], "rb") as bam_in, pysam.AlignmentFile(sys.argv[2], "wb", template=bam_in) as bam_out:
    for read in bam_in:
        read_name = read.query_name
        if read_name.endswith("ATAT"):
            bam_out.write(read)
