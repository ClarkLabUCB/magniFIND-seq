# Trusted configuration for the public SRR5082088/GSM2406675 analysis.
INPUT_BAM="SRR5082088.bam"
CALLED_BARCODES="called_barcodes.tsv"
GTF_FILE="reference/genes.gtf"
STAR_INDEX="reference/star_index"
OUT_DIR="../output/SRR5082088_chimeric"

STAR="STAR"
SAMTOOLS="samtools"
PYTHON="python3"
THREADS=8
READ_LENGTH=96
CHIM_SEGMENT_MIN=12
CHIM_JUNCTION_OVERHANG_MIN=12
CHIM_MAIN_SEGMENT_MULT_NMAX=1
