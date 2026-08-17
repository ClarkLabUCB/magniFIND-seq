# Trusted shell configuration for manuscript raw read/UMI QC.
PIPELINE_OUTPUT="../output/gene_expression"
OUT_DIR="../output/raw_qc"

# Fixed cutoffs reported for the K562-in-EL4 experiment. These are inclusive.
ASSIGNED_READ_THRESHOLD=478877
DEDUPLICATED_UMI_THRESHOLD=12149
ASSIGNED_READ_PERCENTILE=0.65
DEDUPLICATED_UMI_PERCENTILE=0.70
