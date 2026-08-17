# Analysis methods

This document describes the reference implementation in this repository. Exact
software versions and parameters for a completed run are also written into its
`run_metadata/` and `fusion_scan/` output directories.

## Gene-expression quantification

Paired-end magniFIND-seq FASTQ files were supplied with one manifest row per
sample or cell. Prior to processing, both files in every pair were read in full
to validate FASTQ structure, equal record counts, sequence/quality lengths, and
matching paired read identifiers. Sample identifiers were checked for uniqueness
and filesystem-safe characters.

For UMI-containing reads, read 1 was searched for an exact match to the
Smart-seq3 template-switch oligonucleotide sequence `TTGCGCAATG`. The eight
nucleotides immediately following this sequence were extracted as the UMI. The
template-switch sequence and UMI were removed from read 1, and both read names
were annotated with sample and UMI fields in the form
`|CB:<sample>|UB:<UMI>`. Read pairs lacking an exact template-switch match were
written to the non-UMI pair. Matching the production implementation, any
non-empty substring in the eight-base UMI window was written to the UMI-positive
pair, including a short UMI at the end of a read; an empty post-UMI read-1
sequence was also retained. Both UMI and non-UMI pairs were preserved.

Both output pairs were passed to Trim Galore with the production options
`--illumina --paired --fastqc`; only the trimmed UMI-positive pair proceeded to
alignment and quantification. Trimmed UMI-positive pairs were aligned with
production STAR 2.7.2a to a user-specified genome index, retaining attributes
`NH HI AS nM XS`. For gene quantification, alignments were
restricted to MAPQ 255, indicating uniquely mapped STAR alignments. Sample and
UMI values encoded in read names were copied to CB and UB BAM tags without
changing the number of alignment records.

Paired-end reads were assigned to GTF exon features at the gene level with
featureCounts 2.0.6 using `-p -B -R BAM`, with `gene_id` as the default grouping
attribute. The production command did not include `--countReadPairs`.
Gene-assigned alignments were coordinate sorted and UMIs were collapsed with
UMI-tools using its method default and the featureCounts XT/XS, CB, and UB tags;
the production command did not include `--method` or `--paired`. Per-sample UMI
tables were merged in manifest order into a sparse gene-by-sample matrix. Gene ID, gene name, and
gene biotype metadata were parsed from the same GTF.

The sparse count matrix was imported into R 4.3.3 using Matrix 1.6-5 and stored
as a Seurat 5.0.3 object. No minimum-gene, minimum-cell, mitochondrial-content,
or biological doublet filters were imposed during object construction. Sample
identifiers, optional sample-level metadata, assay-level gene ID/name/type
metadata, software provenance, and R session information were retained in the
object.

## Bulk-to-single-cell differential-expression concordance

This workflow started from a raw integer bulk gene-count matrix. Alignment and
gene quantification of the upstream bulk single-end FASTQs were outside this
module; all complete biological-replicate libraries must be quantified against
the same documented reference and annotation before running this workflow.
Downsampled FASTQs are not suitable for reproducing the differential-expression
analysis.

Each condition was required to contain at least the configured number of
biological replicates (default two; the study design contains three).

Raw bulk RNA-seq gene counts and an explicit two-column sample-to-condition
table can be reanalyzed with the pinned edgeR 3.42.4 diagnostic workflow. Genes are retained with
`filterByExpr`, library composition was normalized by the trimmed mean of M
values, and biological dispersion was estimated with robust estimation. A
two-condition generalized linear model was fitted with the edgeR quasi-
likelihood framework, and condition A versus condition B was tested with a
quasi-likelihood F test. Differentially expressed genes were defined by the
configured thresholds, which default to FDR < 0.05 and absolute log2 fold
change >= 1.

Untransformed UMI counts were obtained from the configured assay and count
layer of a Seurat object. A predefined logical QC metadata field was applied
when configured. Counts from the two comparison groups were
summed across cells within each condition and converted to counts per million.
Descriptive single-cell pseudobulk log2 fold change was calculated as
`log2(CPM_A + 1) - log2(CPM_B + 1)` by default. Per-condition cellular
detection rates were calculated as the fraction of retained cells with a
nonzero count. This pseudobulk step is descriptive and does not model
single-cell biological replication.

A fresh Seurat object was then constructed from the same retained raw counts,
log-normalized at scale factor 10,000, and analyzed with `FindMarkers` for
condition A versus condition B. The default uses the Wilcoxon test with
`logfc.threshold=0` and `min.pct=0`; these settings are recorded and
configurable. The `avg_log2FC` (or version-equivalent `avg_logFC`) values from
this analysis, rather than the descriptive pseudobulk fold changes, are used
for manuscript concordance.

For exact manuscript signature reproduction, ranking starts from the archived
production differential-expression table generated with edgeR 4.0.16. Its
reported GroupB/GroupA logFC is multiplied by -1 to give K562-R versus K562-S.
Recomputing with edgeR 3.42.4 is retained as a transparent sensitivity analysis
but is not substituted for that production table.

Bulk genes were matched to single-cell features by gene symbol. For duplicated
bulk symbols, the entry with the lowest FDR and then largest absolute bulk
log2 fold change was retained. Candidate signature genes had to have at least
one single-cell count. Production selection required FDR < 0.05 and absolute
log2 fold change >= 1 and excluded `^MT-` and `^RP[SL]` gene symbols. Genes
were ordered by FDR and then absolute log2 fold change,
and the top 50 condition-A-high and top 50 condition-B-high genes were selected
independently. Thus, “top 50 per direction” denotes 100 genes in the combined
signature rather than 50 genes total.

For manuscript reproduction, failure to recover all 50 genes in either
direction is a hard error rather than a shortened signature.

Bulk and single-cell FindMarkers log2 fold changes were compared across the
combined selected signature using Pearson correlation, Spearman correlation,
direction concordance, and ordinary least-squares regression. The same
descriptive statistics were recorded separately for each direction, all
matched expressed genes, and all eligible bulk DEGs. Because genes are selected
using the bulk test before correlation, the selected-gene Pearson P value is
conditioned on feature selection and is not an independent inferential test.
Genes are also not statistically independent observations. A fully confounded
condition and sequencing-batch design cannot distinguish biological condition
effects from batch effects and should not be used for causal claims.

For the bulk-derived resistance score used in the original K562 analysis, retained raw
single-cell counts were independently normalized with Seurat LogNormalize
(scale factor 10,000). The selected condition-A-high and condition-B-high
genes were passed as separate feature sets to `AddModuleScore`. The default
parameters were `ctrl=50`, `nbin=24`, and `seed=20260729`. Resistance score was
defined per cell as the condition-A-high module score minus the
condition-B-high module score. When condition A is K562-R and condition B is
K562-S, this reproduces the original R-high minus S-high definition. This score
is also calculated for other QC-retained labels, such as PBMC, without using
those groups in the A-versus-B differential-expression test. This score
is descriptive and depends on the selected bulk signature and normalization
context; it is not a clinically validated resistance classifier.

## Fusion-aware mapping and BCR::ABL1 screening

Fusion analysis was performed independently from gene-expression mapping using
the UMI-trimmed paired FASTQs. STAR 2.7.2a was run without two-pass mapping (the
production `FUSION_TWOPASS_MODE=None`) with
`--chimSegmentMin 12`, `--chimJunctionOverhangMin 12`,
`--chimOutType Junctions WithinBAM SoftClip`, `--chimOutJunctionFormat 1`,
`--chimMainSegmentMultNmax 1`, and
`--outSAMunmapped Within KeepPairs`. Alignment attributes
`NH HI AS nM NM MD XS ch` were retained. The coordinate-sorted BAM was not
filtered by mapping quality or alignment flag; unmapped, low-MAPQ,
supplementary, split, and soft-clipped evidence was preserved. A second BAM with
CB/UB tags was created when tags were present in query names, and record counts
were required to equal those in the source BAM.

STAR `Chimeric.out.junction` data records were screened for one coordinate in
the configured BCR interval and the other in the configured ABL1 interval,
accepting either STAR report order. Records were canonicalized to BCR followed
by ABL1 and retained as candidates only when both reported strands matched the
configured gene strands. The reference configuration uses GRCh38-style
intervals BCR `chr22:23179704-23318037 (+)` and ABL1
`chr9:130713946-130887675 (+)`; these values must be changed when another
assembly or contig convention is used.

Evidence was summarized per sample as supporting chimeric records, distinct
read names, distinct literal UB strings, and distinct canonical breakpoint
tuples. The manuscript binary candidate field is one when at least one
expected-orientation record is present. As an additional evidence tier, samples
with no matching records were called `not_detected`; samples
with evidence from fewer than two distinct UBs were called `candidate_low_umi`;
and samples with at least two distinct UBs were called
`candidate_ge_min_umi`. Missing, incomplete, or malformed mapping output was
reported as `not_evaluated` rather than negative. The UB summary does not
perform sequence-error correction and should not be interpreted as a validated
molecule count.

## Raw read- and UMI-based single-bead QC

The author-supplied production `analyze_pipeline.py` generated a long-format
read-tracking table containing raw UMI/non-UMI, trimmed UMI/non-UMI, assignment,
tagged-read, tagged-UMI, and deduplicated-UMI metrics. Its raw FASTQ filename
classifier counted both R1 and R2 despite a comment stating that only R1 should
be counted. The repository therefore emits both
`readtracking_production_legacy.tsv`, which preserves that behavior, and
`readtracking_corrected_r1.tsv`, which counts only R1. They must not be mixed.

For the K562-in-EL4 experiment only, featureCounts' `Assigned` value from the
production `-p -B` command was read from each completed gene-pipeline summary
and the total deduplicated UMI count was calculated by summing the corresponding
UMI-tools gene table. Samples were retained when both values were less than or
equal to 478,877 and 12,149, respectively. Equality therefore passes. The
sample-distribution 65th percentile for assigned alignments and 70th percentile
for deduplicated UMIs were additionally calculated with the R-compatible type-7
definition and written beside the fixed applied values. Newly computed
percentiles do not silently replace the manuscript cutoffs.

## Public K562 intergenic chimeric transcripts

For SRR5082088, secondary and supplementary input alignments were excluded with
SAM flag mask `0x900`. Sequence and quality were restored to original read
orientation, the first 96 transcript bases were written as single-end FASTQ,
and CB, UB, CR, and UR were preserved in an explicit read metadata table and
URL-escaped STAR query name. The reads were realigned to a user-supplied
Cell Ranger GRCh38-2024-A-compatible STAR index.

Both coordinates in every STAR chimeric junction were intersected with GENCODE
v44 `gene` features. A record was intergenic when at least one annotation pair
assigned different gene IDs to the ends. Orientation was ignored when assigning
`whole_intergene`, `BCR_any`, `ABL1_any`, and `BCR_ABL1`. Repeated junctions and
overlapping annotations were deduplicated by original QNAME within category.
Fusion-associated CB-UMI pairs and positive cells were restricted to the
provided called-cell barcode set; their explicit denominators were all unique
called-cell CB-UMI pairs recovered from primary reads and all called barcodes.

## PBMC and K562-r embedding

Public 10x PBMC counts and magniFIND-seq K562-r counts were converted to a
common gene-symbol space and merged without batch correction. The merged raw
counts were log-normalized, 3,000 variable features selected, scaled, reduced by
30-component PCA, and embedded with UMAP using dimensions 1-20. The K562-r-only
object used 2,000 variable features, 30-component PCA, and UMAP dimensions 1-15.
Both used seed 20260730 and a sequential worker plan. For
each Seurat input, the raw-count assay was explicitly named in the input
manifest rather than inferred from the object's mutable default assay. Counts
were required to be finite, non-negative integers. Fusion-candidate status was joined
only after sample identity was retained and was not used for filtering or any
dimensionality-reduction step; every K562-r sample was required to have an
explicit binary candidate entry. A separate K562-r-only object was normalized and
embedded de novo rather than reusing coordinates from the joint UMAP.

Fusion results represent computational screening for transcript evidence.
Positive findings require orthogonal confirmation, and a `not_detected` result
does not establish biological absence when expression or sequencing depth is
insufficient.

Pipeline restart fingerprints were calculated from SHA-256 file-content hashes,
effective options, relevant source files, and tool versions. Fusion output was
fully assembled and validated in a sample-specific temporary directory before
replacement of any previous sample directory. Configuration snapshots,
effective parameters, source-control state, and source snapshot hashes were
retained with completed outputs.
