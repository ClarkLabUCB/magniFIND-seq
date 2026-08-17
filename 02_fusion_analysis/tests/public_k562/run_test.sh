#!/usr/bin/env bash
set -euo pipefail
analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
seq96="$(printf 'ACGT%.0s' {1..24})"
reverse_seq96="$(printf 'A%.0s' {1..96})"
qual96="$(printf 'I%.0s' {1..96})"
{
  printf '@HD\tVN:1.6\tSO:unsorted\n'
  printf 'read1\t0\tchr22\t100\t255\t96M\t*\t0\t0\t%s\t%s\tCB:Z:CELL1\tUB:Z:U1\tCR:Z:RAW1\tUR:Z:RU1\n' "$seq96" "$qual96"
  printf 'read2\t16\tchr22\t100\t255\t96M\t*\t0\t0\t%s\t%s\tCB:Z:CELL1\tUB:Z:U2\tCR:Z:RAW1\tUR:Z:RU2\n' "$reverse_seq96" "$qual96"
  printf 'read3\t0\tchr9\t300\t255\t96M\t*\t0\t0\t%s\t%s\tCB:Z:CELL2\tUB:Z:U3\tCR:Z:RAW2\tUR:Z:RU3\n' "$seq96" "$qual96"
  printf 'secondary\t256\tchr9\t300\t255\t96M\t*\t0\t0\t%s\t%s\tCB:Z:CELL2\tUB:Z:BAD\n' "$seq96" "$qual96"
} > "$tmp/input.sam"
"$analysis_dir/code/public_k562/extract_primary_transcripts.py" \
  --input "$tmp/input.sam" --input-format sam \
  --fastq "$tmp/reads.fastq.gz" --metadata "$tmp/metadata.tsv"
gzip -t "$tmp/reads.fastq.gz"
[[ "$(gzip -cd "$tmp/reads.fastq.gz" | wc -l | tr -d ' ')" == 12 ]]
python3 - "$tmp/reads.fastq.gz" <<'PY'
import gzip, sys
with gzip.open(sys.argv[1], "rt") as handle:
    lines = [line.rstrip() for line in handle]
assert lines[5] == "T" * 96  # reverse-strand SAM sequence restored to read orientation
PY
cat > "$tmp/genes.gtf" <<'EOF'
chr22	test	gene	90	200	.	+	.	gene_id "BCR_ID"; gene_name "BCR";
chr9	test	gene	290	400	.	+	.	gene_id "ABL1_ID"; gene_name "ABL1";
chr1	test	gene	490	600	.	-	.	gene_id "G1_ID"; gene_name "GENE1";
chr2	test	gene	690	800	.	+	.	gene_id "G2_ID"; gene_name "GENE2";
EOF
python3 - "$tmp/metadata.tsv" "$tmp/junctions.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as h:
    names = {r["original_qname"]: r["star_read_name"] for r in csv.DictReader(h, delimiter="\t")}
rows = [
    ("chr22", 100, "+", "chr9", 300, "+", names["read1"]),
    ("chr22", 101, "+", "chr9", 301, "+", names["read1"]),  # same physical read
    ("chr22", 100, "+", "chr1", 500, "-", names["read2"]),
    ("chr9", 300, "+", "chr2", 700, "+", names["read3"]),
]
with open(sys.argv[2], "w") as h:
    for a, p1, s1, b, p2, s2, name in rows:
        h.write(f"{a}\t{p1}\t{s1}\t{b}\t{p2}\t{s2}\t-1\t0\t0\t{name}\t1\t48M48S\t1\t48S48M\n")
PY
printf 'barcode\nCELL1\nCELL2\n' > "$tmp/called.tsv"
"$analysis_dir/code/public_k562/summarize_chimeric_categories.py" \
  --junctions "$tmp/junctions.tsv" --gtf "$tmp/genes.gtf" \
  --read-metadata "$tmp/metadata.tsv" --called-barcodes "$tmp/called.tsv" \
  --out-dir "$tmp/output"
python3 - "$tmp/output/category_summary.tsv" <<'PY'
import csv, math, sys
with open(sys.argv[1], newline="") as h:
    rows = {r["category"]: r for r in csv.DictReader(h, delimiter="\t")}
assert rows["whole_intergene"]["physical_reads_all"] == "3"
assert rows["BCR_any"]["physical_reads_all"] == "2"
assert rows["ABL1_any"]["physical_reads_all"] == "2"
assert rows["BCR_ABL1"]["physical_reads_all"] == "1"
assert rows["BCR_ABL1"]["fusion_cb_umi_pairs"] == "1"
assert rows["BCR_ABL1"]["total_called_cb_umi_pairs"] == "3"
assert math.isclose(float(rows["BCR_ABL1"]["fusion_umi_percent"]), 100 / 3)
assert rows["BCR_ABL1"]["positive_called_cells"] == "1"
assert rows["BCR_ABL1"]["total_called_cells"] == "2"
assert float(rows["BCR_ABL1"]["positive_cell_percent"]) == 50
PY
echo "[PASS] Public K562 chimeric fixture passed"
