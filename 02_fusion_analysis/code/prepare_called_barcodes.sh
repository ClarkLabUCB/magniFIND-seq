#!/usr/bin/env bash
set -euo pipefail

[[ $# -le 1 ]] || { echo "Usage: $0 [output-directory]" >&2; exit 2; }
out_dir="${1:-$PWD/public_k562_inputs}"
[[ -n "$out_dir" && "$out_dir" != / ]] || { echo "[ERROR] Unsafe output directory" >&2; exit 2; }
mkdir -p "$out_dir"
compressed="$out_dir/GSM2406675_10X001_barcodes.tsv.gz"
barcodes="$out_dir/GSM2406675_called_barcodes.tsv"
temporary="${compressed}.tmp.$$"
trap 'rm -f "$temporary" "${barcodes}.tmp.$$"' EXIT
url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2406675&file=GSM2406675_10X001_barcodes.tsv.gz&format=file'
expected_sha256='b44a3d466b03cafe7479ed98b3f6aa9becc9ce3651f17b49f7f77ca03d72cd32'

curl -L --fail --retry 3 -o "$temporary" "$url"
actual_sha256="$(sha256sum "$temporary" | awk '{print $1}')"
[[ "$actual_sha256" == "$expected_sha256" ]] || {
  echo "[ERROR] GEO barcode checksum mismatch: $actual_sha256" >&2
  exit 1
}
gzip -t "$temporary"
gzip -cd "$temporary" > "${barcodes}.tmp.$$"
awk '
  $0 !~ /^[ACGTN]+-[0-9]+$/ { exit 2 }
  { count++ }
  END { if (count != 5768) exit 3 }
' "${barcodes}.tmp.$$"
mv "$temporary" "$compressed"
mv "${barcodes}.tmp.$$" "$barcodes"
trap - EXIT
printf '[DONE] %s called barcodes: %s\n' "$(wc -l < "$barcodes" | tr -d ' ')" "$barcodes"
