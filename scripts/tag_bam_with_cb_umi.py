#!/usr/bin/env python3
"""Copy CB/UB values encoded in query names into BAM tags."""

import argparse
import re

import pysam


CB_RE = re.compile(r"(?:^|\|)CB:([^|]+)")
UB_RE = re.compile(r"(?:^|\|)UB:([^|]+)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Copy CB/UB values encoded in read names into BAM tags."
    )
    parser.add_argument("--input", required=True, help="Input BAM")
    parser.add_argument("--output", required=True, help="Output BAM")
    parser.add_argument(
        "--require-any-tag",
        action="store_true",
        help="Fail if no alignment read name contains both CB and UB values.",
    )
    parser.add_argument(
        "--require-all-tags",
        action="store_true",
        help="Fail unless every alignment read name contains both CB and UB values.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    total = 0
    tagged = 0
    with pysam.AlignmentFile(args.input, "rb") as in_bam:
        with pysam.AlignmentFile(args.output, "wb", template=in_bam) as out_bam:
            for read in in_bam:
                total += 1
                cb = CB_RE.search(read.query_name)
                ub = UB_RE.search(read.query_name)
                if cb:
                    read.set_tag("CB", cb.group(1), value_type="Z")
                if ub:
                    read.set_tag("UB", ub.group(1), value_type="Z")
                if cb and ub:
                    tagged += 1
                out_bam.write(read)

    print(f"[INFO] Alignments copied: {total}")
    print(f"[INFO] Alignments with CB and UB: {tagged}")
    if args.require_any_tag and tagged == 0:
        raise SystemExit("No alignments contained both CB and UB in the read name")
    if args.require_all_tags and (tagged == 0 or tagged != total):
        raise SystemExit(
            f"Expected CB/UB on every alignment, but tagged {tagged} of {total}"
        )


if __name__ == "__main__":
    main()
