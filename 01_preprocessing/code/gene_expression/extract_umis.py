#!/usr/bin/env python3
import argparse
import time

from umi_methods import SeqExperiment


def main():
    parser = argparse.ArgumentParser(
        description="Extract Smart-seq3 UMI reads from paired-end FASTQ files."
    )
    parser.add_argument("--fastq", "-f1", required=True, help="Read 1 FASTQ(.gz)")
    parser.add_argument("--fastq2", "-f2", required=True, help="Read 2 FASTQ(.gz)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--sample_name", required=True, help="Sample/cell barcode label")
    parser.add_argument("--tso", default="TTGCGCAATG", help="TSO sequence")
    parser.add_argument("--umi-length", type=int, default=8, help="UMI length after TSO")
    args = parser.parse_args()

    if not args.tso:
        parser.error("--tso must not be empty")
    if args.umi_length < 1:
        parser.error("--umi-length must be at least 1")

    start = time.time()
    exp = SeqExperiment(args.fastq, args.fastq2, tso=args.tso, umi_length=args.umi_length)
    umi_r1, umi_r2, non_umi_r1, non_umi_r2 = exp.write_umi_fastqs(
        args.output, args.sample_name
    )

    for path in (umi_r1, umi_r2, non_umi_r1, non_umi_r2):
        print(f"[INFO] Wrote: {path}")
    print(f"[INFO] Total read pairs: {exp.total_reads}")
    print(f"[INFO] Read pairs with an exact TSO match: {exp.tso_reads}")
    print(f"[INFO] UMI-positive read pairs: {exp.umi_reads}")
    print(f"[INFO] Non-UMI read pairs: {exp.non_umi_reads}")
    if exp.total_reads:
        print(f"[INFO] UMI-positive fraction: {exp.umi_reads / exp.total_reads:.6f}")
    print(f"[INFO] Elapsed seconds: {time.time() - start:.1f}")


if __name__ == "__main__":
    main()
