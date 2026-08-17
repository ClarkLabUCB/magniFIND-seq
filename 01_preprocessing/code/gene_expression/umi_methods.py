import gzip
import io
import os
from contextlib import contextmanager

from xopen import xopen


@contextmanager
def deterministic_gzip_text(path):
    """Write gzip text with no filename or wall-clock timestamp in its header."""
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as text:
                yield text


class SeqExperiment:
    """Reproduce the production Smart-seq3 UMI/non-UMI FASTQ split."""

    def __init__(self, file1, file2, tso="TTGCGCAATG", umi_length=8):
        self.file1 = file1
        self.file2 = file2
        self.tso = tso
        self.umi_length = int(umi_length)
        if not self.tso:
            raise ValueError("TSO sequence must not be empty")
        if self.umi_length < 1:
            raise ValueError("UMI length must be at least 1")
        self.total_reads = 0
        self.tso_reads = 0
        self.umi_reads = 0
        self.non_umi_reads = 0

    @staticmethod
    def _read_fastq_record(handle):
        header = handle.readline()
        if not header:
            return None
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("Incomplete FASTQ record")
        if not header.startswith("@"):
            raise ValueError(f"FASTQ header does not start with '@': {header.rstrip()}")
        if not plus.startswith("+"):
            raise ValueError(f"FASTQ plus line does not start with '+': {plus.rstrip()}")
        read_id = header.rstrip()[1:]
        sequence = seq.rstrip()
        quality = qual.rstrip()
        if len(sequence) != len(quality):
            raise ValueError(f"Sequence/quality length mismatch for read: {read_id}")
        return read_id, sequence, quality

    @staticmethod
    def _normalize_read_id(read_id):
        token = read_id.split()[0]
        if token.endswith(("/1", "/2")):
            token = token[:-2]
        return token

    def iter_reads(self):
        with xopen(self.file1, "rt") as f1, xopen(self.file2, "rt") as f2:
            while True:
                rec1 = self._read_fastq_record(f1)
                rec2 = self._read_fastq_record(f2)
                if rec1 is None and rec2 is None:
                    return
                if rec1 is None or rec2 is None:
                    raise ValueError("FASTQ pair has different number of records")

                id1, seq1, qual1 = rec1
                id2, seq2, qual2 = rec2
                read_id = self._normalize_read_id(id1)
                if read_id != self._normalize_read_id(id2):
                    raise ValueError(f"Read IDs do not match: {id1} vs {id2}")

                self.total_reads += 1
                tso_pos = seq1.find(self.tso)
                umi = ""
                if tso_pos != -1:
                    self.tso_reads += 1
                    after_tso = tso_pos + len(self.tso)
                    umi = seq1[after_tso : after_tso + self.umi_length]
                    trimmed_start = after_tso + self.umi_length
                    seq1 = seq1[trimmed_start:]
                    qual1 = qual1[trimmed_start:]

                # The production implementation classified any non-empty UMI
                # substring as UMI-positive. This intentionally retains a UMI
                # shorter than umi_length when a read terminates just after TSO,
                # and it retains an empty post-UMI R1 sequence.
                if umi:
                    self.umi_reads += 1
                else:
                    self.non_umi_reads += 1
                yield read_id, umi, seq1, qual1, seq2, qual2

    def iter_umi_reads(self):
        """Yield only the UMI-positive records from the production split."""
        for record in self.iter_reads():
            if record[1]:
                yield record

    def write_umi_fastqs(self, out_dir, sample_name):
        os.makedirs(out_dir, exist_ok=True)

        umi_r1_out = os.path.join(out_dir, f"{sample_name}_R1_umi_out.fastq.gz")
        umi_r2_out = os.path.join(out_dir, f"{sample_name}_R2_umi_out.fastq.gz")
        non_umi_r1_out = os.path.join(out_dir, f"{sample_name}_R1_out.fastq.gz")
        non_umi_r2_out = os.path.join(out_dir, f"{sample_name}_R2_out.fastq.gz")

        with (
            deterministic_gzip_text(umi_r1_out) as umi_out1,
            deterministic_gzip_text(umi_r2_out) as umi_out2,
            deterministic_gzip_text(non_umi_r1_out) as non_umi_out1,
            deterministic_gzip_text(non_umi_r2_out) as non_umi_out2,
        ):
            for read_id, umi, seq1, qual1, seq2, qual2 in self.iter_reads():
                annotated = f"{read_id}|CB:{sample_name}|UB:{umi}"
                out1, out2 = (
                    (umi_out1, umi_out2) if umi else (non_umi_out1, non_umi_out2)
                )
                out1.write(f"@{annotated}\n{seq1}\n+\n{qual1}\n")
                out2.write(f"@{annotated}\n{seq2}\n+\n{qual2}\n")

        return umi_r1_out, umi_r2_out, non_umi_r1_out, non_umi_r2_out
