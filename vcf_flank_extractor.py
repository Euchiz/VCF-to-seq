#!/usr/bin/env python3
"""Extract 200 bp sequences around variants in a VCF from a reference FASTA."""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Iterator, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read variants from a VCF, extract sequence windows from a genome FASTA, "
            "and write them to an output FASTA."
        )
    )
    parser.add_argument("vcf", help="Input VCF file.")
    parser.add_argument("genome_fasta", help="Reference genome FASTA file.")
    parser.add_argument("output_fasta", help="Output FASTA file.")
    parser.add_argument(
        "--window-size",
        type=int,
        default=200,
        help="Total length of sequence to extract around each variant (default: 200).",
    )
    return parser.parse_args()


def read_fasta(path: str) -> Dict[str, str]:
    """Load a FASTA file into memory as {contig: sequence}."""
    sequences: Dict[str, list[str]] = {}
    current_name: str | None = None

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                if current_name in sequences:
                    raise ValueError(f"Duplicate FASTA contig name: {current_name}")
                sequences[current_name] = []
            else:
                if current_name is None:
                    raise ValueError("FASTA sequence data encountered before any header line")
                sequences[current_name].append(line.upper())

    return {name: "".join(chunks) for name, chunks in sequences.items()}


def build_contig_lookup(genome: Dict[str, str]) -> Dict[str, str]:
    """Build a lookup that maps common chromosome aliases to FASTA contig names."""
    lookup: Dict[str, str] = {}

    for contig in genome:
        candidates = {contig}

        if contig.startswith("chr"):
            candidates.add(contig[3:])
        else:
            candidates.add(f"chr{contig}")

        if contig in {"MT", "chrMT"}:
            candidates.update({"M", "chrM"})
        if contig in {"M", "chrM"}:
            candidates.update({"MT", "chrMT"})

        for candidate in candidates:
            if candidate in lookup and lookup[candidate] != contig:
                raise ValueError(
                    "Ambiguous contig alias mapping for "
                    f"'{candidate}' ({lookup[candidate]} vs {contig})"
                )
            lookup[candidate] = contig

    return lookup


def parse_vcf_records(path: str) -> Iterator[Tuple[int, str, int, str, str]]:
    """Yield (line_no, chrom, pos, ref, alt) for each ALT allele in the VCF."""
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Malformed VCF line (expected >=5 fields): {line.rstrip()}")

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")

            for alt in alts:
                if alt == ".":
                    continue
                yield line_no, chrom, pos, ref, alt


def is_valid_allele(allele: str) -> bool:
    """Return True when the allele uses standard DNA/IUPAC N bases."""
    return bool(allele) and set(allele.upper()) <= {"A", "C", "G", "T", "N"}


def reference_matches(sequence: str, pos_1_based: int, ref: str) -> bool:
    """Check whether REF matches the reference genome at POS."""
    start = pos_1_based - 1
    end = start + len(ref)
    if start < 0 or end > len(sequence):
        return False
    return sequence[start:end].upper() == ref.upper()


def extract_window(sequence: str, pos_1_based: int, window_size: int) -> Tuple[int, int, str]:
    """Return (start_1_based, end_1_based, subsequence) around the variant position."""
    if window_size <= 0:
        raise ValueError("window_size must be a positive integer")

    half = window_size // 2
    start = pos_1_based - half
    end = start + window_size - 1

    if start < 1:
        start = 1
        end = min(window_size, len(sequence))
    if end > len(sequence):
        end = len(sequence)
        start = max(1, end - window_size + 1)

    subseq = sequence[start - 1 : end]
    return start, end, subseq


def write_variant_windows(
    vcf_path: str,
    genome_fasta_path: str,
    output_fasta_path: str,
    window_size: int,
) -> int:
    genome = read_fasta(genome_fasta_path)
    contig_lookup = build_contig_lookup(genome)
    records_written = 0

    skipped_records = 0

    with open(output_fasta_path, "w", encoding="utf-8") as out_handle:
        for line_no, chrom, pos, ref, alt in parse_vcf_records(vcf_path):
            resolved_chrom = contig_lookup.get(chrom)
            if resolved_chrom is None:
                raise ValueError(f"Chromosome '{chrom}' not found in FASTA")

            if not is_valid_allele(ref) or not is_valid_allele(alt):
                print(
                    f"Skipping VCF line {line_no}: invalid REF/ALT bases ({ref}>{alt})",
                    file=sys.stderr,
                )
                skipped_records += 1
                continue

            if not reference_matches(genome[resolved_chrom], pos, ref):
                print(
                    f"Skipping VCF line {line_no}: REF '{ref}' does not match "
                    f"{resolved_chrom}:{pos} in FASTA",
                    file=sys.stderr,
                )
                skipped_records += 1
                continue

            start, end, seq = extract_window(genome[resolved_chrom], pos, window_size)
            if seq and set(seq) == {"N"}:
                print(
                    f"Skipping VCF line {line_no}: extracted window is all Ns "
                    f"({resolved_chrom}:{start}-{end})",
                    file=sys.stderr,
                )
                skipped_records += 1
                continue

            header = f"{resolved_chrom}:{pos}:{ref}>{alt}|window={start}-{end}|len={len(seq)}"
            out_handle.write(f">{header}\n{seq}\n")
            records_written += 1

    if skipped_records:
        print(f"Skipped {skipped_records} problematic VCF record(s)", file=sys.stderr)

    return records_written


def main() -> None:
    args = parse_args()
    total = write_variant_windows(
        args.vcf,
        args.genome_fasta,
        args.output_fasta,
        args.window_size,
    )
    print(f"Wrote {total} variant sequence(s) to {args.output_fasta}")


if __name__ == "__main__":
    main()
