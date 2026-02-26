#!/usr/bin/env python3
"""Extract 200 bp sequences around variants in a VCF from a reference FASTA."""

from __future__ import annotations

import argparse
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


def parse_vcf_records(path: str) -> Iterator[Tuple[str, int, str, str]]:
    """Yield (chrom, pos, ref, alt) for each ALT allele in the VCF."""
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
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
                yield chrom, pos, ref, alt


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
    records_written = 0

    with open(output_fasta_path, "w", encoding="utf-8") as out_handle:
        for chrom, pos, ref, alt in parse_vcf_records(vcf_path):
            if chrom not in genome:
                raise ValueError(f"Chromosome '{chrom}' not found in FASTA")

            start, end, seq = extract_window(genome[chrom], pos, window_size)
            header = f"{chrom}:{pos}:{ref}>{alt}|window={start}-{end}|len={len(seq)}"
            out_handle.write(f">{header}\n{seq}\n")
            records_written += 1

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
