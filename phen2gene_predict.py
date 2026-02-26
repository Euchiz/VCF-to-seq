#!/usr/bin/env python3
"""Rank candidate causal genes from patient HPO terms using Phen2Gene-style scoring.

This script implements a simple gene-centric ranking strategy inspired by Phen2Gene:
for each patient HPO term, collect all associated genes from an HPO annotation database,
then sum per-term weights into a final score for each gene.

Accepted annotation database formats
------------------------------------
1) Required columns in any order (header required):
   - hpo_id (or HP, HPO, hpo_term)
   - gene (or gene_symbol, symbol)
   - weight (or score) [optional; defaults to 1.0]

2) Delimiters: tab, comma, or whitespace (auto-detected).

Example:
    hpo_id\tgene\tweight
    HP:0001249\tMECP2\t5.43
    HP:0001511\tNSD1\t2.10
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

HPO_ALIASES = {"hpo_id", "hpo", "hp", "hpo_term", "term", "phenotype"}
GENE_ALIASES = {"gene", "gene_symbol", "symbol", "gene_name"}
WEIGHT_ALIASES = {"weight", "score", "beta", "term_gene_weight", "h2g_score"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Predict disease causal genes from patient HPO terms by summing HPO term–gene "
            "weights (Phen2Gene-style gene-centric ranking)."
        )
    )
    parser.add_argument("hpo_terms", help="Path to a file with one HPO term (HP:xxxxxxx) per line.")
    parser.add_argument(
        "annotation_db",
        help=(
            "Path to HPO term–gene annotation table with header. Must include HPO and gene "
            "columns; optional weight column."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default="predicted_genes.tsv",
        help="Output TSV path for ranked genes (default: predicted_genes.tsv).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=50,
        help="Number of top-ranked genes to print to stdout (default: 50).",
    )
    return parser.parse_args()


def normalize(text: str) -> str:
    return text.strip().lower()


def detect_delimiter(header_line: str) -> str:
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    return " "


def load_patient_terms(path: Path) -> List[str]:
    terms: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            terms.append(line)
    if not terms:
        raise ValueError(f"No HPO terms found in: {path}")
    return terms


def resolve_column_indices(fieldnames: Sequence[str]) -> Tuple[int, int, int | None]:
    hpo_idx = gene_idx = weight_idx = None
    for idx, name in enumerate(fieldnames):
        key = normalize(name)
        if key in HPO_ALIASES and hpo_idx is None:
            hpo_idx = idx
        elif key in GENE_ALIASES and gene_idx is None:
            gene_idx = idx
        elif key in WEIGHT_ALIASES and weight_idx is None:
            weight_idx = idx

    if hpo_idx is None or gene_idx is None:
        raise ValueError(
            "Annotation DB must include header columns for HPO term and gene "
            "(e.g., hpo_id, gene, weight)."
        )

    return hpo_idx, gene_idx, weight_idx


def parse_annotation_db(path: Path) -> Dict[str, Dict[str, float]]:
    term_gene_weights: Dict[str, Dict[str, float]] = defaultdict(dict)

    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"Annotation DB is empty: {path}")
        delimiter = detect_delimiter(first)
        handle.seek(0)

        if delimiter == " ":
            rows = [line.strip().split() for line in handle if line.strip()]
            header = rows[0]
            hpo_idx, gene_idx, weight_idx = resolve_column_indices(header)
            data_rows = rows[1:]
        else:
            reader = csv.reader(handle, delimiter=delimiter)
            rows = [row for row in reader if row]
            header = rows[0]
            hpo_idx, gene_idx, weight_idx = resolve_column_indices(header)
            data_rows = rows[1:]

    for row in data_rows:
        if len(row) <= max(hpo_idx, gene_idx):
            continue
        hpo = row[hpo_idx].strip()
        gene = row[gene_idx].strip()
        if not hpo or not gene:
            continue
        weight = 1.0
        if weight_idx is not None and weight_idx < len(row):
            raw_weight = row[weight_idx].strip()
            if raw_weight:
                try:
                    weight = float(raw_weight)
                except ValueError:
                    continue
        term_gene_weights[hpo][gene] = weight

    if not term_gene_weights:
        raise ValueError(f"No valid HPO-gene entries parsed from: {path}")

    return term_gene_weights


def score_genes(patient_terms: Iterable[str], term_gene_weights: Dict[str, Dict[str, float]]) -> Tuple[Dict[str, float], Dict[str, int], List[str]]:
    gene_scores: Dict[str, float] = defaultdict(float)
    matched_counts: Dict[str, int] = defaultdict(int)
    missing_terms: List[str] = []

    for term in patient_terms:
        associations = term_gene_weights.get(term)
        if not associations:
            missing_terms.append(term)
            continue

        for gene, weight in associations.items():
            gene_scores[gene] += weight
            matched_counts[gene] += 1

    return gene_scores, matched_counts, missing_terms


def write_output(path: Path, ranked: Sequence[Tuple[str, float, int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["rank", "gene", "score", "matched_hpo_terms"])
        for i, (gene, score, matched_terms) in enumerate(ranked, start=1):
            writer.writerow([i, gene, f"{score:.6f}", matched_terms])


def main() -> None:
    args = parse_args()

    hpo_terms = load_patient_terms(Path(args.hpo_terms))
    term_gene_weights = parse_annotation_db(Path(args.annotation_db))
    gene_scores, matched_counts, missing_terms = score_genes(hpo_terms, term_gene_weights)

    ranked = sorted(
        ((gene, score, matched_counts[gene]) for gene, score in gene_scores.items()),
        key=lambda item: (-item[1], -item[2], item[0]),
    )

    write_output(Path(args.output), ranked)

    print(f"Loaded {len(hpo_terms)} patient HPO terms")
    print(f"Matched {len(hpo_terms) - len(missing_terms)} terms; {len(missing_terms)} had no annotation")
    print(f"Ranked {len(ranked)} genes")
    print(f"Saved ranking to: {args.output}\n")

    if missing_terms:
        print("Unmatched terms:", ", ".join(missing_terms))
        print()

    top_n = max(0, args.top)
    if top_n == 0:
        return

    print(f"Top {min(top_n, len(ranked))} genes")
    print("rank\tgene\tscore\tmatched_hpo_terms")
    for i, (gene, score, matched_terms) in enumerate(ranked[:top_n], start=1):
        print(f"{i}\t{gene}\t{score:.6f}\t{matched_terms}")


if __name__ == "__main__":
    main()
