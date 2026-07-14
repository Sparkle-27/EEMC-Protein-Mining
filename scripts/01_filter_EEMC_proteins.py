#!/usr/bin/env python3
"""Filter and normalize Prodigal protein FASTA files for the EEMC workflow.

Filtering rules
---------------
1. Protein length must be within the inclusive user-defined range.
2. By default, the Prodigal header must contain partial=00.
3. Sequence characters must be valid protein letters.
4. If an external completeness table is supplied, score must be >= threshold.

The program prefixes each original protein ID with its genome ID:
    GenomeID|OriginalProteinID

GenomeID is derived from the input FAA filename.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, TextIO, Tuple


PARTIAL_RE = re.compile(r"(?:^|[;\s])partial=([01]{2})(?:;|\s|$)")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter EMCC Prodigal proteins and create a normalized FASTA."
    )
    parser.add_argument("--file-list", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--kept-table", required=True, type=Path)
    parser.add_argument("--discarded-table", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--min-length", type=int, default=40)
    parser.add_argument("--max-length", type=int, default=1000)
    parser.add_argument(
        "--require-partial00",
        action="store_true",
        help="Keep only proteins whose Prodigal header contains partial=00.",
    )
    parser.add_argument(
        "--completeness-table",
        type=Path,
        default=None,
        help=(
            "Optional TSV with columns protein_id and completeness. "
            "protein_id must use GenomeID|OriginalProteinID."
        ),
    )
    parser.add_argument("--min-completeness", type=float, default=90.0)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def genome_id_from_path(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in (".faa", ".fasta", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def fasta_records(handle: TextIO) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    seq_parts = []

    for raw in handle:
        line = raw.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line[1:]
            seq_parts = []
        else:
            if header is None:
                raise ValueError("FASTA sequence appeared before the first header.")
            seq_parts.append("".join(line.split()))

    if header is not None:
        yield header, "".join(seq_parts)


def read_completeness(path: Optional[Path]) -> Dict[str, float]:
    if path is None:
        return {}

    values: Dict[str, float] = {}
    with path.open("rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        lowered = [x.strip().lower() for x in header]

        try:
            id_idx = lowered.index("protein_id")
            score_idx = lowered.index("completeness")
        except ValueError as exc:
            raise ValueError(
                "Completeness TSV must contain columns: protein_id and completeness"
            ) from exc

        for line_no, raw in enumerate(handle, start=2):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            try:
                protein_id = fields[id_idx]
                score = float(fields[score_idx])
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"Invalid completeness row at line {line_no}: {raw.rstrip()}"
                ) from exc
            values[protein_id] = score

    return values


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    for start in range(0, len(seq), width):
        yield seq[start : start + width]


def main() -> int:
    args = parse_args()

    if args.min_length < 1 or args.max_length < args.min_length:
        raise ValueError("Invalid protein length range.")
    if not args.file_list.is_file():
        raise FileNotFoundError(args.file_list)

    completeness = read_completeness(args.completeness_table)
    use_external_completeness = args.completeness_table is not None

    input_paths = [
        Path(line.strip())
        for line in args.file_list.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not input_paths:
        raise ValueError("Input FAA file list is empty.")

    for output in (
        args.output_fasta,
        args.kept_table,
        args.discarded_table,
        args.summary,
    ):
        output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen_ids = set()

    with args.output_fasta.open("wt", encoding="utf-8") as fasta_out, \
         args.kept_table.open("wt", encoding="utf-8") as kept_out, \
         args.discarded_table.open("wt", encoding="utf-8") as discarded_out:
        kept_out.write(
            "protein_id\tgenome_id\toriginal_id\tlength\tpartial\t"
            "external_completeness\tcontains_X\tsource_faa\n"
        )
        discarded_out.write(
            "protein_id\tgenome_id\toriginal_id\tlength\tpartial\t"
            "external_completeness\treasons\tsource_faa\n"
        )

        for faa_path in input_paths:
            if not faa_path.is_file():
                raise FileNotFoundError(faa_path)

            genome_id = genome_id_from_path(faa_path)
            counts["input_files"] += 1

            with open_text(faa_path) as handle:
                for header, raw_seq in fasta_records(handle):
                    counts["input_proteins"] += 1

                    original_id = header.split(None, 1)[0]
                    protein_id = f"{genome_id}|{original_id}"

                    if protein_id in seen_ids:
                        raise ValueError(f"Duplicated normalized protein ID: {protein_id}")
                    seen_ids.add(protein_id)

                    seq = raw_seq.upper().replace(" ", "").replace("\t", "")
                    seq = seq.rstrip("*")
                    length = len(seq)

                    partial_match = PARTIAL_RE.search(header)
                    partial = partial_match.group(1) if partial_match else "NA"

                    external_score = completeness.get(protein_id)
                    score_text = (
                        f"{external_score:g}" if external_score is not None else "NA"
                    )

                    reasons = []
                    if length < args.min_length:
                        reasons.append("length_below_min")
                    if length > args.max_length:
                        reasons.append("length_above_max")

                    if args.require_partial00:
                        if partial == "NA":
                            reasons.append("missing_partial_tag")
                        elif partial != "00":
                            reasons.append(f"partial_{partial}")

                    illegal = sorted(set(seq) - VALID_AA)
                    if illegal:
                        reasons.append("illegal_residue:" + ",".join(illegal))
                    if not seq:
                        reasons.append("empty_sequence")

                    if use_external_completeness:
                        if external_score is None:
                            reasons.append("missing_external_completeness")
                        elif external_score < args.min_completeness:
                            reasons.append("external_completeness_below_threshold")

                    if reasons:
                        counts["discarded_proteins"] += 1
                        for reason in reasons:
                            counts["discard_" + reason.split(":", 1)[0]] += 1
                        discarded_out.write(
                            f"{protein_id}\t{genome_id}\t{original_id}\t{length}\t"
                            f"{partial}\t{score_text}\t{';'.join(reasons)}\t"
                            f"{faa_path}\n"
                        )
                        continue

                    counts["retained_proteins"] += 1
                    if "X" in seq:
                        counts["retained_with_X"] += 1

                    remainder = header[len(original_id) :]
                    fasta_out.write(f">{protein_id}{remainder}\n")
                    for chunk in wrap_fasta(seq):
                        fasta_out.write(chunk + "\n")

                    kept_out.write(
                        f"{protein_id}\t{genome_id}\t{original_id}\t{length}\t"
                        f"{partial}\t{score_text}\t{int('X' in seq)}\t"
                        f"{faa_path}\n"
                    )
    summary_order = [
        "input_files",
        "input_proteins",
        "retained_proteins",
        "discarded_proteins",
        "discard_length_below_min",
        "discard_length_above_max",
        "discard_missing_partial_tag",
        "discard_partial_01",
        "discard_partial_10",
        "discard_partial_11",
        "discard_illegal_residue",
        "discard_empty_sequence",
        "discard_missing_external_completeness",
        "discard_external_completeness_below_threshold",
        "retained_with_X",
    ]

    with args.summary.open("wt", encoding="utf-8") as summary_out:
        summary_out.write("metric\tvalue\n")
        for metric in summary_order:
            summary_out.write(f"{metric}\t{counts.get(metric, 0)}\n")
        summary_out.write(
            f"external_completeness_filter_applied\t"
            f"{int(use_external_completeness)}\n"
        )
        summary_out.write(
            f"external_completeness_threshold\t"
            f"{args.min_completeness if use_external_completeness else 'NA'}\n"
        )

    print(f"Input FAA files:     {counts['input_files']}")
    print(f"Input proteins:      {counts['input_proteins']}")
    print(f"Retained proteins:   {counts['retained_proteins']}")
    print(f"Discarded proteins:  {counts['discarded_proteins']}")
    print(f"Filtered FASTA:      {args.output_fasta}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
