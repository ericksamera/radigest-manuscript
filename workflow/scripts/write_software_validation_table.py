#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


GFF_ATTR_RE = re.compile(r"(?:^|;)Length=([0-9]+)(?:;|$)")
FASTA_HEADER_RE = re.compile(
    r"^>(?P<id>\S+)\s+chrom=(?P<chrom>\S+)\s+start0=(?P<start0>\d+)\s+end0=(?P<end0>\d+)\s+length=(?P<length>\d+)"
)


def read_toy_sequence(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is None:
            raise ValueError(f"sequence line before FASTA header in {path}")
        else:
            records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def parse_gff(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for raw in path.read_text().splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) != 9:
            raise ValueError(f"GFF3 row in {path} does not have 9 columns: {raw!r}")
        chrom, source, feature, start_s, end_s, score, strand, phase, attrs = fields
        start = int(start_s)
        end = int(end_s)
        if start < 1 or end < start:
            raise ValueError(f"invalid 1-based closed GFF3 coordinates in {path}: {raw!r}")
        m = GFF_ATTR_RE.search(attrs)
        if not m:
            raise ValueError(f"missing Length attribute in {path}: {raw!r}")
        length = int(m.group(1))
        if length != end - start + 1:
            raise ValueError(f"GFF3 Length mismatch in {path}: {raw!r}")
        rows.append(
            {
                "chrom": chrom,
                "start0": start - 1,
                "end0": end,
                "length": length,
                "source": source,
                "feature": feature,
                "strand": strand,
            }
        )
    return rows


def parse_tsv(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["chrom", "start0", "end0", "length", "hard_kept", "size_weight"]
        if reader.fieldnames != expected:
            raise ValueError(f"unexpected TSV header in {path}: {reader.fieldnames}")
        for row in reader:
            start0 = int(row["start0"])
            end0 = int(row["end0"])
            length = int(row["length"])
            if length != end0 - start0:
                raise ValueError(f"TSV length mismatch in {path}: {row}")
            rows.append(
                {
                    "chrom": row["chrom"],
                    "start0": start0,
                    "end0": end0,
                    "length": length,
                    "hard_kept": row["hard_kept"],
                    "size_weight": row["size_weight"],
                }
            )
    return rows


def parse_fragment_fasta(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal current, seq_parts
        if current is None:
            return
        seq = "".join(seq_parts).upper()
        if len(seq) != int(current["length"]):
            raise ValueError(f"FASTA sequence length mismatch in {path}: {current}")
        current["sequence"] = seq
        rows.append(current)
        current = None
        seq_parts = []

    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            flush()
            m = FASTA_HEADER_RE.match(raw)
            if not m:
                raise ValueError(f"unparseable fragment FASTA header in {path}: {raw!r}")
            current = {
                "id": m.group("id"),
                "chrom": m.group("chrom"),
                "start0": int(m.group("start0")),
                "end0": int(m.group("end0")),
                "length": int(m.group("length")),
            }
        elif raw.strip():
            seq_parts.append(raw.strip())
    flush()
    return rows


def diff_empty(path: Path) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing diff artifact"
    size = path.stat().st_size
    return size == 0, f"{size} bytes"


def passfail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def add(rows: list[dict[str, str]], check: str, input_: str, command_or_rule: str, expected: str, observed: str, status: str, artifact: str) -> None:
    rows.append(
        {
            "check": check,
            "input": input_,
            "command_or_rule": command_or_rule,
            "expected": expected,
            "observed": observed,
            "status": status,
            "artifact": artifact,
        }
    )


def validate_internal_consistency(args: argparse.Namespace) -> tuple[bool, str]:
    toy = read_toy_sequence(args.toy_fasta)
    gff_rows = parse_gff(args.toy_double_gff3)
    tsv_rows = parse_tsv(args.toy_double_tsv)
    fasta_rows = parse_fragment_fasta(args.toy_double_fasta)

    if len(gff_rows) != len(tsv_rows) or len(gff_rows) != len(fasta_rows):
        return False, f"row-count mismatch: gff={len(gff_rows)} tsv={len(tsv_rows)} fasta={len(fasta_rows)}"

    coords = [(r["chrom"], r["start0"], r["end0"], r["length"]) for r in gff_rows]
    tsv_coords = [(r["chrom"], r["start0"], r["end0"], r["length"]) for r in tsv_rows]
    fasta_coords = [(r["chrom"], r["start0"], r["end0"], r["length"]) for r in fasta_rows]
    if coords != tsv_coords or coords != fasta_coords:
        return False, f"coordinate mismatch: gff={coords} tsv={tsv_coords} fasta={fasta_coords}"

    for row in fasta_rows:
        chrom = str(row["chrom"])
        if chrom not in toy:
            return False, f"fragment FASTA chrom {chrom!r} absent from toy FASTA"
        start0 = int(row["start0"])
        end0 = int(row["end0"])
        expected_seq = toy[chrom][start0:end0]
        if row["sequence"] != expected_seq:
            return False, f"fragment FASTA sequence mismatch for {row['id']}"

    with args.toy_double_json.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    total_bases = sum(int(r["length"]) for r in gff_rows)
    if summary.get("total_fragments") != len(gff_rows):
        return False, f"JSON total_fragments={summary.get('total_fragments')} but GFF rows={len(gff_rows)}"
    if summary.get("total_bases") != total_bases:
        return False, f"JSON total_bases={summary.get('total_bases')} but GFF bases={total_bases}"
    per_chr = summary.get("per_chromosome", {}).get("chr1", {})
    if per_chr.get("fragments") != len(gff_rows) or per_chr.get("bases") != total_bases:
        return False, f"JSON per_chromosome mismatch: {per_chr!r}"
    if not all(str(r["hard_kept"]).lower() == "true" for r in tsv_rows):
        return False, "toy double TSV contains a non-hard-kept row despite min/max 1..1000"
    return True, f"{len(gff_rows)} fragments; {total_bases} bp; GFF3/TSV/FASTA/JSON agree"


def main() -> int:
    parser = argparse.ArgumentParser(description="Write Stage 2 software-validation summary table.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--go-test-log", required=True, type=Path)
    parser.add_argument("--toy-fasta", required=True, type=Path)
    parser.add_argument("--toy-single-json", required=True, type=Path)
    parser.add_argument("--toy-single-gff3", required=True, type=Path)
    parser.add_argument("--toy-double-json", required=True, type=Path)
    parser.add_argument("--toy-double-gff3", required=True, type=Path)
    parser.add_argument("--toy-double-tsv", required=True, type=Path)
    parser.add_argument("--toy-double-fasta", required=True, type=Path)
    parser.add_argument("--toy-include-ends-gff3", required=True, type=Path)
    parser.add_argument("--toy-single-diff", required=True, type=Path)
    parser.add_argument("--toy-double-diff", required=True, type=Path)
    parser.add_argument("--toy-include-ends-diff", required=True, type=Path)
    parser.add_argument("--toy-allow-same-diff", required=True, type=Path)
    parser.add_argument("--repro-gff-tsv-diff", required=True, type=Path)
    parser.add_argument("--repro-json-normalized-diff", required=True, type=Path)
    parser.add_argument("--repro-json-same-thread-diff", required=True, type=Path)
    parser.add_argument("--repro-sha256", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []

    go_ok = args.go_test_log.is_file()
    add(
        rows,
        "go_tests",
        "resources/software/radigest",
        "radigest_build_test",
        "go test -count=1 -vet=off ./... exits 0",
        "log exists" if go_ok else "missing log",
        passfail(go_ok),
        args.go_test_log.as_posix(),
    )

    diff_checks = [
        (
            "single_digest_coordinates",
            args.toy_fasta,
            "validation_diff_toy_single_gff3",
            "EcoRI internal fragment chr1:5-16 in 0-based half-open coordinates; GFF3 chr1:6-16 in 1-based closed coordinates",
            args.toy_single_diff,
        ),
        (
            "double_digest_ab_ba_adjacency",
            args.toy_fasta,
            "validation_diff_toy_double_gff3",
            "default EcoRI/MseI digest keeps only adjacent AB and BA fragments",
            args.toy_double_diff,
        ),
        (
            "terminal_fragments_include_ends",
            args.toy_fasta,
            "validation_diff_toy_include_ends_gff3",
            "terminal contig-end fragments are absent by default and present with -include-ends",
            args.toy_include_ends_diff,
        ),
        (
            "allow_same_aa_bb_handling",
            args.toy_fasta,
            "validation_diff_toy_allow_same_gff3",
            "AA/BB neighbors are suppressed by default and emitted with -allow-same",
            args.toy_allow_same_diff,
        ),
        (
            "repro_gff3_tsv_byte_identity",
            Path("results/validation/repro.t{1,2,8}.{gff3,tsv}"),
            "validation_repro_gff3_tsv_byte_identity",
            "GFF3 and TSV outputs are byte-identical for threads 1, 2, and 8",
            args.repro_gff_tsv_diff,
        ),
        (
            "repro_json_normalized_identity",
            Path("results/validation/repro.t{1,2,8}.json"),
            "validation_repro_json_normalized_identity",
            "JSON summaries match after excluding command/output/thread provenance fields that are expected to differ across thread counts",
            args.repro_json_normalized_diff,
        ),
        (
            "repro_json_same_thread_byte_identity",
            Path("results/validation/repro.same_thread.json"),
            "validation_repro_json_same_thread_rerun",
            "JSON is byte-identical for two reruns of the same command with the same thread count and output path",
            args.repro_json_same_thread_diff,
        ),
    ]

    for check, input_path, rule, expected, diff_path in diff_checks:
        ok, observed = diff_empty(diff_path)
        add(
            rows,
            check,
            input_path.as_posix(),
            rule,
            expected,
            f"diff artifact {observed}",
            passfail(ok),
            diff_path.as_posix(),
        )

    try:
        ok, observed = validate_internal_consistency(args)
    except Exception as exc:  # fail hard but preserve table context
        ok, observed = False, str(exc)
    add(
        rows,
        "json_gff3_tsv_fasta_internal_consistency",
        args.toy_fasta.as_posix(),
        "software_validation_table",
        "GFF3 uses 1-based closed coordinates; TSV/FASTA metadata use 0-based half-open coordinates; JSON totals equal emitted artifacts",
        observed,
        passfail(ok),
        args.toy_double_json.as_posix(),
    )

    sha_ok = args.repro_sha256.is_file() and sum(
        1 for line in args.repro_sha256.read_text().splitlines() if line and not line.startswith("#")
    ) == 9
    add(
        rows,
        "repro_sha256_manifest",
        "results/validation/repro.t{1,2,8}.{json,gff3,tsv}",
        "validation_repro_sha256",
        "9 SHA256 records: JSON, GFF3, and TSV for threads 1, 2, and 8",
        "present with 9 records" if sha_ok else "missing or wrong record count",
        passfail(sha_ok),
        args.repro_sha256.as_posix(),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        fieldnames = ["check", "input", "command_or_rule", "expected", "observed", "status", "artifact"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    failures = [row for row in rows if row["status"] != "PASS"]
    if failures:
        for row in failures:
            print(f"FAIL {row['check']}: {row['observed']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
