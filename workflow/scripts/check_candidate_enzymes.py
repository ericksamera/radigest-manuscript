#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_names(path: Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text().splitlines():
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        names.append(value)
    return names


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check candidate enzyme names against a radigest supported-enzyme snapshot."
    )
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--supported", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    candidates = read_names(args.candidate)
    supported = set(read_names(args.supported))

    args.out.parent.mkdir(parents=True, exist_ok=True)

    unsupported: list[str] = []
    duplicated: set[str] = set()
    seen: set[str] = set()

    with args.out.open("w") as out:
        out.write("enzyme\tstatus\n")
        for enzyme in candidates:
            if enzyme in seen:
                duplicated.add(enzyme)
            seen.add(enzyme)

            status = "supported" if enzyme in supported else "unsupported"
            out.write(f"{enzyme}\t{status}\n")

            if status == "unsupported":
                unsupported.append(enzyme)

    if duplicated:
        print(
            "warning: duplicated candidate enzyme(s): "
            + ", ".join(sorted(duplicated)),
            file=sys.stderr,
        )

    if unsupported:
        print(
            "unsupported candidate enzyme(s): "
            + ", ".join(unsupported),
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
