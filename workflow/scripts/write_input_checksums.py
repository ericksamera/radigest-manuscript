#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write SHA256 checksums for expected immutable workflow inputs."
    )
    parser.add_argument(
        "--inputs",
        required=True,
        type=Path,
        help="Text file containing one project-relative input path per line.",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output checksum manifest.",
    )
    args = parser.parse_args()

    root = Path.cwd()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# generated_utc\t" + dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "# root\t" + str(root),
        "# format\tsha256_or_status  relative_path",
    ]

    for raw in args.inputs.read_text().splitlines():
        rel = raw.strip()
        if not rel or rel.startswith("#"):
            continue

        path = root / rel
        if path.is_file():
            lines.append(f"{sha256_file(path)}  {rel}")
        elif path.exists():
            lines.append(f"NONFILE  {rel}")
        else:
            lines.append(f"MISSING  {rel}")

    args.out.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
