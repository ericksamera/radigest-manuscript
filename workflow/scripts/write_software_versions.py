#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import platform
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def command_exists(cmd: list[str]) -> bool:
    if not cmd:
        return False
    first = cmd[0]
    return Path(first).exists() or shutil.which(first) is not None


def capture(command: str, timeout: int = 60) -> tuple[str, str, str]:
    cmd = shlex.split(command)
    if not command_exists(cmd):
        return ("MISSING", "NA", command)

    try:
        proc = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except Exception as exc:
        return ("ERROR", str(exc).replace("\n", "\\n"), command)

    output = proc.stdout.strip().replace("\n", "\\n")
    if not output:
        output = "<no output>"

    return (f"exit={proc.returncode}", output, command)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Write software version and binary checksum provenance."
    )
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--radigest", required=True, type=Path)
    parser.add_argument("--go-command", default="go")
    parser.add_argument("--snakemake-command", default="snakemake")
    args = parser.parse_args()

    rows: list[tuple[str, str, str, str]] = []

    commands = {
        "python_current": f"{shlex.quote(sys.executable)} --version",
        "mamba": "mamba --version",
        "conda": "conda --version",
        "git": "git --version",
        "make": "make --version",
        "go": f"{args.go_command} version",
        "snakemake": f"{args.snakemake_command} --version",
        "radigest": f"{shlex.quote(str(args.radigest))} -version",
        "samtools": "samtools --version",
        "bwa_mem2": "bwa-mem2 version",
        "fastp": "fastp --version",
        "bbmap": "bbversion.sh",
        "hyperfine": "hyperfine --version",
        "jq": "jq --version",
    }

    for name, command in commands.items():
        status, output, rendered = capture(command)
        rows.append((name, status, output, rendered))

    binary_paths = [
        args.radigest,
        args.radigest.parent / "radigest-screen-pairs",
        args.radigest.parent / "radigest-rank-pairs",
        args.radigest.parent / "radigest-fit-size-model",
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open("w") as out:
        out.write("# generated_utc\t" + dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ") + "\n")
        out.write("# platform\t" + platform.platform() + "\n")
        out.write("# python_executable\t" + sys.executable + "\n")
        out.write("record_type\tname\tstatus_or_sha256\tvalue\tcommand_or_path\n")

        for name, status, output, command in rows:
            out.write(f"command\t{name}\t{status}\t{output}\t{command}\n")

        for path in binary_paths:
            if path.is_file():
                out.write(f"sha256\t{path.name}\t{sha256_file(path)}\tNA\t{path}\n")
            else:
                out.write(f"sha256\t{path.name}\tMISSING\tNA\t{path}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
