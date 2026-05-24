from pathlib import Path

configfile: "config/config.yaml"

RADIGEST_DIR = config["software"]["radigest"]["source_dir"]
RADIGEST = config["software"]["radigest"]["executable"]
RADIGEST_VERSION = str(config["software"]["radigest"]["version"])
RADIGEST_COMMIT = str(config["software"]["radigest"]["commit"])
RADIGEST_SOURCE_REVISION = config["software"]["radigest"]["source_revision_file"]

RADIGEST_HELPERS = [
    config["software"]["radigest"]["helper_screen_pairs"],
    config["software"]["radigest"]["helper_rank_pairs"],
    config["software"]["radigest"]["helper_fit_size_model"],
]

RADIGEST_SOURCE_TREE_SHA256 = "results/provenance/radigest.source_tree.sha256"
SOFTWARE_VALIDATION_TABLE = "results/tables/software_validation.tsv"

PROVENANCE = [
    config["provenance"]["go_test_log"],
    config["provenance"]["radigest_source_revision"],
    RADIGEST_SOURCE_TREE_SHA256,
    config["provenance"]["radigest_version"],
    config["provenance"]["supported_enzymes"],
    "results/provenance/candidate_enzymes.check.tsv",
    config["provenance"]["input_checksums"],
    config["provenance"]["workflow_config_snapshot"],
    config["provenance"]["software_versions"],
    config["provenance"]["directory_manifest"],
]

OPTIONAL_LOCK = [
    "conda-lock.yml",
    "results/provenance/conda_lock.log",
]


def expected_existing_inputs(wildcards):
    """Return existing files listed in metadata/expected_inputs.txt.

    Missing biological inputs are intentionally not returned as Snakemake inputs;
    the checksum script records them as MISSING instead.
    """
    paths = []
    manifest = Path("metadata/expected_inputs.txt")
    if not manifest.exists():
        return paths

    for raw in manifest.read_text().splitlines():
        rel = raw.strip()
        if not rel or rel.startswith("#"):
            continue
        if Path(rel).is_file():
            paths.append(rel)

    return paths


def radigest_source_files(wildcards):
    """Return vendored radigest source files, excluding build outputs and caches."""
    root = Path(RADIGEST_DIR)
    paths = []

    for path in root.rglob("*"):
        if not path.is_file():
            continue

        parts = path.relative_to(root).parts
        if ".git" in parts:
            continue
        if "bin" in parts:
            continue
        if "__pycache__" in parts:
            continue

        paths.append(path.as_posix())

    return sorted(paths)


rule all:
    input:
        PROVENANCE + OPTIONAL_LOCK + [SOFTWARE_VALIDATION_TABLE]


rule radigest_source_revision:
    input:
        RADIGEST_SOURCE_REVISION
    output:
        config["provenance"]["radigest_source_revision"]
    params:
        version=RADIGEST_VERSION,
        commit=RADIGEST_COMMIT
    shell:
        r"""
        set -euo pipefail
        {{
          date -u +"generated_utc=%Y-%m-%dT%H:%M:%SZ"
          printf 'configured_version\t%s\n' "{params.version}"
          printf 'configured_commit\t%s\n' "{params.commit}"
          printf 'source_revision_file\t%s\n' "{input}"
          echo
          cat "{input}"
        }} > "{output}"
        """


rule radigest_source_tree_sha256:
    input:
        source_files=radigest_source_files
    output:
        RADIGEST_SOURCE_TREE_SHA256
    shell:
        r"""
        set -euo pipefail
        {{
          date -u +"# generated_utc	%Y-%m-%dT%H:%M:%SZ"
          printf '# format\tsha256  project_relative_path\n'
          for f in {input.source_files}; do
            sha256sum "$f"
          done | LC_ALL=C sort
        }} > "{output}"
        """


rule radigest_build_test:
    input:
        source_files=radigest_source_files
    output:
        log=config["provenance"]["go_test_log"],
        radigest=RADIGEST,
        screen=config["software"]["radigest"]["helper_screen_pairs"],
        rank=config["software"]["radigest"]["helper_rank_pairs"],
        fit=config["software"]["radigest"]["helper_fit_size_model"]
    params:
        root=RADIGEST_DIR,
        version=RADIGEST_VERSION,
        commit=RADIGEST_COMMIT
    shell:
        r"""
        set -euo pipefail
        OUT="$(pwd)/{output.log}"

        cd "{params.root}"

        {{
          date -u +"generated_utc=%Y-%m-%dT%H:%M:%SZ"
          echo "snapshot_dir=$(pwd)"
          echo "vendored_source_version={params.version}"
          echo "vendored_source_commit={params.commit}"
          echo "conda_prefix=${{CONDA_PREFIX:-NA}}"
          echo "python=$(command -v python3)"
          echo "go=$(command -v go)"
          echo "snakemake=$(command -v snakemake)"

          echo
          echo "[tool_versions]"
          python3 --version
          go version
          snakemake --version

          echo
          echo "[go_env]"
          go env

          echo
          echo "[go_generate]"
          TMP_GEN="$(mktemp)"
          trap 'rm -f "$TMP_GEN"' EXIT
          cp -p internal/enzyme/enzymes_generated.go "$TMP_GEN"
          go generate ./...
          if ! cmp -s "$TMP_GEN" internal/enzyme/enzymes_generated.go; then
            echo "go generate changed internal/enzyme/enzymes_generated.go" >&2
            diff -u "$TMP_GEN" internal/enzyme/enzymes_generated.go || true
            exit 1
          fi
          touch -r "$TMP_GEN" internal/enzyme/enzymes_generated.go
          echo "go_generate_clean=OK"

          echo
          echo "[make_clean]"
          make clean

          echo
          echo "[make_build]"
          make VERSION="{params.version}" build

          echo
          echo "[go_test]"
          go test -count=1 -vet=off ./...
        }} 2>&1 | tee "$OUT"
        """


rule radigest_version:
    input:
        radigest=RADIGEST
    output:
        config["provenance"]["radigest_version"]
    shell:
        "{input.radigest} -version > {output}"


rule supported_enzymes:
    input:
        radigest=RADIGEST
    output:
        config["provenance"]["supported_enzymes"]
    shell:
        "{input.radigest} -list-enzymes > {output}"


rule candidate_enzymes_check:
    input:
        script="workflow/scripts/check_candidate_enzymes.py",
        candidate=config["inputs"]["candidate_enzymes"],
        supported=config["provenance"]["supported_enzymes"]
    output:
        "results/provenance/candidate_enzymes.check.tsv"
    shell:
        "python3 {input.script} --candidate {input.candidate} --supported {input.supported} --out {output}"


rule input_checksums:
    input:
        expected=config["inputs"]["expected_inputs"],
        script="workflow/scripts/write_input_checksums.py",
        tracked=expected_existing_inputs
    output:
        config["provenance"]["input_checksums"]
    shell:
        "python3 {input.script} --inputs {input.expected} --out {output}"


rule workflow_config_snapshot:
    input:
        "config/config.yaml"
    output:
        config["provenance"]["workflow_config_snapshot"]
    shell:
        "cp {input} {output}"


rule software_versions:
    input:
        script="workflow/scripts/write_software_versions.py",
        radigest=RADIGEST,
        screen=config["software"]["radigest"]["helper_screen_pairs"],
        rank=config["software"]["radigest"]["helper_rank_pairs"],
        fit=config["software"]["radigest"]["helper_fit_size_model"]
    output:
        config["provenance"]["software_versions"]
    shell:
        "python3 {input.script} --radigest {input.radigest} --go-command go --snakemake-command snakemake --out {output}"


rule directory_manifest:
    output:
        config["provenance"]["directory_manifest"]
    shell:
        r"""
        set -euo pipefail
        {{
          date -u +"# generated_utc	%Y-%m-%dT%H:%M:%SZ"
          printf '# root\t%s\n' "$PWD"
          printf '# format\tproject_relative_directory\n'
          find config envs metadata workflow resources results \
            -path '*/.git' -prune -o -type d -print | LC_ALL=C sort
        }} > "{output}"
        """


rule conda_lock:
    input:
        "envs/environment.yml"
    output:
        lockfile="conda-lock.yml",
        log="results/provenance/conda_lock.log"
    shell:
        "conda-lock -f {input} -p linux-64 --lockfile {output.lockfile} > {output.log} 2>&1"


include: "workflow/rules/validation.smk"
