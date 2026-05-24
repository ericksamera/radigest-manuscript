# Stage 2 — digest correctness and determinism validation.
# This file is deliberately limited to software validation. It does not perform
# enzyme screening, empirical TLEN fitting, transfer analysis, benchmarking, or
# manuscript figure generation.

VALIDATION = config.get("validation", {})
VALIDATION_DIR = "results/validation"
VALIDATION_CHECK_DIR = f"{VALIDATION_DIR}/checks"
VALIDATION_TOY_FASTA = VALIDATION.get("toy_fasta", "resources/validation/toy.fa")
VALIDATION_EXPECTED_TOY_DOUBLE_GFF3 = VALIDATION.get(
    "expected_toy_double_gff3", "resources/validation/expected.toy.double.gff3"
)
VALIDATION_REPRO_FASTA = VALIDATION.get(
    "reproducibility_reference_fasta", config["systems"]["salmon"]["reference_fasta"]
)
VALIDATION_REPRO_ENZYMES = VALIDATION.get("reproducibility_enzymes", "PstI,MspI")
VALIDATION_HARD_MIN = int(VALIDATION.get("hard_size_window_bp", {}).get("min", 300))
VALIDATION_HARD_MAX = int(VALIDATION.get("hard_size_window_bp", {}).get("max", 600))
VALIDATION_SCORE_MIN = int(VALIDATION.get("score_range_bp", {}).get("min", 1))
VALIDATION_SCORE_MAX = int(VALIDATION.get("score_range_bp", {}).get("max", 2000))
VALIDATION_THREADS = [1, 2, 8]

TOY_SINGLE_JSON = f"{VALIDATION_DIR}/toy.single.json"
TOY_SINGLE_GFF3 = f"{VALIDATION_DIR}/toy.single.gff3"
TOY_DOUBLE_JSON = f"{VALIDATION_DIR}/toy.double.json"
TOY_DOUBLE_GFF3 = f"{VALIDATION_DIR}/toy.double.gff3"
TOY_DOUBLE_TSV = f"{VALIDATION_DIR}/toy.double.tsv"
TOY_DOUBLE_FASTA = f"{VALIDATION_DIR}/toy.double.fa"
TOY_INCLUDE_ENDS_GFF3 = f"{VALIDATION_DIR}/toy.include_ends.gff3"

VALIDATION_DIFFS = [
    f"{VALIDATION_CHECK_DIR}/toy.single.gff3.diff",
    f"{VALIDATION_CHECK_DIR}/toy.double.gff3.diff",
    f"{VALIDATION_CHECK_DIR}/toy.include_ends.gff3.diff",
    f"{VALIDATION_CHECK_DIR}/toy.allow_same.gff3.diff",
    f"{VALIDATION_CHECK_DIR}/repro.gff3_tsv.byte_identity.diff",
    f"{VALIDATION_CHECK_DIR}/repro.json.normalized_identity.diff",
    f"{VALIDATION_CHECK_DIR}/repro.json.same_thread_rerun.diff",
]

REPRO_JSON = expand(f"{VALIDATION_DIR}/repro.t{{thread}}.json", thread=VALIDATION_THREADS)
REPRO_GFF3 = expand(f"{VALIDATION_DIR}/repro.t{{thread}}.gff3", thread=VALIDATION_THREADS)
REPRO_TSV = expand(f"{VALIDATION_DIR}/repro.t{{thread}}.tsv", thread=VALIDATION_THREADS)


rule validation_toy_single_digest:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_TOY_FASTA
    output:
        json=TOY_SINGLE_JSON,
        gff=TOY_SINGLE_GFF3
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.json})"
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes EcoRI \
          -min 1 -max 1000 \
          -json "{output.json}" \
          -gff "{output.gff}" \
          -threads 1
        """


rule validation_toy_double_digest:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_TOY_FASTA
    output:
        json=TOY_DOUBLE_JSON,
        gff=TOY_DOUBLE_GFF3,
        tsv=TOY_DOUBLE_TSV,
        fasta=TOY_DOUBLE_FASTA
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.json})"
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes EcoRI,MseI \
          -min 1 -max 1000 \
          -json "{output.json}" \
          -gff "{output.gff}" \
          -fragments-tsv "{output.tsv}" \
          -fragments-fasta "{output.fasta}" \
          -threads 1
        """


rule validation_toy_include_ends_digest:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_TOY_FASTA
    output:
        gff=TOY_INCLUDE_ENDS_GFF3
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.gff})"
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes EcoRI,MseI \
          -include-ends \
          -min 1 -max 1000 \
          -gff "{output.gff}" \
          -threads 1
        """


rule validation_diff_toy_single_gff3:
    input:
        observed=TOY_SINGLE_GFF3
    output:
        f"{VALIDATION_CHECK_DIR}/toy.single.gff3.diff"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        expected="$(mktemp)"
        trap 'rm -f "$expected"' EXIT
        cat > "$expected" <<'EOF'
##gff-version 3
chr1	radigest	fragment	6	16	.	+	.	ID=chr1_1;Length=11
EOF
        diff -u "$expected" "{input.observed}" > "{output}"
        """


rule validation_diff_toy_double_gff3:
    input:
        expected=VALIDATION_EXPECTED_TOY_DOUBLE_GFF3,
        observed=TOY_DOUBLE_GFF3
    output:
        f"{VALIDATION_CHECK_DIR}/toy.double.gff3.diff"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        diff -u "{input.expected}" "{input.observed}" > "{output}"
        """


rule validation_diff_toy_include_ends_gff3:
    input:
        observed=TOY_INCLUDE_ENDS_GFF3
    output:
        f"{VALIDATION_CHECK_DIR}/toy.include_ends.gff3.diff"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        expected="$(mktemp)"
        trap 'rm -f "$expected"' EXIT
        cat > "$expected" <<'EOF'
##gff-version 3
chr1	radigest	fragment	1	5	.	+	.	ID=chr1_1;Length=5
chr1	radigest	fragment	6	11	.	+	.	ID=chr1_2;Length=6
chr1	radigest	fragment	12	16	.	+	.	ID=chr1_3;Length=5
chr1	radigest	fragment	17	21	.	+	.	ID=chr1_4;Length=5
EOF
        diff -u "$expected" "{input.observed}" > "{output}"
        """


rule validation_diff_toy_allow_same_gff3:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_TOY_FASTA
    output:
        f"{VALIDATION_CHECK_DIR}/toy.allow_same.gff3.diff"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        : > "{output}"
        tmpdir="$(mktemp -d)"
        trap 'rm -rf "$tmpdir"' EXIT

        # AA: EcoRI/EcoRI adjacency is suppressed by default when the second enzyme is absent.
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes EcoRI,NcoI \
          -min 1 -max 1000 \
          -gff "$tmpdir/aa.default.gff3" \
          -threads 1
        cat > "$tmpdir/expected.default.gff3" <<'EOF'
##gff-version 3
EOF
        diff -u "$tmpdir/expected.default.gff3" "$tmpdir/aa.default.gff3" >> "{output}"

        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes EcoRI,NcoI \
          -allow-same \
          -min 1 -max 1000 \
          -gff "$tmpdir/aa.allow_same.gff3" \
          -threads 1
        cat > "$tmpdir/expected.aa.gff3" <<'EOF'
##gff-version 3
chr1	radigest	fragment	6	16	.	+	.	ID=chr1_1;Length=11
EOF
        diff -u "$tmpdir/expected.aa.gff3" "$tmpdir/aa.allow_same.gff3" >> "{output}"

        # BB: same-type fragments are symmetric if the absent enzyme is listed first.
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes NcoI,EcoRI \
          -allow-same \
          -min 1 -max 1000 \
          -gff "$tmpdir/bb.allow_same.gff3" \
          -threads 1
        diff -u "$tmpdir/expected.aa.gff3" "$tmpdir/bb.allow_same.gff3" >> "{output}"
        """


rule validation_repro_digest:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_REPRO_FASTA
    output:
        json=f"{VALIDATION_DIR}/repro.t{{thread}}.json",
        gff=f"{VALIDATION_DIR}/repro.t{{thread}}.gff3",
        tsv=f"{VALIDATION_DIR}/repro.t{{thread}}.tsv"
    params:
        enzymes=VALIDATION_REPRO_ENZYMES,
        hard_min=VALIDATION_HARD_MIN,
        hard_max=VALIDATION_HARD_MAX,
        score_min=VALIDATION_SCORE_MIN,
        score_max=VALIDATION_SCORE_MAX
    threads: lambda wildcards: int(wildcards.thread)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.json})"
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes "{params.enzymes}" \
          -min {params.hard_min} \
          -max {params.hard_max} \
          -score-min {params.score_min} \
          -score-max {params.score_max} \
          -size-model hard \
          -json "{output.json}" \
          -gff "{output.gff}" \
          -fragments-tsv "{output.tsv}" \
          -threads {wildcards.thread}
        """


rule validation_repro_gff3_tsv_byte_identity:
    input:
        gff1=f"{VALIDATION_DIR}/repro.t1.gff3",
        gff2=f"{VALIDATION_DIR}/repro.t2.gff3",
        gff8=f"{VALIDATION_DIR}/repro.t8.gff3",
        tsv1=f"{VALIDATION_DIR}/repro.t1.tsv",
        tsv2=f"{VALIDATION_DIR}/repro.t2.tsv",
        tsv8=f"{VALIDATION_DIR}/repro.t8.tsv"
    output:
        f"{VALIDATION_CHECK_DIR}/repro.gff3_tsv.byte_identity.diff"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        : > "{output}"
        for pair in \
          "{input.gff1} {input.gff2}" \
          "{input.gff1} {input.gff8}" \
          "{input.tsv1} {input.tsv2}" \
          "{input.tsv1} {input.tsv8}"; do
          set -- $pair
          if ! cmp -s "$1" "$2"; then
            echo "byte mismatch: $1 vs $2" >> "{output}"
            sha256sum "$1" "$2" >> "{output}"
            exit 1
          fi
        done
        """


rule validation_repro_json_normalized_identity:
    input:
        json1=f"{VALIDATION_DIR}/repro.t1.json",
        json2=f"{VALIDATION_DIR}/repro.t2.json",
        json8=f"{VALIDATION_DIR}/repro.t8.json"
    output:
        f"{VALIDATION_CHECK_DIR}/repro.json.normalized_identity.diff"
    run:
        import difflib
        import json
        from pathlib import Path

        output_path = Path(str(output[0]))
        output_path.parent.mkdir(parents=True, exist_ok=True)

        def normalized(path):
            with open(path, "r", encoding="utf-8") as handle:
                doc = json.load(handle)
            doc.pop("command", None)
            doc.pop("outputs", None)
            for key in ("gff", "fragments_tsv", "fragments_fasta"):
                doc.pop(key, None)
            params = doc.get("parameters")
            if isinstance(params, dict):
                params.pop("threads", None)
            return json.dumps(doc, sort_keys=True, indent=2).splitlines(keepends=True)

        baseline = normalized(str(input.json1))
        failures = []
        for candidate in (str(input.json2), str(input.json8)):
            candidate_norm = normalized(candidate)
            diff = list(difflib.unified_diff(
                baseline,
                candidate_norm,
                fromfile=str(input.json1) + " (normalized)",
                tofile=candidate + " (normalized)",
            ))
            failures.extend(diff)
        output_path.write_text("".join(failures), encoding="utf-8")
        if failures:
            raise ValueError("normalized JSON summaries differ across thread counts")


rule validation_repro_json_same_thread_rerun:
    input:
        radigest=RADIGEST,
        fasta=VALIDATION_REPRO_FASTA
    output:
        f"{VALIDATION_CHECK_DIR}/repro.json.same_thread_rerun.diff"
    params:
        enzymes=VALIDATION_REPRO_ENZYMES,
        hard_min=VALIDATION_HARD_MIN,
        hard_max=VALIDATION_HARD_MAX,
        score_min=VALIDATION_SCORE_MIN,
        score_max=VALIDATION_SCORE_MAX
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        : > "{output}"
        tmpdir="$(mktemp -d)"
        trap 'rm -rf "$tmpdir"' EXIT
        json_path="$tmpdir/repro.same_thread.json"
        first_path="$tmpdir/repro.same_thread.first.json"

        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes "{params.enzymes}" \
          -min {params.hard_min} \
          -max {params.hard_max} \
          -score-min {params.score_min} \
          -score-max {params.score_max} \
          -size-model hard \
          -json "$json_path" \
          -threads 1
        cp "$json_path" "$first_path"
        {input.radigest} \
          -fasta "{input.fasta}" \
          -enzymes "{params.enzymes}" \
          -min {params.hard_min} \
          -max {params.hard_max} \
          -score-min {params.score_min} \
          -score-max {params.score_max} \
          -size-model hard \
          -json "$json_path" \
          -threads 1
        if ! cmp -s "$first_path" "$json_path"; then
          echo "same-thread JSON rerun byte mismatch" >> "{output}"
          diff -u "$first_path" "$json_path" >> "{output}" || true
          exit 1
        fi
        """


rule validation_repro_sha256:
    input:
        json=REPRO_JSON,
        gff=REPRO_GFF3,
        tsv=REPRO_TSV
    output:
        f"{VALIDATION_DIR}/repro.sha256"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        {{
          printf '# JSON contains command/output/thread provenance, so cross-thread JSON hashes are expected to differ.\n'
          printf '# GFF3 and TSV hashes must match across thread counts; normalized JSON identity is checked separately.\n'
          sha256sum {input.json} {input.gff} {input.tsv}
        }} > "{output}"
        """


rule software_validation_table:
    input:
        script="workflow/scripts/write_software_validation_table.py",
        go_test_log=config["provenance"]["go_test_log"],
        toy_fasta=VALIDATION_TOY_FASTA,
        expected_toy_double_gff3=VALIDATION_EXPECTED_TOY_DOUBLE_GFF3,
        toy_single_json=TOY_SINGLE_JSON,
        toy_single_gff3=TOY_SINGLE_GFF3,
        toy_double_json=TOY_DOUBLE_JSON,
        toy_double_gff3=TOY_DOUBLE_GFF3,
        toy_double_tsv=TOY_DOUBLE_TSV,
        toy_double_fasta=TOY_DOUBLE_FASTA,
        toy_include_ends_gff3=TOY_INCLUDE_ENDS_GFF3,
        toy_single_diff=f"{VALIDATION_CHECK_DIR}/toy.single.gff3.diff",
        toy_double_diff=f"{VALIDATION_CHECK_DIR}/toy.double.gff3.diff",
        toy_include_ends_diff=f"{VALIDATION_CHECK_DIR}/toy.include_ends.gff3.diff",
        toy_allow_same_diff=f"{VALIDATION_CHECK_DIR}/toy.allow_same.gff3.diff",
        repro_json=REPRO_JSON,
        repro_gff=REPRO_GFF3,
        repro_tsv=REPRO_TSV,
        repro_gff_tsv_diff=f"{VALIDATION_CHECK_DIR}/repro.gff3_tsv.byte_identity.diff",
        repro_json_normalized_diff=f"{VALIDATION_CHECK_DIR}/repro.json.normalized_identity.diff",
        repro_json_same_thread_diff=f"{VALIDATION_CHECK_DIR}/repro.json.same_thread_rerun.diff",
        repro_sha256=f"{VALIDATION_DIR}/repro.sha256"
    output:
        SOFTWARE_VALIDATION_TABLE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        python3 "{input.script}" \
          --out "{output}" \
          --go-test-log "{input.go_test_log}" \
          --toy-fasta "{input.toy_fasta}" \
          --toy-single-json "{input.toy_single_json}" \
          --toy-single-gff3 "{input.toy_single_gff3}" \
          --toy-double-json "{input.toy_double_json}" \
          --toy-double-gff3 "{input.toy_double_gff3}" \
          --toy-double-tsv "{input.toy_double_tsv}" \
          --toy-double-fasta "{input.toy_double_fasta}" \
          --toy-include-ends-gff3 "{input.toy_include_ends_gff3}" \
          --toy-single-diff "{input.toy_single_diff}" \
          --toy-double-diff "{input.toy_double_diff}" \
          --toy-include-ends-diff "{input.toy_include_ends_diff}" \
          --toy-allow-same-diff "{input.toy_allow_same_diff}" \
          --repro-gff-tsv-diff "{input.repro_gff_tsv_diff}" \
          --repro-json-normalized-diff "{input.repro_json_normalized_diff}" \
          --repro-json-same-thread-diff "{input.repro_json_same_thread_diff}" \
          --repro-sha256 "{input.repro_sha256}"
        """
