#!/bin/bash
# Unit tests for workflow/make_dag.sh
# Tests argument validation, DAG file structure, and variable substitution.

SCRIPT="$(cd "$(dirname "$0")/.." && pwd)/workflow/make_dag.sh"
PASS=0
FAIL=0

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

assert_eq() {
    local desc="$1" expected="$2" actual="$3"
    if [[ "$actual" == "$expected" ]]; then
        echo "  PASS: $desc"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $desc"
        echo "        Expected: '$expected'"
        echo "        Actual:   '$actual'"
        FAIL=$((FAIL + 1))
    fi
}

assert_contains() {
    local desc="$1" pattern="$2" file="$3"
    if grep -qF "$pattern" "$file" 2>/dev/null; then
        echo "  PASS: $desc"
        PASS=$((PASS + 1))
    else
        echo "  FAIL: $desc (pattern '$pattern' not found)"
        FAIL=$((FAIL + 1))
    fi
}

assert_exit_nonzero() {
    local desc="$1"
    shift
    if "$@" >/dev/null 2>&1; then
        echo "  FAIL: $desc (expected non-zero exit, got 0)"
        FAIL=$((FAIL + 1))
    else
        echo "  PASS: $desc"
        PASS=$((PASS + 1))
    fi
}

# ---------------------------------------------------------------------------
# Setup: temporary working directory mimicking expected layout
# ---------------------------------------------------------------------------

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$TMPDIR/input_files" "$TMPDIR/workflow"

# Dummy input files
touch "$TMPDIR/input_files/sample.ms2"
touch "$TMPDIR/input_files/sample.mgf"
touch "$TMPDIR/input_files/sample.mzML"

# Symlink the workflow scripts so make_dag.sh can reference its .sub files
cp "$SCRIPT" "$TMPDIR/workflow/make_dag.sh"
for sub in endogenius_round1.sub endogenius_round2.sub casanovo_target.sub \
           casanovo_decoy.sub generate_decoy.sub filter_casanovo.sub; do
    touch "$TMPDIR/workflow/$sub"
done

RUN_SCRIPT="bash $TMPDIR/workflow/make_dag.sh"

# ---------------------------------------------------------------------------
# Argument validation tests
# ---------------------------------------------------------------------------

echo "=== Argument validation ==="

assert_exit_nonzero "exits non-zero with no arguments" \
    bash "$TMPDIR/workflow/make_dag.sh"

assert_exit_nonzero "exits non-zero when -o is missing" \
    bash "$TMPDIR/workflow/make_dag.sh" \
    -s sample -j sample_round1.json -k sample_round2.json \
    -a EG_database.fasta -m motif.csv -f 0.05 -d filter.fasta

assert_exit_nonzero "exits non-zero when -s is missing" \
    bash "$TMPDIR/workflow/make_dag.sh" \
    -j sample_round1.json -k sample_round2.json \
    -a EG_database.fasta -m motif.csv -f 0.05 -d filter.fasta -o out

assert_exit_nonzero "exits non-zero when -f (FDR) is missing" \
    bash "$TMPDIR/workflow/make_dag.sh" \
    -s sample -j sample_round1.json -k sample_round2.json \
    -a EG_database.fasta -m motif.csv -d filter.fasta -o out

# ---------------------------------------------------------------------------
# DAG generation tests
# ---------------------------------------------------------------------------

echo "=== DAG generation ==="

pushd "$TMPDIR" > /dev/null
bash "$TMPDIR/workflow/make_dag.sh" \
    -s sample \
    -j sample_round1.json \
    -k sample_round2.json \
    -a EG_database.fasta \
    -m motif.csv \
    -f 0.05 \
    -d filter.fasta \
    -o test_pipeline \
    > /dev/null 2>&1
popd > /dev/null

DAG="$TMPDIR/test_pipeline.dag"

assert_eq "DAG file is created" "0" "$([ -f "$DAG" ] && echo 0 || echo 1)"

# ---------------------------------------------------------------------------
# Job presence tests
# ---------------------------------------------------------------------------

echo "=== Job presence ==="

assert_contains "contains ENDOGENIUS_R1 job"   "JOB ENDOGENIUS_R1"   "$DAG"
assert_contains "contains GENERATE_DECOY job"  "JOB GENERATE_DECOY"  "$DAG"
assert_contains "contains CASANOVO_TARGET job" "JOB CASANOVO_TARGET" "$DAG"
assert_contains "contains CASANOVO_DECOY job"  "JOB CASANOVO_DECOY"  "$DAG"
assert_contains "contains ANALYZE_FDR job"     "JOB ANALYZE_FDR"     "$DAG"
assert_contains "contains ENDOGENIUS_R2 job"   "JOB ENDOGENIUS_R2"   "$DAG"

JOB_COUNT=$(grep -c "^JOB " "$DAG")
assert_eq "exactly 6 jobs in DAG" "6" "$JOB_COUNT"

# ---------------------------------------------------------------------------
# Dependency tests
# ---------------------------------------------------------------------------

echo "=== Dependencies ==="

assert_contains "EG_R1 -> CASANOVO_TARGET"              "PARENT ENDOGENIUS_R1 CHILD CASANOVO_TARGET"           "$DAG"
assert_contains "GENERATE_DECOY -> CASANOVO_DECOY"      "PARENT GENERATE_DECOY CHILD CASANOVO_DECOY"           "$DAG"
assert_contains "TARGET+DECOY -> ANALYZE_FDR"           "PARENT CASANOVO_TARGET CASANOVO_DECOY CHILD ANALYZE_FDR" "$DAG"
assert_contains "ANALYZE_FDR -> ENDOGENIUS_R2"          "PARENT ANALYZE_FDR CHILD ENDOGENIUS_R2"               "$DAG"

# ---------------------------------------------------------------------------
# Variable substitution tests
# ---------------------------------------------------------------------------

echo "=== Variable substitution ==="

assert_contains "FDR value written to DAG"              'fdr="0.05"'            "$DAG"
assert_contains "sample name in MS2 path"               'input_files/sample.ms2' "$DAG"
assert_contains "sample name in MGF path"               'input_files/sample.mgf' "$DAG"
assert_contains "EG FASTA in round1 VARS"               'EG_database.fasta'     "$DAG"
assert_contains "motif CSV in VARS"                     'motif.csv'             "$DAG"
assert_contains "filter FASTA in ANALYZE_FDR VARS"      'filter.fasta'          "$DAG"
assert_contains "output rooted under all_results"       'all_results/sample'    "$DAG"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]]
