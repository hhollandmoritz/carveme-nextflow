#!/usr/bin/env bash

# runs a dry run of several pipeline calls; also known as "stub tests"

set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"

test_root="test-results/automated"

fail() {
    printf 'TEST FAILED: %s\n' "$1" >&2
    exit 1
}

check_command() {
    command -v "$1" >/dev/null 2>&1 ||
        fail "Required command not found: $1"
}

count_outputs() {
    local directory="$1"
    local pattern="$2"

    if [[ ! -d "$directory" ]]; then
        printf 'Expected output directory does not exist: %s\n' \
            "$directory" >&2
        printf '0\n'
        return
    fi

    find -L "$directory" \
        -maxdepth 1 \
        -type f \
        -name "$pattern" \
        -print |
        wc -l
}

run_stub_test() {
    local test_name="$1"
    local params_file="$2"
    local expected_models="$3"
    local expected_logs="$4"
    local output_dir="${test_root}/${test_name}"

    printf '\nRunning stub test: %s\n' "$test_name"

    nextflow run main.nf \
        -profile apptainer \
        -params-file "$params_file" \
        --outdir "$output_dir" \
        -stub-run

    model_count="$(count_outputs "$output_dir/models" '*.xml')"
    log_count="$(count_outputs "$output_dir/logs" '*.log')"

    if [[ "$model_count" -ne "$expected_models" ]]; then
        printf 'TEST FAILED: %s: expected %s model(s), found %s\n' \
            "$test_name" \
            "$expected_models" \
            "$model_count" >&2

        find "$output_dir/models" \
            -maxdepth 1 \
            -printf '%y %p -> %l\n' >&2

        return 1
    fi

    if [[ "$log_count" -ne "$expected_logs" ]]; then
        printf 'TEST FAILED: %s: expected %s log(s), found %s\n' \
            "$test_name" \
            "$expected_logs" \
            "$log_count" >&2

        find "$output_dir/logs" \
            -maxdepth 1 \
            -printf '%y %p -> %l\n' >&2

        return 1
    fi

    printf 'PASSED: %s\n' "$test_name"
    return 0

}

check_command nextflow
check_command apptainer

test -f containers/carveme-1.6.6.sif ||
    fail "Container not found: containers/carveme-1.6.6.sif"

test -f test-data/*.faa ||
    fail "Test genome not found"

test -f assets/SDMM_media.tsv ||
    fail "SDMM media database not found"

printf 'Checking container...\n'

apptainer exec containers/carveme-1.6.6.sif \
    python -c \
    "import carveme, reframed, pyscipopt, libsbml"

printf 'Checking Nextflow configuration...\n'

nextflow inspect -profile apptainer . >/dev/null

failures=0

run_stub_test \
    "no-gapfill" \
    "test-params/no-gapfill.yml" \
    1 \
    1 ||
    failures=$((failures + 1))

run_stub_test \
    "m9" \
    "test-params/M9.yml" \
    1 \
    1 ||
    failures=$((failures + 1))

run_stub_test \
    "sdmm" \
    "test-params/sdmm.yml" \
    1 \
    1 ||
    failures=$((failures + 1))

run_stub_test \
    "combined" \
    "test-params/combined_media.yml" \
    3 \
    3 ||
    failures=$((failures + 1))

if [[ "$failures" -gt 0 ]]; then
    printf '\n%s stub test(s) failed.\n' "$failures" >&2
    exit 1
fi

printf '\nAll stub tests passed.\n'
