#!/usr/bin/env bash

# Runs a full CarveMe test using the test data.

set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"

output_dir="test-results/integration"
models_dir="${output_dir}/models"

printf 'Running real CarveMe integration test...\n'

nextflow run main.nf \
    -profile apptainer \
    -params-file test-params/no-gapfill.yml \
    --outdir "$output_dir"

if [[ ! -d "$models_dir" ]]; then
    printf 'TEST FAILED: model directory does not exist: %s\n' \
        "$models_dir" >&2
    exit 1
fi

model_count="$(
    find -L "$models_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.sbml' \
        -print |
        wc -l
)"

if [[ "$model_count" -ne 1 ]]; then
    printf 'TEST FAILED: expected one model, found %s\n' \
        "$model_count" >&2

    find "$models_dir" \
        -maxdepth 1 \
        -printf '%y %p -> %l\n' >&2

    exit 1
fi

model="$(
    find -L "$models_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.sbml' \
        -print |
        head -n 1
)"

if [[ ! -s "$model" ]]; then
    printf 'TEST FAILED: model is empty: %s\n' "$model" >&2
    exit 1
fi

printf 'Integration test passed: %s\n' "$model"
