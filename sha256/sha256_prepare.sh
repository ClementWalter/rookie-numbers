#!/usr/bin/env bash
# Prepare SHA256 STWO prover state for csp-benchmarks
set -eo pipefail

# Check required environment variables
if [[ -z "${INPUT_SIZE:-}" ]]; then
    echo "Error: INPUT_SIZE environment variable is required" >&2
    echo "Usage: INPUT_SIZE=128 STATE_JSON=/tmp/state.json $0" >&2
    exit 1
fi

if [[ -z "${STATE_JSON:-}" ]]; then
    echo "Error: STATE_JSON environment variable is required" >&2
    echo "Usage: INPUT_SIZE=128 STATE_JSON=/tmp/state.json $0" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

# Build if not present
if [[ ! -f "$BINARY" ]]; then
    echo "Building sha256_prover..."
    cargo build --release -p sha256 --bin sha256_prover --manifest-path "$SCRIPT_DIR/../Cargo.toml"
fi

exec "$BINARY" prepare --input-size "$INPUT_SIZE" --state-json "$STATE_JSON"
