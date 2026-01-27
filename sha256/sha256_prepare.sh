#!/usr/bin/env bash
# Prepare SHA256 STWO prover state for csp-benchmarks
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

# Build if not present
if [[ ! -f "$BINARY" ]]; then
    echo "Building sha256_prover..."
    cargo build --release -p sha256 --bin sha256_prover --manifest-path "$SCRIPT_DIR/../Cargo.toml"
fi

exec "$BINARY" prepare --input-size "$INPUT_SIZE" --state-json "$STATE_JSON"
