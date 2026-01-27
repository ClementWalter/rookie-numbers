#!/usr/bin/env bash
# Report SHA256 STWO proof sizes for csp-benchmarks
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

if [[ ! -f "$BINARY" ]]; then
    echo "Error: Binary not found. Run sha256_prepare.sh first." >&2
    exit 1
fi

exec "$BINARY" measure --state-json "$STATE_JSON" --sizes-json "$SIZES_JSON"
