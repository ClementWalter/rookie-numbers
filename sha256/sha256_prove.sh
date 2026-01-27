#!/usr/bin/env bash
# Run SHA256 STWO prover for csp-benchmarks
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

if [[ ! -f "$BINARY" ]]; then
    echo "Error: Binary not found. Run sha256_prepare.sh first." >&2
    exit 1
fi

exec "$BINARY" prove --state-json "$STATE_JSON"
