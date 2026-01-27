#!/usr/bin/env bash
# Run SHA256 STWO prover for csp-benchmarks
set -eo pipefail

# Check required environment variables
if [[ -z "${STATE_JSON:-}" ]]; then
    echo "Error: STATE_JSON environment variable is required" >&2
    echo "Usage: STATE_JSON=/tmp/state.json $0" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

if [[ ! -f "$BINARY" ]]; then
    echo "Error: Binary not found. Run sha256_prepare.sh first." >&2
    exit 1
fi

exec "$BINARY" prove --state-json "$STATE_JSON"
