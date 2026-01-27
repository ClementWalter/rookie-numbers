#!/usr/bin/env bash
# Report SHA256 STWO proof sizes for csp-benchmarks
set -eo pipefail

# Check required environment variables
if [[ -z "${STATE_JSON:-}" ]]; then
    echo "Error: STATE_JSON environment variable is required" >&2
    echo "Usage: STATE_JSON=/tmp/state.json SIZES_JSON=/tmp/sizes.json $0" >&2
    exit 1
fi

if [[ -z "${SIZES_JSON:-}" ]]; then
    echo "Error: SIZES_JSON environment variable is required" >&2
    echo "Usage: STATE_JSON=/tmp/state.json SIZES_JSON=/tmp/sizes.json $0" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARY="$SCRIPT_DIR/../target/release/sha256_prover"

if [[ ! -f "$BINARY" ]]; then
    echo "Error: Binary not found. Run sha256_prepare.sh first." >&2
    exit 1
fi

exec "$BINARY" measure --state-json "$STATE_JSON" --sizes-json "$SIZES_JSON"
