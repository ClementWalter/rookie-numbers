# Wat

Last year, [Shahar Papini tweet](https://x.com/PapiniShahar/status/1831402791400812624) suggested that it's possible to reach 10Mhz using Stwo. This repo is for benchmarking different configuration of constraints and interactions traces to see what can actually be achieved.

## How to run

```bash
RUSTFLAGS="-C target-cpu=native" cargo test -r -- test_prove_rookie
```
