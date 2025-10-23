# Wat

Last year,
[Shahar Papini tweet](https://x.com/PapiniShahar/status/1831402791400812624)
suggested that it's possible to reach 10Mhz using Stwo. This repo is for
benchmarking different configuration of stwo, as well as some utils and actual
meaningful AIR implementations.

## How to run

### Theoretical maximum frequency benchmarks

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench frequency
```

Some results can also be found in
[this Google Sheet](https://docs.google.com/spreadsheets/d/1MEiPB4X7zjQgXYMV5Uk0t0JzbnBf024zYWQTREIyj8Q/edit?usp=sharing).

### Actual AIR implementations benchmarks

#### Sha256

To bench several configurations:

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench sha256
```

To run a single test:

```bash
LOG_N_INSTANCES=17 N_ITER=3 RUSTFLAGS="-C target-cpu=native" cargo t -r test_prove_sha256
```
