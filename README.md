# Wat

Last year,
[Shahar Papini tweet](https://x.com/PapiniShahar/status/1831402791400812624)
suggested that it's possible to reach 10Mhz using Stwo. This repo is for
benchmarking different configuration of stwo, as well as some utils and actual
meaningful AIR implementations.

## How to run

### Setup

- Make sure you have the correct toolchain installed.

```bash
rustup toolchain install nightly-2025-07-14
```

- Linter and formatter [trunk](https://trunk.io/)

- Stwo submodule

```bash
git submodule update --init --recursive
```

- Check everything is working

```bash
cargo b -r
trunk check --all
```

Note: `trunk check --all` generates as of today few warnings due to `stwo`
crate:

```text
  ISSUES

Cargo.lock:352:0
 352:0  high  `derivative` is unmaintained; consider      osv-scanner/RUSTSEC-2024-0388
              using an alternative. Current version is
              vulnerable: 2.2.0.
 689:0  high  'paste' - no longer maintained. Current     osv-scanner/RUSTSEC-2024-0436
              version is vulnerable: 1.0.15.

external/stwo/crates/constraint-framework/src/lib.rs:304:32
 304:32  medium  you are deriving `PartialEq` and    clippy
                 can implement `Eq`                   /derive_partial_eq_without_eq

external/stwo/crates/stwo/src/prover/backend/simd/grind.rs:137:17
 137:17  medium  unused variable: `prefixed_digest`  clippy/unused_variables

Checked 133 files
✖ 2 security issues
✖ 2 lint issues
```

### Theoretical maximum frequency benchmarks

This benchmarks lets you estimate the theoretical maximum frequency of the
prover on your machine based on the trace size.

It fills a random trace of the given size and enforces no constraints.

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench frequency
```

Some results can also be found in
[this Google Sheet](https://docs.google.com/spreadsheets/d/1MEiPB4X7zjQgXYMV5Uk0t0JzbnBf024zYWQTREIyj8Q/edit?usp=sharing).

### Actual AIR implementations benchmarks

#### Sha256

To run a single test:

```bash
RUST_LOG=sha256::tests=info,stwo=off RUST_LOG_SPAN_EVENTS=close LOG_N_INSTANCES=13 N_ITER=3 RUSTFLAGS="-C target-cpu=native" cargo t -r --package sha256 --lib -- tests::test_prove_verify_long_input --exact --nocapture
```

Depending on your machine (especially the number of cores and available memory),
you may find OOM errors. It is a good idea to start low instance count and one
single thread:

```bash
RUST_LOG=sha256::tests=info,stwo=off RUST_LOG_SPAN_EVENTS=close LOG_N_INSTANCES=13 N_ITER=1 RUSTFLAGS="-C target-cpu=native" cargo t -r --package sha256 --lib -- tests::test_prove_verify_long_input --exact --nocapture
```

to check memory usage and performance, and then increase the number of
iterations and log size. Best results should be obtained with `N_ITER` set to
the actual number of cores on your machine. For example, on macOs

```bash
RUST_LOG=sha256::tests=info,stwo=off RUST_LOG_SPAN_EVENTS=close LOG_N_INSTANCES=13 N_ITER=$(sysctl -n hw.logicalcpu) RUSTFLAGS="-C target-cpu=native" cargo t -r --package sha256 --lib -- tests::test_prove_verify_long_input --exact --nocapture
```

The global allocator can be changed by using the `peak-alloc` or `jemalloc`
features:

- `peak-alloc`: uses the [peak-alloc](https://crates.io/crates/peak-alloc) crate
  to track the peak memory usage.
- `jemalloc`: uses the [jemalloc](https://crates.io/crates/tikv-jemallocator)
  allocator that seems to provide better performance.

To bench several configurations:

```bash
RUSTFLAGS="-C target-cpu=native" cargo bench --bench sha256
```
