use peak_alloc::PeakAlloc;
use prover::prover::prove_rookie;
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    divan::main();
}

const N: &[usize] = &[3, 3 * 10, 3 * 50, 3 * 150, 3 * 500, 3 * 1000, 3 * 2000];

#[divan::bench(
    consts = N,
    args = [15, 16, 17, 18, 19, 20],
    sample_count = 5
)]
fn bench_rookie<const N: usize>(bencher: divan::Bencher, log_size: u32) {
    bencher.bench(|| {
        PEAK_ALLOC.reset_peak_usage();
        let result = prove_rookie::<Blake2sMerkleChannel, N>(log_size);
        assert!(result.is_ok());
        let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
        println!("Peak memory: {} MB", peak_bytes);
        divan::black_box(peak_bytes);
    });
}
