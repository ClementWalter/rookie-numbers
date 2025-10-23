use peak_alloc::PeakAlloc;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rookie::prover::prove_dummy;
use stwo::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    divan::main();
}

const N: &[usize] = &[2_usize.pow(12), 2_usize.pow(13)];

#[divan::bench(
    consts = N,
    args = [13, 14],
    sample_count = 1
)]
fn bench_frequency<const N: usize>(bencher: divan::Bencher, log_size: u32) {
    bencher.bench(|| {
        PEAK_ALLOC.reset_peak_usage();
        (0..10)
            .into_par_iter()
            .map(|_| prove_dummy::<Blake2sMerkleChannel, N>(log_size))
            .collect::<Vec<_>>();
        let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
        println!("Peak memory: {peak_bytes} MB");
        divan::black_box(peak_bytes);
    });
}
