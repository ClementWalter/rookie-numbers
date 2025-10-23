use peak_alloc::PeakAlloc;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sha256::prove_sha256;
use stwo::core::pcs::PcsConfig;
use tracing::info;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() {
    divan::main();
}

const N_ITER: &[usize] = &[6, 7, 8];

#[divan::bench(
    consts = N_ITER,
    args = [13, 14],
    sample_count = 1
)]
fn bench_sha256<const N_ITER: usize>(bencher: divan::Bencher, log_size: u32) {
    #[cfg(feature = "parallel")]
    info!("Stwo Parallel");
    #[cfg(not(feature = "parallel"))]
    info!("Stwo Non-parallel");

    bencher.bench(|| {
        PEAK_ALLOC.reset_peak_usage();
        (0..N_ITER)
            .into_par_iter()
            .map(|_| prove_sha256(log_size, PcsConfig::default()))
            .collect::<Vec<_>>();
        let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
        println!("Peak memory: {peak_bytes} MB");
        divan::black_box(peak_bytes);
    });
}
