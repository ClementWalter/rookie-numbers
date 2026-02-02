use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sha256::{preprocess_sha256, print_enabled_features, prove_sha256, PreprocessedData};
use stwo::core::pcs::PcsConfig;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

fn main() {
    tracing_subscriber::registry()
        .with(EnvFilter::from_default_env())
        .with(tracing_subscriber::fmt::layer())
        .init();
    divan::main();
}

const N_ITER: &[usize] = &[6, 7, 8];

#[divan::bench(
    consts = N_ITER,
    args = [13, 14],
    sample_count = 1
)]
fn bench_sha256<const N_ITER: usize>(bencher: divan::Bencher, log_size: u32) {
    print_enabled_features();

    let config = PcsConfig::default();
    let input = vec![0xab_u8; (1 << log_size) * 512];

    // Preprocess once before the benchmark (not timed)
    tracing::info!("Preprocessing for log_size={}...", log_size);
    let preprocessed: PreprocessedData = preprocess_sha256(log_size, config);
    tracing::info!("Preprocessing complete");

    bencher.bench(|| {
        #[cfg(feature = "peak-alloc")]
        PEAK_ALLOC.reset_peak_usage();
        (0..N_ITER)
            .into_par_iter()
            .map(|_| prove_sha256(&input, config, &preprocessed))
            .collect::<Vec<_>>();
        #[cfg(feature = "peak-alloc")]
        {
            let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
            tracing::info!("Peak memory: {peak_bytes} MB");
            divan::black_box(peak_bytes);
        }
    });
}
