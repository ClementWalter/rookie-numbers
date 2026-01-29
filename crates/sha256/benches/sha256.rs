use rayon::iter::{IntoParallelIterator, ParallelIterator};
use sha256::{print_enabled_features, prove_sha256};
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

    let input = vec![0xab_u8; (1 << log_size) * 512];
    bencher.bench(|| {
        #[cfg(feature = "peak-alloc")]
        PEAK_ALLOC.reset_peak_usage();
        (0..N_ITER)
            .into_par_iter()
            .map(|_| prove_sha256(&input, PcsConfig::default()))
            .collect::<Vec<_>>();
        #[cfg(feature = "peak-alloc")]
        {
            let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
            tracing::info!("Peak memory: {peak_bytes} MB");
            divan::black_box(peak_bytes);
        }
    });
}
