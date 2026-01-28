#![allow(non_camel_case_types)]
#![feature(
    portable_simd,
    array_chunks,
    iter_array_chunks,
    macro_metavar_expr_concat
)]
pub mod components;
pub mod macros;
pub mod partitions;
pub mod preprocessed;
pub mod relations;
pub mod sha256;

// Re-export stwo types for external consumers
#[cfg(feature = "smalloc")]
use smalloc::Smalloc;
pub use stwo::{
    core::{pcs::PcsConfig, proof::StarkProof, vcs::blake2_merkle::Blake2sMerkleHasher},
    prover::backend::Column,
};
#[cfg(feature = "smalloc")]
#[global_allocator]
static GLOBAL: Smalloc = Smalloc::new();

#[cfg(feature = "smalloc")]
#[ctor::ctor]
unsafe fn init_smalloc() {
    GLOBAL.init();
}

#[cfg(feature = "peak-alloc")]
use peak_alloc::PeakAlloc;
#[cfg(feature = "peak-alloc")]
#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

#[cfg(feature = "jemalloc")]
use tikv_jemallocator::Jemalloc;
#[cfg(feature = "jemalloc")]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use num_traits::Zero;
use stwo::{
    core::{
        channel::Blake2sChannel,
        fields::qm31::SecureField,
        pcs::CommitmentSchemeVerifier,
        poly::circle::CanonicCoset,
        vcs::blake2_merkle::Blake2sMerkleChannel,
        verifier::{verify, VerificationError},
    },
    prover::{backend::simd::SimdBackend, poly::circle::PolyOps, prove, CommitmentSchemeProver},
};
use stwo_constraint_framework::TraceLocationAllocator;
use tracing::{debug, info, span, Level};

use crate::{
    components::{gen_interaction_trace, gen_trace},
    preprocessed::PreProcessedTrace,
    relations::Relations,
};

/// Maximum log_size when `dynamic-preprocessed-shape`.
pub const MAX_PREPROCESSED_LOG_SIZE: u32 = 21;

/// Returns the effective log_size for preprocessed column chunking.
///
/// When `dynamic-preprocessed-shape` feature is enabled, returns the actual log_size,
/// which optimizes for large proofs by matching preprocessing to the exact proof size.
///
/// When disabled (default), returns max(log_size, MAX_PREPROCESSED_LOG_SIZE) so small proofs
/// use a fixed preprocessed size (one chunk, unsuffixed IDs) while larger proofs scale.
#[cfg(feature = "dynamic-preprocessed-shape")]
pub fn preprocessed_log_size(log_size: u32) -> u32 {
    log_size
}

#[cfg(not(feature = "dynamic-preprocessed-shape"))]
pub fn preprocessed_log_size(log_size: u32) -> u32 {
    MAX_PREPROCESSED_LOG_SIZE.max(log_size)
}

/// Log size for twiddles when dynamic-preprocessed-shape is disabled.
/// This matches the original hardcoded value before the dynamic reshape feature.
#[cfg(not(feature = "dynamic-preprocessed-shape"))]
const TWIDDLES_LOG_SIZE: u32 = 21;

pub fn prove_sha256(
    log_size: u32,
    config: PcsConfig,
) -> (StarkProof<Blake2sMerkleHasher>, components::ClaimedSum) {
    // Precompute twiddles.
    let span = span!(Level::INFO, "Precompute twiddles").entered();
    #[cfg(feature = "dynamic-preprocessed-shape")]
    let twiddles_log_size = preprocessed_log_size(log_size).max(log_size);
    #[cfg(not(feature = "dynamic-preprocessed-shape"))]
    let twiddles_log_size = TWIDDLES_LOG_SIZE.max(log_size);
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(twiddles_log_size + config.fri_config.log_blowup_factor + 2)
            .circle_domain()
            .half_coset,
    );
    span.exit();

    // Setup protocol.
    let channel = &mut Blake2sChannel::default();
    let mut commitment_scheme =
        CommitmentSchemeProver::<_, Blake2sMerkleChannel>::new(config, &twiddles);

    // Preprocessed trace.
    let span = span!(Level::INFO, "Constant").entered();
    let span_1 = span!(Level::INFO, "Simd generation").entered();
    let preprocessed_trace = PreProcessedTrace::new(log_size);
    span_1.exit();
    let span_2 = span!(Level::INFO, "Extend evals").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.trace);
    tree_builder.commit(channel);
    span_2.exit();
    span.exit();

    // Trace.
    let span = span!(Level::INFO, "Trace").entered();
    let (trace, lookup_data) = gen_trace(log_size);
    let span_1 = span!(Level::INFO, "Extend evals").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);
    span_1.exit();
    span.exit();

    // Draw lookup elements.
    let relations = Relations::draw(channel);

    // Interaction trace.
    let span = span!(Level::INFO, "Interaction").entered();
    let (trace, claimed_sum) = gen_interaction_trace(lookup_data, &relations);
    let span_1 = span!(Level::INFO, "Extend evals").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);
    span_1.exit();
    span.exit();

    debug!(
        "Columns count: {:?}",
        commitment_scheme
            .trees
            .as_ref()
            .map(|tree| tree.polynomials.len())
    );
    debug!(
        "Columns length: {:?}",
        commitment_scheme.trees.as_ref().map(|tree| {
            let max_len = tree
                .polynomials
                .iter()
                .map(|poly| poly.evals.values.len().ilog2())
                .collect::<Vec<_>>()
                .iter()
                .copied()
                .max()
                .unwrap();
            assert!(max_len <= log_size + 1);
            max_len
        })
    );

    // Prove constraints.
    let span = span!(Level::INFO, "Prove").entered();
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preprocessed_columns(&preprocessed_trace.ids);
    let components =
        components::Components::new(log_size, trace_allocator, &relations, &claimed_sum);

    #[cfg(feature = "track-relations")]
    println!(
        "Trace log degree bounds: {:?}",
        components.trace_log_degree_bounds()
    );

    if claimed_sum.scheduling + claimed_sum.compression + claimed_sum.preprocessed.sum()
        != SecureField::zero()
    {
        #[cfg(feature = "track-relations")]
        println!(
            "Relation summary: {:?}",
            components.track_relations(&commitment_scheme)
        );
        panic!(
            "Relation summary is not zero: {}",
            claimed_sum.scheduling + claimed_sum.compression + claimed_sum.preprocessed.sum()
        );
    }

    let proof = prove(&components.provers(), channel, commitment_scheme);
    if let Err(e) = proof {
        panic!("Proof error: {e:?}");
    }
    span.exit();

    (proof.unwrap(), claimed_sum)
}

/// Verify a SHA256 proof.
///
/// # Arguments
/// * `proof` - The STARK proof to verify
/// * `log_size` - The log2 of the number of SHA256 instances
/// * `claimed_sum` - The claimed sums from the prover
///
/// # Returns
/// * `Ok(())` if the proof is valid
/// * `Err(VerificationError)` if the proof is invalid
pub fn verify_sha256(
    proof: StarkProof<Blake2sMerkleHasher>,
    log_size: u32,
    claimed_sum: &components::ClaimedSum,
) -> Result<(), VerificationError> {
    // Verify that claimed sums sum to zero (logup protocol constraint)
    if claimed_sum.scheduling + claimed_sum.compression + claimed_sum.preprocessed.sum()
        != SecureField::zero()
    {
        return Err(VerificationError::OodsNotMatching);
    }

    // Setup verification protocol (mirror of prover setup)
    let channel = &mut Blake2sChannel::default();
    let commitment_scheme =
        &mut CommitmentSchemeVerifier::<Blake2sMerkleChannel>::new(proof.config);

    // Get preprocessed column IDs for trace allocation
    let preprocessed_trace = PreProcessedTrace::new(log_size);

    // Commit preprocessed trace (tree 0)
    // Get sizes from the preprocessed trace
    let preprocessed_sizes: Vec<u32> = preprocessed_trace
        .trace
        .iter()
        .map(|eval| eval.domain.log_size())
        .collect();
    commitment_scheme.commit(proof.commitments[0], &preprocessed_sizes, channel);

    // Create components to get trace sizes (need relations for component creation)
    // Use dummy values for claimed_sum fields since we only need sizes
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preprocessed_columns(&preprocessed_trace.ids);
    let dummy_relations = Relations::draw(&mut Blake2sChannel::default());
    let temp_components =
        components::Components::new(log_size, trace_allocator, &dummy_relations, claimed_sum);

    // Get trace log degree bounds from all components
    let all_bounds = temp_components.trace_log_degree_bounds();

    // Collect sizes for trace tree (tree 1) - all columns from all components
    let mut trace_sizes: Vec<u32> = Vec::new();
    for component_bounds in &all_bounds {
        if component_bounds.len() > 1 {
            trace_sizes.extend(&component_bounds[1]);
        }
    }

    // Commit main trace (tree 1)
    commitment_scheme.commit(proof.commitments[1], &trace_sizes, channel);

    // Draw relations AFTER trace commit (must match prover order)
    let relations = Relations::draw(channel);

    // Recreate components with correct relations
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preprocessed_columns(&preprocessed_trace.ids);
    let components =
        components::Components::new(log_size, trace_allocator, &relations, claimed_sum);
    let all_bounds = components.trace_log_degree_bounds();

    // Collect sizes for interaction tree (tree 2)
    let mut interaction_sizes: Vec<u32> = Vec::new();
    for component_bounds in &all_bounds {
        if component_bounds.len() > 2 {
            interaction_sizes.extend(&component_bounds[2]);
        }
    }

    // Commit interaction trace (tree 2)
    commitment_scheme.commit(proof.commitments[2], &interaction_sizes, channel);

    // Get all components for verification
    let component_refs: Vec<&dyn stwo::core::air::Component> = components.components();

    // Verify the proof
    verify(&component_refs, channel, commitment_scheme, proof)
}

pub fn print_enabled_features() {
    let features: Vec<&str> = vec![
        #[cfg(feature = "parallel")]
        "Stwo parallel",
        #[cfg(not(feature = "parallel"))]
        "Stwo non-parallel",
        #[cfg(feature = "peak-alloc")]
        "peak-alloc",
        #[cfg(feature = "jemalloc")]
        "jemalloc",
        #[cfg(feature = "dynamic-preprocessed-shape")]
        "dynamic-preprocessed-shape",
    ];

    if features.is_empty() {
        info!("Features: (none)");
    } else {
        info!("Features: {}", features.join(", "));
    }
}
#[cfg(test)]
mod tests {
    use std::{env, time::Instant};

    use rayon::iter::{IntoParallelIterator, ParallelIterator};

    use super::*;

    /// Print all enabled features

    #[test_log::test]
    fn test_prove_sha256() {
        print_enabled_features();

        // Get from environment variable:
        let log_n_instances = env::var("LOG_N_INSTANCES")
            .unwrap_or_else(|_| "13".to_string())
            .parse::<u32>()
            .unwrap();
        let n_iter = env::var("N_ITER")
            .unwrap_or_else(|_| "1".to_string())
            .parse::<u32>()
            .unwrap();
        let log_size = log_n_instances;

        info!("Log size: {}", log_size);
        info!("Number of iterations: {}", n_iter);

        #[cfg(feature = "peak-alloc")]
        PEAK_ALLOC.reset_peak_usage();
        let span = span!(Level::INFO, "Prove").entered();

        let start = Instant::now();
        (0..n_iter)
            .into_par_iter()
            .map(|_| prove_sha256(log_size, PcsConfig::default()))
            .collect::<Vec<_>>();
        span.exit();
        info!(
            "Throughput {:?}",
            (1 << log_n_instances) as f32 * n_iter as f32 / start.elapsed().as_secs() as f32
        );

        #[cfg(feature = "peak-alloc")]
        {
            let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
            info!("Peak memory: {} MB", peak_bytes);
        }
    }

    #[test_log::test]
    fn test_prove_verify_roundtrip() {
        print_enabled_features();

        let log_size = 8; // Small size for fast testing
        info!("Testing prove/verify roundtrip with log_size={}", log_size);

        // Prove
        let (proof, claimed_sum) = prove_sha256(log_size, PcsConfig::default());
        info!("Proof generated successfully");

        // Verify
        let result = verify_sha256(proof, log_size, &claimed_sum);
        assert!(result.is_ok(), "Verification failed: {:?}", result.err());
        info!("Proof verified successfully");
    }
}
