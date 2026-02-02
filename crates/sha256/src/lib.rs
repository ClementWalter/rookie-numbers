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
    core::{
        fri::FriConfig, pcs::PcsConfig, proof::StarkProof, vcs::blake2_merkle::Blake2sMerkleHasher,
    },
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
use serde::{Deserialize, Serialize};
use stwo::{
    core::{
        channel::{Blake2sChannel, MerkleChannel},
        fields::qm31::SecureField,
        pcs::CommitmentSchemeVerifier,
        poly::circle::CanonicCoset,
        vcs::{blake2_hash::Blake2sHash, blake2_merkle::Blake2sMerkleChannel},
        verifier::{verify, VerificationError},
    },
    prover::{
        backend::simd::{column::BaseColumn, SimdBackend},
        poly::circle::{CircleEvaluation, PolyOps},
        prove,
        vcs::prover::MerkleProver,
        CommitmentSchemeProver, CommitmentTreeProver, Poly,
    },
};
use stwo_constraint_framework::{
    preprocessed_columns::PreProcessedColumnId, TraceLocationAllocator,
};
use tracing::{debug, info, span, Level};
use utils::simd::aligned_vec_from_slice;

use crate::{
    components::{gen_interaction_trace, gen_trace},
    preprocessed::PreProcessedTrace,
    relations::Relations,
    sha256::pad_message,
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

/// Serializable preprocessed data for SHA-256 proving.
///
/// This struct caches the expensive preprocessing computation including:
/// - Extended polynomial evaluations (after interpolation and blowup)
/// - Merkle tree layers (the commitment structure)
/// - Column IDs for trace allocation
///
/// The data can be serialized to disk and reused across multiple proofs,
/// avoiding the need to regenerate the preprocessing each time.
#[derive(Clone, Serialize, Deserialize)]
pub struct PreprocessedData {
    /// The log_size this preprocessing was generated for.
    /// In non-dynamic mode, preprocessing for log_size=21 covers all log_size <= 21.
    pub log_size: u32,
    /// Column IDs for trace allocation
    pub ids: Vec<String>,
    /// Domain log sizes for each extended evaluation column
    pub domain_log_sizes: Vec<u32>,
    /// Extended polynomial evaluations (after blowup) - raw u32 data per column.
    /// Each inner Vec contains the flattened PackedBaseField data as u32 values.
    pub extended_evals: Vec<Vec<u32>>,
    /// Merkle tree layers (hashes) - Vec<Blake2sHash> per layer.
    /// First layer is the root (single hash), last layer is the largest.
    pub merkle_layers: Vec<Vec<Blake2sHash>>,
}

impl PreprocessedData {
    /// Get the preprocessed column IDs as PreProcessedColumnId objects.
    pub fn column_ids(&self) -> Vec<PreProcessedColumnId> {
        self.ids
            .iter()
            .map(|id| PreProcessedColumnId { id: id.clone() })
            .collect()
    }

    /// Reconstruct the CommitmentTreeProver from the cached data.
    ///
    /// This converts the serialized data back into the stwo types needed for proving.
    /// Uses 64-byte aligned allocation for SIMD compatibility.
    pub fn to_commitment_tree(
        &self,
    ) -> (
        Vec<Poly<SimdBackend>>,
        MerkleProver<SimdBackend, Blake2sMerkleHasher>,
    ) {
        // Reconstruct polynomials from extended evaluations
        let polynomials: Vec<Poly<SimdBackend>> = self
            .extended_evals
            .iter()
            .zip(self.domain_log_sizes.iter())
            .map(|(data, &domain_log_size)| {
                // Allocate 64-byte aligned memory and copy data
                let aligned_data: Vec<u32> = aligned_vec_from_slice(data);

                // Convert to Vec<PackedBaseField> using bytemuck
                // Safety: aligned_vec_from_slice ensures 64-byte alignment
                let packed_data: Vec<_> = bytemuck::cast_slice(&aligned_data).to_vec();

                let values = BaseColumn::from_simd(packed_data);
                let domain = CanonicCoset::new(domain_log_size).circle_domain();
                let evals = CircleEvaluation::new(domain, values);

                Poly::new(None, evals)
            })
            .collect();

        // Reconstruct MerkleProver from layers
        let merkle_prover = MerkleProver {
            layers: self.merkle_layers.clone(),
        };

        (polynomials, merkle_prover)
    }
}

/// Generate preprocessed data for SHA-256 proving.
///
/// This function performs the expensive preprocessing computation once, including:
/// - Generating the preprocessed trace columns
/// - Interpolating and extending polynomials (with blowup factor)
/// - Building the Merkle tree commitment
///
/// The result can be serialized to disk and reused for multiple proofs with the same
/// or smaller log_size (in non-dynamic mode).
///
/// # Arguments
/// * `log_size` - The log2 of the trace size to preprocess for
/// * `config` - The PCS configuration (needed for blowup factor)
///
/// # Returns
/// A `PreprocessedData` struct that can be serialized and reused.
pub fn preprocess_sha256(log_size: u32, config: PcsConfig) -> PreprocessedData {
    // 1. Generate PreProcessedTrace (raw columns + ids)
    let span = span!(Level::INFO, "Preprocess").entered();

    let span_1 = span!(Level::INFO, "Generate trace").entered();
    let preprocessed_trace = PreProcessedTrace::new(log_size);
    span_1.exit();

    // 2. Compute twiddles
    let span_2 = span!(Level::INFO, "Precompute twiddles").entered();
    #[cfg(feature = "dynamic-preprocessed-shape")]
    let twiddles_log_size = preprocessed_log_size(log_size);
    #[cfg(not(feature = "dynamic-preprocessed-shape"))]
    let twiddles_log_size = TWIDDLES_LOG_SIZE;

    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(twiddles_log_size + config.fri_config.log_blowup_factor + 2)
            .circle_domain()
            .half_coset,
    );
    span_2.exit();

    // 3. Create commitment scheme and commit
    let span_3 = span!(Level::INFO, "Commit").entered();
    let channel = &mut Blake2sChannel::default();
    let mut commitment_scheme =
        CommitmentSchemeProver::<_, Blake2sMerkleChannel>::new(config, &twiddles);

    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.trace);
    tree_builder.commit(channel);
    span_3.exit();

    // 4. Extract the committed tree data
    let span_4 = span!(Level::INFO, "Extract data").entered();
    let tree = &commitment_scheme.trees[0];

    // Extract column IDs
    let ids: Vec<String> = preprocessed_trace
        .ids
        .iter()
        .map(|id| id.id.clone())
        .collect();

    // Extract domain log sizes and extended evaluations
    let domain_log_sizes: Vec<u32> = tree
        .polynomials
        .iter()
        .map(|poly| poly.evals.domain.log_size())
        .collect();

    let extended_evals: Vec<Vec<u32>> = tree
        .polynomials
        .iter()
        .map(|poly| {
            // Cast PackedBaseField data to u32 slice (zero-copy)
            let packed_slice: &[u32] = bytemuck::cast_slice(&poly.evals.values.data);
            packed_slice.to_vec()
        })
        .collect();

    // Extract Merkle tree layers
    let merkle_layers: Vec<Vec<Blake2sHash>> = tree.commitment.layers.clone();

    span_4.exit();
    span.exit();

    PreprocessedData {
        log_size,
        ids,
        domain_log_sizes,
        extended_evals,
        merkle_layers,
    }
}

/// Prove SHA-256 hash computation for the given input.
///
/// # Arguments
/// * `input` - The message bytes to hash
/// * `config` - The PCS configuration for the proof
/// * `preprocessed` - The preprocessed data (from `preprocess_sha256`)
///
/// # Returns
/// A tuple of (proof, log_size, claimed_sum) where log_size is the computed trace size.
pub fn prove_sha256(
    input: &[u8],
    config: PcsConfig,
    preprocessed: &PreprocessedData,
) -> (StarkProof<Blake2sMerkleHasher>, u32, components::ClaimedSum) {
    // Pad the input message according to SHA-256 rules
    let chunks = pad_message(input);

    // Generate trace and get log_size
    let span = span!(Level::INFO, "Trace").entered();
    let (log_size, trace, lookup_data) = gen_trace(&chunks);
    span.exit();

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

    // Preprocessed trace - reconstruct from cached data.
    let span = span!(Level::INFO, "Load preprocessed").entered();
    let (polynomials, merkle_prover) = preprocessed.to_commitment_tree();

    // Get the root before pushing (for channel mixing)
    let root = merkle_prover.layers[0][0];

    // Push the commitment tree directly (skips interpolation, extension, and Merkle tree building)
    commitment_scheme.trees.push(CommitmentTreeProver {
        polynomials,
        commitment: merkle_prover,
    });

    // Mix the root into the channel (same as what commit() does)
    Blake2sMerkleChannel::mix_root(channel, root);
    span.exit();

    // Commit trace.
    let span = span!(Level::INFO, "Commit trace").entered();
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

    claimed_sum.mix_into(channel);

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
            max_len
        })
    );

    // Prove constraints.
    let span = span!(Level::INFO, "Prove").entered();
    let preprocessed_ids = preprocessed.column_ids();
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preprocessed_columns(&preprocessed_ids);
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

    (proof.expect("proof should be valid"), log_size, claimed_sum)
}

/// Verify a SHA256 proof.
///
/// # Arguments
/// * `proof` - The STARK proof to verify
/// * `log_size` - The log2 of the trace size (returned by `prove_sha256`)
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

    claimed_sum.mix_into(channel);

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
    use std::time::Instant;

    use super::*;

    /// Default log_size for preprocessing in tests.
    /// In non-dynamic mode, this value doesn't matter as preprocessing uses MAX_PREPROCESSED_LOG_SIZE.
    /// In dynamic mode, we use a small value suitable for test inputs.
    const TEST_PREPROCESS_LOG_SIZE: u32 = 4;

    #[test_log::test]
    fn test_prove_sha256_hello_world() {
        print_enabled_features();

        let input = b"hello world";
        info!("Testing prove with input: {:?}", std::str::from_utf8(input));

        let config = PcsConfig::default();

        // Preprocess (in non-dynamic mode, this works for any log_size <= MAX_PREPROCESSED_LOG_SIZE)
        let preprocessed = preprocess_sha256(TEST_PREPROCESS_LOG_SIZE, config);

        let start = Instant::now();
        let (_proof, log_size, _claimed_sum) = prove_sha256(input, config, &preprocessed);
        info!(
            "Proof generated in {:?}, log_size={}",
            start.elapsed(),
            log_size
        );
    }

    #[test_log::test]
    fn test_prove_verify_roundtrip() {
        print_enabled_features();

        let input = b"hello world";
        info!(
            "Testing prove/verify roundtrip with input: {:?}",
            std::str::from_utf8(input)
        );

        let config = PcsConfig::default();

        // Preprocess
        let preprocessed = preprocess_sha256(TEST_PREPROCESS_LOG_SIZE, config);

        // Prove
        let (proof, log_size, claimed_sum) = prove_sha256(input, config, &preprocessed);
        info!("Proof generated successfully, log_size={}", log_size);

        // Verify
        let result = verify_sha256(proof, log_size, &claimed_sum);
        assert!(result.is_ok(), "Verification failed: {:?}", result.err());
        info!("Proof verified successfully");
    }

    #[test_log::test]
    fn test_prove_verify_empty_input() {
        print_enabled_features();

        let input = b"";
        info!("Testing prove/verify with empty input");

        let config = PcsConfig::default();

        // Preprocess
        let preprocessed = preprocess_sha256(TEST_PREPROCESS_LOG_SIZE, config);

        // Prove
        let (proof, log_size, claimed_sum) = prove_sha256(input, config, &preprocessed);
        info!("Proof generated successfully, log_size={}", log_size);

        // Verify
        let result = verify_sha256(proof, log_size, &claimed_sum);
        assert!(result.is_ok(), "Verification failed: {:?}", result.err());
        info!("Proof verified successfully");
    }

    #[test_log::test]
    fn test_prove_verify_long_input() {
        print_enabled_features();

        // Create input that requires 2 chunks (100 bytes)
        let input = [0xab_u8; 100];
        info!("Testing prove/verify with {} byte input", input.len());

        let config = PcsConfig::default();

        // Preprocess
        let preprocessed = preprocess_sha256(TEST_PREPROCESS_LOG_SIZE, config);

        // Prove
        let (proof, log_size, claimed_sum) = prove_sha256(&input, config, &preprocessed);
        info!("Proof generated successfully, log_size={}", log_size);

        // Verify
        let result = verify_sha256(proof, log_size, &claimed_sum);
        assert!(result.is_ok(), "Verification failed: {:?}", result.err());
        info!("Proof verified successfully");
    }
}
