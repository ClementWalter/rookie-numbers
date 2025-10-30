#![allow(non_camel_case_types)]
#![feature(portable_simd, array_chunks, iter_array_chunks)]
pub mod components;
pub mod macros;
pub mod partitions;
pub mod preprocessed;
pub mod relations;
pub mod sha256;

use num_traits::Zero;
use stwo::{
    core::{
        channel::Blake2sChannel,
        fields::qm31::SecureField,
        pcs::PcsConfig,
        poly::circle::CanonicCoset,
        proof::StarkProof,
        vcs::blake2_merkle::{Blake2sMerkleChannel, Blake2sMerkleHasher},
    },
    prover::{backend::simd::SimdBackend, poly::circle::PolyOps, prove, CommitmentSchemeProver},
};
use stwo_constraint_framework::TraceLocationAllocator;
use tracing::{info, span, Level};

use crate::{
    components::{gen_interaction_trace, gen_trace},
    preprocessed::PreProcessedTrace,
    relations::Relations,
};

pub fn prove_sha256(log_size: u32, config: PcsConfig) -> StarkProof<Blake2sMerkleHasher> {
    // Precompute twiddles.
    let span = span!(Level::INFO, "Precompute twiddles").entered();
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(log_size + config.fri_config.log_blowup_factor + 2)
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

    info!(
        "Columns count: {:?}",
        commitment_scheme
            .trees
            .as_ref()
            .map(|tree| tree.evaluations.len())
    );
    info!(
        "Columns length: {:?}",
        commitment_scheme.trees.as_ref().map(|tree| {
            let max_len = tree
                .evaluations
                .iter()
                .map(|eval| eval.values.length.ilog2())
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

    #[cfg(debug_assertions)]
    println!(
        "Trace log degree bounds: {:?}",
        components.trace_log_degree_bounds()
    );

    if claimed_sum.scheduling + claimed_sum.compression + claimed_sum.preprocessed.sum()
        != SecureField::zero()
    {
        #[cfg(debug_assertions)]
        tracing::info!(
            "Relation summary: {:?}",
            components.track_relations(&commitment_scheme)
        );
        panic!("Relation summary is not zero");
    }

    let proof = prove(&components.provers(), channel, commitment_scheme);
    if let Err(e) = proof {
        panic!("Proof error: {e:?}");
    }
    span.exit();

    proof.unwrap()
}

#[cfg(test)]
mod tests {
    use std::{env, time::Instant};

    use peak_alloc::PeakAlloc;
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use tracing::info;

    use super::*;

    #[global_allocator]
    static PEAK_ALLOC: PeakAlloc = PeakAlloc;

    #[test_log::test]
    fn test_prove_sha256() {
        #[cfg(feature = "parallel")]
        info!("Stwo Parallel");
        #[cfg(not(feature = "parallel"))]
        info!("Stwo Non-parallel");

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

        let peak_bytes = PEAK_ALLOC.peak_usage_as_mb();
        info!("Peak memory: {} MB", peak_bytes);
    }
}
