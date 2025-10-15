#![allow(non_camel_case_types)]
#![feature(portable_simd)]
pub mod components;
pub mod partitions;
pub mod preprocessed;
pub mod relations;
pub mod sha256;

use stwo_prover::constraint_framework::TraceLocationAllocator;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::channel::Blake2sChannel;
use stwo_prover::core::pcs::{CommitmentSchemeProver, PcsConfig};
use stwo_prover::core::poly::circle::{CanonicCoset, PolyOps};
use stwo_prover::core::prover::{prove, StarkProof};
use stwo_prover::core::vcs::blake2_merkle::{Blake2sMerkleChannel, Blake2sMerkleHasher};

use crate::components::gen_interaction_trace;
use crate::components::gen_trace;
use crate::preprocessed::PreProcessedTrace;
use crate::relations::Relations;

use tracing::{info, span, Level};

const CHUNK_SIZE: usize = 32; // 16 u32 = 32 u16

pub fn prove_sha256(log_size: u32, config: PcsConfig) -> StarkProof<Blake2sMerkleHasher> {
    // Precompute twiddles.
    let span = span!(Level::INFO, "Precompute twiddles").entered();
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(21 + 1 + config.fri_config.log_blowup_factor)
            .circle_domain()
            .half_coset,
    );
    span.exit();

    // Setup protocol.
    let channel = &mut Blake2sChannel::default();
    let mut commitment_scheme =
        CommitmentSchemeProver::<_, Blake2sMerkleChannel>::new(config, &twiddles);

    // Preprocessed trace.
    let preprocessed_trace = PreProcessedTrace::default();
    let span = span!(Level::INFO, "Constant").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.gen_trace());
    tree_builder.commit(channel);
    span.exit();

    // Trace.
    let span = span!(Level::INFO, "Trace").entered();
    let (trace, lookup_data) = gen_trace(log_size);
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);
    span.exit();

    // Draw lookup elements.
    let relations = Relations::draw(channel);

    // Interaction trace.
    let span = span!(Level::INFO, "Interaction").entered();
    let (trace, claimed_sum) = gen_interaction_trace(lookup_data, &relations);
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);
    span.exit();

    // Prove constraints.
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preproccessed_columns(&preprocessed_trace.ids());
    let scheduling_component = components::scheduling::air::Component::new(
        trace_allocator,
        components::scheduling::air::Eval {
            log_size,
            relations: relations.clone(),
        },
        claimed_sum.scheduling,
    );

    let compression_component = components::compression::air::Component::new(
        trace_allocator,
        components::compression::air::Eval {
            log_size,
            relations,
        },
        claimed_sum.compression,
    );

    info!("Scheduling component info:\n{}", scheduling_component);
    prove(
        &[&scheduling_component, &compression_component],
        channel,
        commitment_scheme,
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prove_sha256() {
        use std::time::Instant;

        let log_size = 18;
        let start = Instant::now();
        prove_sha256(log_size, PcsConfig::default());
        let duration = start.elapsed();
        println!("Time spent in prove_sha256: {:?}", duration);
        println!("Log size: {}", log_size);
    }
}
