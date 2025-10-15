#![allow(non_camel_case_types)]
#![feature(portable_simd)]
pub mod components;
pub mod partitions;
pub mod preprocessed;
pub mod relations;
pub mod sha256;

use stwo::{
    core::{
        channel::Blake2sChannel,
        pcs::PcsConfig,
        poly::circle::CanonicCoset,
        proof::StarkProof,
        vcs::blake2_merkle::{Blake2sMerkleChannel, Blake2sMerkleHasher},
    },
    prover::{backend::simd::SimdBackend, poly::circle::PolyOps, prove, CommitmentSchemeProver},
};
use stwo_constraint_framework::TraceLocationAllocator;
use tracing::{span, Level};

use crate::{
    components::{gen_interaction_trace, gen_trace},
    preprocessed::PreProcessedTrace,
    relations::Relations,
};

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
    let span = span!(Level::INFO, "Prove").entered();
    let trace_allocator =
        &mut TraceLocationAllocator::new_with_preprocessed_columns(&preprocessed_trace.ids());
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

    let proof = prove(
        &[&scheduling_component, &compression_component],
        channel,
        commitment_scheme,
    )
    .unwrap();
    span.exit();

    proof
}

#[cfg(test)]
mod tests {
    use tracing::info;

    use super::*;

    #[test_log::test]
    fn test_prove_sha256() {
        let log_size = 18;
        info!("Log size: {}", log_size);
        prove_sha256(log_size, PcsConfig::default());
    }
}
