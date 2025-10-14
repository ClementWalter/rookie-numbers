#![allow(non_camel_case_types)]
#![feature(portable_simd)]
pub mod components;
pub mod partitions;
pub mod preprocessed;
pub mod relations;
pub mod sha256;

use stwo_prover::constraint_framework::EvalAtRow;
use stwo_prover::constraint_framework::FrameworkEval;
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
const LOG_EXPAND: u32 = 1;

#[derive(Clone)]
pub struct Eval {
    pub log_size: u32,
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        self.log_size
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        self.log_size() + LOG_EXPAND
    }
    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        eval_sha256_constraints(&mut eval, &self.relations);
        eval
    }
}

pub const fn eval_sha256_constraints<E: EvalAtRow>(_eval: &mut E, _lookup_elements: &Relations) {}

pub fn prove_sha256(log_size: u32, config: PcsConfig) -> StarkProof<Blake2sMerkleHasher> {
    // Precompute twiddles.
    let span = span!(Level::INFO, "Precompute twiddles").entered();
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(log_size + LOG_EXPAND + config.fri_config.log_blowup_factor)
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
    let scheduling_component = components::scheduling::air::Component::new(
        &mut TraceLocationAllocator::new_with_preproccessed_columns(&preprocessed_trace.ids()),
        components::scheduling::air::Eval {
            log_size,
            relations: relations.clone(),
        },
        claimed_sum.scheduling,
    );

    let compression_component = components::compression::air::Component::new(
        &mut TraceLocationAllocator::new_with_preproccessed_columns(&preprocessed_trace.ids()),
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
