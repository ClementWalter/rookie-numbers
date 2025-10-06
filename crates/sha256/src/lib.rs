#![allow(non_camel_case_types)]
pub mod relations;

use num_traits::Zero;
use stwo_prover::constraint_framework::EvalAtRow;
use stwo_prover::constraint_framework::FrameworkComponent;
use stwo_prover::constraint_framework::FrameworkEval;
use stwo_prover::constraint_framework::TraceLocationAllocator;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::channel::Blake2sChannel;
use stwo_prover::core::fields::m31::BaseField;
use stwo_prover::core::fields::qm31::SecureField;
use stwo_prover::core::pcs::{CommitmentSchemeProver, PcsConfig};
use stwo_prover::core::poly::circle::{CanonicCoset, CircleEvaluation, PolyOps};
use stwo_prover::core::poly::BitReversedOrder;
use stwo_prover::core::prover::{prove, StarkProof};
use stwo_prover::core::vcs::blake2_merkle::{Blake2sMerkleChannel, Blake2sMerkleHasher};
use stwo_prover::core::ColumnVec;

use crate::relations::{LookupData, Relations};

use tracing::{info, span, Level};

const CHUNK_SIZE: usize = 32; // 16 u32 = 32 u16
const H_SIZE: usize = 16; // 8 u32 = 16 u16

pub type Component = FrameworkComponent<Eval>;

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
        self.log_size() + 1
    }
    fn evaluate<E: EvalAtRow>(&self, eval: E) -> E {
        eval
    }
}

pub fn gen_trace(
    log_size: u32,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    LookupData,
) {
    let trace = ColumnVec::new();
    let lookup_data = LookupData::new(log_size);

    (trace, lookup_data)
}

pub fn gen_interaction_trace(
    _log_size: u32,
    _lookup_data: LookupData,
    _relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    SecureField,
) {
    let trace = ColumnVec::new();
    let claimed_sum = SecureField::zero();
    (trace, claimed_sum)
}

pub fn prove_sha256(
    log_size: u32,
    config: PcsConfig,
) -> (Component, StarkProof<Blake2sMerkleHasher>) {
    // Precompute twiddles.
    let span = span!(Level::INFO, "Precompute twiddles").entered();
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(log_size + 1 + config.fri_config.log_blowup_factor)
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
    let mut tree_builder = commitment_scheme.tree_builder();
    let constant_trace = vec![];
    tree_builder.extend_evals(constant_trace);
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
    let (trace, claimed_sum) = gen_interaction_trace(log_size, lookup_data, &relations);
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);
    span.exit();

    // Prove constraints.
    let component = Component::new(
        &mut TraceLocationAllocator::default(),
        Eval {
            log_size,
            relations,
        },
        claimed_sum,
    );
    info!("Poseidon component info:\n{}", component);
    let proof = prove(&[&component], channel, commitment_scheme).unwrap();

    (component, proof)
}

#[cfg(test)]
mod tests {}
