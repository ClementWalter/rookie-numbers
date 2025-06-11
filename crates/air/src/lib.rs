#![feature(raw_slice_split)]

pub mod components;
pub mod preprocessed;
pub mod relations;

use crate::components::{Claim, InteractionClaim};
use stwo_prover::core::prover::StarkProof;
use stwo_prover::core::vcs::ops::MerkleHasher;

pub struct Proof<const N: usize, H: MerkleHasher> {
    pub claim: Claim<N>,
    pub interaction_claim: InteractionClaim,
    pub stark_proof: StarkProof<H>,
    pub interaction_pow: u64,
}
