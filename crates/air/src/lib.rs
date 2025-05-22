pub mod preprocessed;
pub mod components;

use stwo_prover::core::prover::StarkProof;
use stwo_prover::core::vcs::ops::MerkleHasher;
use crate::components::single_constraint;

pub struct Proof<H: MerkleHasher> {
    pub claim: single_constraint::Claim,
    pub stark_proof: StarkProof<H>,
}
