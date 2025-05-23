pub mod components;
pub mod preprocessed;

use crate::components::single_constraint;
use stwo_prover::core::prover::StarkProof;
use stwo_prover::core::vcs::ops::MerkleHasher;

pub struct Proof<const N: usize, H: MerkleHasher> {
    pub claim: single_constraint::Claim<N>,
    pub stark_proof: StarkProof<H>,
}
