pub mod components;
pub mod preprocessed;

use crate::components::Claim;
use stwo_prover::core::prover::StarkProof;
use stwo_prover::core::vcs::ops::MerkleHasher;

pub struct Proof<const N: usize, H: MerkleHasher> {
    pub claim: Claim<N>,
    pub stark_proof: StarkProof<H>,
}
