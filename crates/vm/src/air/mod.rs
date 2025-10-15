pub mod components;
pub mod preprocessed;
pub mod relations;

use components::{Claim, InteractionClaim};
use serde::{Deserialize, Serialize};
use stwo::core::{proof::StarkProof, vcs::MerkleHasher};

#[derive(Serialize, Deserialize)]
pub struct Proof<const N: usize, H: MerkleHasher> {
    pub claim: Claim<N>,
    pub interaction_claim: InteractionClaim<N>,
    pub stark_proof: StarkProof<H>,
    pub interaction_pow: u64,
}
