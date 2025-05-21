pub mod interactions;
pub mod preprocessed;
pub mod components;

use stwo_prover::constraint_framework::TraceLocationAllocator;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::air::ComponentProver;
use stwo_prover::core::prover::StarkProof;
use stwo_prover::core::vcs::ops::MerkleHasher;
use crate::components::{ Claim, single_constraint };
use crate::interactions::{ InteractionClaim, InteractionElements };

pub struct Components {
    pub single_constraint: single_constraint::Component,
}

impl Components {
    pub fn new(
        tree_span_provider: &mut TraceLocationAllocator,
        claim: &Claim,
        _interaction_elements: &InteractionElements,
        interaction_claim: &InteractionClaim
    ) -> Self {
        Self {
            single_constraint: single_constraint::Component::new(
                tree_span_provider,
                single_constraint::Eval { claim: claim.single_constraint },
                interaction_claim.single_constraint.claimed_sum
            ),
        }
    }

    pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
        vec![&self.single_constraint as &dyn ComponentProver<SimdBackend>]
    }
}

pub struct Proof<H: MerkleHasher> {
    pub claim: Claim,
    pub interaction_pow: u64,
    pub interaction_claim: InteractionClaim,
    pub stark_proof: StarkProof<H>,
}
