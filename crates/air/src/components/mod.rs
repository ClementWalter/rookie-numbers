pub mod single_constraint;

use stwo_prover::core::{
    backend::{ simd::SimdBackend, BackendForChannel },
    channel::{ Channel, MerkleChannel },
    pcs::TreeBuilder,
};

use crate::interactions::InteractionClaimGenerator;

pub struct Claim {
    pub single_constraint: single_constraint::Claim,
}

impl Claim {
    pub fn mix_into(&self, channel: &mut impl Channel) {
        self.single_constraint.mix_into(channel);
    }
}

pub struct ClaimGenerator {
    pub single_constraint: single_constraint::ClaimGenerator,
}

impl ClaimGenerator {
    pub fn new(size: u32) -> Self {
        Self { single_constraint: single_constraint::ClaimGenerator::new(size) }
    }

    pub fn write_trace<MC: MerkleChannel>(
        &self,
        tree_builder: &mut TreeBuilder<SimdBackend, MC>
    ) -> (Claim, InteractionClaimGenerator)
        where SimdBackend: BackendForChannel<MC>
    {
        let (single_constraint_claim, single_constraint_interaction_claim_generator) =
            self.single_constraint.write_trace(tree_builder);
        (
            Claim { single_constraint: single_constraint_claim },
            InteractionClaimGenerator {
                single_constraint: single_constraint_interaction_claim_generator,
            },
        )
    }
}
