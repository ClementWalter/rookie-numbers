#![allow(non_camel_case_types)]
use stwo_prover::relation;
use stwo_prover::core::channel::{ Channel, MerkleChannel };
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::backend::BackendForChannel;
use stwo_prover::core::pcs::TreeBuilder;
use crate::components::single_constraint;

relation!(SimpleRelation, 1);

/// Logup security is defined by the `QM31` space (~124 bits) + `INTERACTION_POW_BITS` -
/// log2(number of relation terms).
/// E.g. assuming a 100-bit security target, the witness may contain up to
/// 1 << (24 + INTERACTION_POW_BITS) relation terms.
pub const INTERACTION_POW_BITS: u32 = 24;

pub struct InteractionElements {
    pub simple_relation: SimpleRelation,
}
impl InteractionElements {
    pub fn draw(channel: &mut impl Channel) -> InteractionElements {
        InteractionElements {
            simple_relation: SimpleRelation::draw(channel),
        }
    }
}

pub struct InteractionClaim {
    pub single_constraint: single_constraint::InteractionClaim,
}

impl InteractionClaim {
    pub fn mix_into(&self, channel: &mut impl Channel) {
        self.single_constraint.mix_into(channel);
    }
}

pub struct InteractionClaimGenerator {
    pub single_constraint: single_constraint::InteractionClaimGenerator,
}

impl InteractionClaimGenerator {
    pub fn write_interaction_trace<MC: MerkleChannel>(
        &self,
        _tree_builder: &mut TreeBuilder<SimdBackend, MC>,
        _interaction_elements: &InteractionElements
    ) -> InteractionClaim
        where SimdBackend: BackendForChannel<MC>
    {
        InteractionClaim {
            single_constraint: self.single_constraint.write_interaction_trace(),
        }
    }
}
