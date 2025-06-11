pub mod multiple_constraints;
pub mod single_constraint;
pub mod single_constraint_with_relation;

use num_traits::Zero;
pub use stwo_air_utils::trace::component_trace::ComponentTrace;
pub use stwo_air_utils_derive::{IterMut, ParIterMut, Uninitialized};
pub use stwo_prover::core::backend::simd::m31::PackedM31;
use stwo_prover::core::{air::Component as ComponentVerifier, fields::qm31::QM31};
use stwo_prover::{
    constraint_framework::TraceLocationAllocator,
    core::{
        air::ComponentProver,
        backend::{simd::SimdBackend, BackendForChannel},
        channel::{Channel, MerkleChannel},
        fields::{m31::M31, qm31::SecureField},
        pcs::TreeVec,
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};

use crate::relations;

pub struct Claim<const N: usize> {
    pub single_constraint: single_constraint::Claim<N>,
    pub multiple_constraints: multiple_constraints::Claim<N>,
    pub single_constraint_with_relation: single_constraint_with_relation::Claim<N>,
}

pub struct LookupData {
    pub single_constraint_with_relation: single_constraint_with_relation::LookupData,
}

pub struct Relations {
    pub memory: relations::Memory,
}

pub struct InteractionClaim {
    pub single_constraint_with_relation: single_constraint_with_relation::InteractionClaim,
}

impl<const N: usize> Claim<N> {
    pub fn new(log_size: u32) -> Self {
        Self {
            single_constraint: single_constraint::Claim::new(log_size),
            multiple_constraints: multiple_constraints::Claim::new(log_size),
            single_constraint_with_relation: single_constraint_with_relation::Claim::new(log_size),
        }
    }

    pub fn log_sizes(&self) -> TreeVec<Vec<u32>> {
        let trees = vec![
            self.single_constraint.log_sizes(),
            self.multiple_constraints.log_sizes(),
            self.single_constraint_with_relation.log_sizes(),
        ];
        TreeVec::concat_cols(trees.into_iter())
    }

    pub fn mix_into(&self, channel: &mut impl Channel) {
        self.single_constraint.mix_into(channel);
        self.multiple_constraints.mix_into(channel);
        self.single_constraint_with_relation.mix_into(channel);
    }

    pub fn write_trace<MC: MerkleChannel>(
        &self,
    ) -> (
        impl IntoIterator<Item = CircleEvaluation<SimdBackend, M31, BitReversedOrder>>,
        LookupData,
    )
    where
        SimdBackend: BackendForChannel<MC>,
    {
        let single_trace = self.single_constraint.write_trace();
        let multiple_trace = self.multiple_constraints.write_trace();
        let (single_constraint_with_relation_trace, lookup_data) =
            self.single_constraint_with_relation.write_trace();
        (
            [
                single_trace.to_evals(),
                multiple_trace.to_evals(),
                single_constraint_with_relation_trace.to_evals(),
            ]
            .into_iter()
            .flatten(),
            LookupData {
                single_constraint_with_relation: lookup_data,
            },
        )
    }
}

impl InteractionClaim {
    pub fn write_interaction_trace(
        relations: &Relations,
        lookup_data: &LookupData,
    ) -> (
        impl IntoIterator<Item = CircleEvaluation<SimdBackend, M31, BitReversedOrder>>,
        InteractionClaim,
    ) {
        let (
            single_constraint_with_relation_trace,
            single_constraint_with_relation_trace_claimed_sum,
        ) = single_constraint_with_relation::InteractionClaim::write_interaction_trace(
            &relations.memory,
            &lookup_data.single_constraint_with_relation,
        );

        (
            [single_constraint_with_relation_trace]
                .into_iter()
                .flatten(),
            InteractionClaim {
                single_constraint_with_relation: single_constraint_with_relation_trace_claimed_sum,
            },
        )
    }

    pub fn claimed_sum(&self) -> SecureField {
        let mut sum = QM31::zero();
        sum += self.single_constraint_with_relation.claimed_sum;
        sum
    }

    pub fn mix_into(&self, channel: &mut impl Channel) {
        self.single_constraint_with_relation.mix_into(channel);
    }
}

impl Relations {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            memory: relations::Memory::draw(channel),
        }
    }
}

pub struct Components<const N: usize> {
    pub single_constraint: single_constraint::Component<N>,
    pub multiple_constraints: multiple_constraints::Component<N>,
    pub single_constraint_with_relation: single_constraint_with_relation::Component<N>,
}

impl<const N: usize> Components<N> {
    pub fn new(
        location_allocator: &mut TraceLocationAllocator,
        claim: &Claim<N>,
        interaction_claim: &InteractionClaim,
        relations: &Relations,
    ) -> Self {
        Self {
            single_constraint: single_constraint::Component::new(
                location_allocator,
                single_constraint::Eval {
                    claim: claim.single_constraint,
                },
                SecureField::default(),
            ),
            multiple_constraints: multiple_constraints::Component::new(
                location_allocator,
                multiple_constraints::Eval {
                    claim: claim.multiple_constraints,
                },
                SecureField::default(),
            ),
            single_constraint_with_relation: single_constraint_with_relation::Component::new(
                location_allocator,
                single_constraint_with_relation::Eval {
                    claim: claim.single_constraint_with_relation,
                    memory: relations.memory.clone(),
                },
                interaction_claim
                    .single_constraint_with_relation
                    .claimed_sum,
            ),
        }
    }

    pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
        vec![
            &self.single_constraint,
            &self.multiple_constraints,
            &self.single_constraint_with_relation,
        ]
    }

    pub fn verifiers(&self) -> Vec<&dyn ComponentVerifier> {
        vec![
            &self.single_constraint,
            &self.multiple_constraints,
            &self.single_constraint_with_relation,
        ]
    }
}
