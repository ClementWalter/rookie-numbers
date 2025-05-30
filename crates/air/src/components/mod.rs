pub mod multiple_constraints;
pub mod single_constraint;

pub use stwo_air_utils::trace::component_trace::ComponentTrace;
pub use stwo_air_utils_derive::{IterMut, ParIterMut, Uninitialized};
use stwo_prover::core::air::Component as ComponentVerifier;
pub use stwo_prover::core::backend::simd::m31::PackedM31;
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

pub struct Claim<const N: usize> {
    pub single_constraint: single_constraint::Claim<N>,
    pub multiple_constraints: multiple_constraints::Claim<N>,
}

impl<const N: usize> Claim<N> {
    pub fn new(log_size: u32) -> Self {
        Self {
            single_constraint: single_constraint::Claim::new(log_size),
            multiple_constraints: multiple_constraints::Claim::new(log_size),
        }
    }

    pub fn log_sizes(&self) -> TreeVec<Vec<u32>> {
        let trees = vec![
            self.single_constraint.log_sizes(),
            self.multiple_constraints.log_sizes(),
        ];
        TreeVec::concat_cols(trees.into_iter())
    }

    pub fn mix_into(&self, channel: &mut impl Channel) {
        self.single_constraint.mix_into(channel);
        self.multiple_constraints.mix_into(channel);
    }

    pub fn write_trace<MC: MerkleChannel>(
        &self,
    ) -> impl IntoIterator<Item = CircleEvaluation<SimdBackend, M31, BitReversedOrder>>
    where
        SimdBackend: BackendForChannel<MC>,
    {
        let single_trace = self.single_constraint.write_trace().to_evals();
        let multiple_trace = self.multiple_constraints.write_trace().to_evals();
        [single_trace, multiple_trace].into_iter().flatten()
    }
}

pub struct Components<const N: usize> {
    pub single_constraint: single_constraint::Component<N>,
    pub multiple_constraints: multiple_constraints::Component<N>,
}

impl<const N: usize> Components<N> {
    pub fn new(location_allocator: &mut TraceLocationAllocator, claim: &Claim<N>) -> Self {
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
        }
    }

    pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
        vec![&self.single_constraint, &self.multiple_constraints]
    }

    pub fn verifiers(&self) -> Vec<&dyn ComponentVerifier> {
        vec![&self.single_constraint, &self.multiple_constraints]
    }
}
