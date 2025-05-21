use num_traits::identities::Zero;
use stwo_prover::{
    constraint_framework::{ EvalAtRow, FrameworkComponent, FrameworkEval },
    core::{
        backend::{ simd::SimdBackend, BackendForChannel },
        channel::{ Channel, MerkleChannel },
        fields::{ m31::M31, qm31::{ SecureField } },
        pcs::{ TreeBuilder, TreeVec },
    },
};
pub use stwo_air_utils::trace::component_trace::ComponentTrace;
pub use stwo_air_utils_derive::{ IterMut, ParIterMut, Uninitialized };
use rayon::iter::{ ParallelIterator };
pub use stwo_prover::core::backend::simd::m31::{ PackedM31 };

pub const N_TRACE_COLUMNS: usize = 3;

#[derive(Copy, Clone)]
pub struct Claim {
    pub log_size: u32,
}

impl Claim {
    pub fn log_sizes(&self) -> TreeVec<Vec<u32>> {
        let trace_log_sizes = vec![self.log_size; N_TRACE_COLUMNS];
        TreeVec::new(vec![vec![], trace_log_sizes, vec![]])
    }

    pub fn mix_into(&self, channel: &mut impl Channel) {
        channel.mix_u64(self.log_size as u64);
    }
}

pub struct InteractionClaim {
    pub claimed_sum: SecureField,
}

impl InteractionClaim {
    pub fn mix_into(&self, channel: &mut impl Channel) {
        channel.mix_felts(&[self.claimed_sum]);
    }
}

pub struct Eval {
    pub claim: Claim,
}

pub type Component = FrameworkComponent<Eval>;

impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        self.claim.log_size
    }

    fn max_constraint_log_degree_bound(&self) -> u32 {
        self.log_size() + 1
    }

    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        let col0 = eval.next_trace_mask();
        let col1 = eval.next_trace_mask();
        let col2 = eval.next_trace_mask();

        eval.add_constraint(col0.clone() + col1.clone() - col2.clone());

        eval
    }
}

pub struct ClaimGenerator {
    pub size: u32,
}

impl ClaimGenerator {
    pub fn new(size: u32) -> Self {
        Self { size }
    }

    #[allow(non_snake_case)]
    pub fn write_trace<MC: MerkleChannel>(
        &self,
        tree_builder: &mut TreeBuilder<SimdBackend, MC>
    ) -> (Claim, InteractionClaimGenerator)
        where SimdBackend: BackendForChannel<MC>
    {
        let mut trace = unsafe {
            ComponentTrace::<N_TRACE_COLUMNS>::uninitialized(self.size.ilog2())
        };
        let M31_0 = PackedM31::broadcast(M31::from(0));
        let M31_1 = PackedM31::broadcast(M31::from(1));
        trace.par_iter_mut().for_each(|mut row| {
            *row[0] = M31_0;
            *row[1] = M31_1;
            *row[2] = M31_1;
        });
        tree_builder.extend_evals(trace.to_evals());

        let claim = Claim { log_size: self.size.ilog2() };
        let interaction_claim_generator = InteractionClaimGenerator {};
        (claim, interaction_claim_generator)
    }
}

pub struct InteractionClaimGenerator {}

impl InteractionClaimGenerator {
    pub fn write_interaction_trace(&self) -> InteractionClaim {
        InteractionClaim {
            claimed_sum: SecureField::zero(),
        }
    }
}
