use stwo::{
    core::{fields::m31::BaseField, ColumnVec},
    prover::{
        backend::simd::SimdBackend,
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};

use crate::relations::Relations;

mod sigma_0_i0_i1;

pub fn gen_trace(
    scheduling: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let _sigma_0 = sigma_0_i0_i1::witness::gen_trace(scheduling);

    vec![]
}

pub fn gen_interaction_trace(
    scheduling: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    relations: &Relations,
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let _sigma_0 = sigma_0_i0_i1::witness::gen_interaction_trace(scheduling, relations);

    vec![]
}
