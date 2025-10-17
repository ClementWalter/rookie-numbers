use std::simd::u32x16;

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
    scheduling_lookup_data: &[Vec<u32x16>],
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let _sigma_0 = sigma_0_i0_i1::witness::gen_trace(scheduling_lookup_data);

    vec![]
}

pub fn gen_interaction_trace(
    trace: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    relations: &Relations,
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let _sigma_0 = sigma_0_i0_i1::witness::gen_interaction_trace(trace, relations);

    vec![]
}
