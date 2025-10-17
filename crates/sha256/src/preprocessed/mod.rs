//! A collection of preprocessed columns, whose values are publicly acknowledged, and independent of
//! the proof.
//!
//! They are similar to regular components but are entirely known by the verifier.
use stwo::{
    core::fields::m31::BaseField,
    prover::{
        backend::simd::SimdBackend,
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId;

pub mod big_sigma_0;
pub mod big_sigma_1;
pub mod ch_left;
pub mod ch_right;
pub mod maj;
pub mod range_check_add;
pub mod sigma_0;
pub mod sigma_1;

pub struct PreProcessedTrace;

impl PreProcessedTrace {
    pub fn gen_trace(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        let mut traces = vec![];

        traces.extend(big_sigma_0::gen_column_simd());
        traces.extend(big_sigma_1::gen_column_simd());
        traces.extend(ch_left::gen_column_simd());
        traces.extend(ch_right::gen_column_simd());
        traces.extend(maj::gen_column_simd());
        traces.extend(range_check_add::gen_column_simd());
        traces.extend(sigma_0::gen_column_simd());
        traces.extend(sigma_1::gen_column_simd());

        traces
    }

    pub fn ids(&self) -> Vec<PreProcessedColumnId> {
        let mut ids = vec![];

        ids.extend(big_sigma_0::Columns::to_ids());
        ids.extend(big_sigma_1::Columns::to_ids());
        ids.extend(ch_left::Columns::to_ids());
        ids.extend(ch_right::Columns::to_ids());
        ids.extend(maj::Columns::to_ids());
        ids.extend(range_check_add::Columns::to_ids());
        ids.extend(sigma_0::Columns::to_ids());
        ids.extend(sigma_1::Columns::to_ids());

        ids
    }
}
