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
use utils::circle_evaluation_u32x16;

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
        big_sigma_0::gen_column_simd()
            .into_iter()
            .chain(big_sigma_1::gen_column_simd())
            .chain(ch_left::gen_column_simd())
            .chain(ch_right::gen_column_simd())
            .chain(maj::gen_column_simd())
            .chain(range_check_add::gen_column_simd())
            .chain(sigma_0::gen_column_simd())
            .chain(sigma_1::gen_column_simd())
            .map(|c| circle_evaluation_u32x16!(c))
            .collect()
    }

    pub fn ids(&self) -> Vec<PreProcessedColumnId> {
        let mut ids = vec![];

        ids.extend(big_sigma_0::BigSigma0I0I1Columns::to_ids(None));
        ids.extend(big_sigma_0::BigSigma0O2Columns::to_ids(None));
        ids.extend(big_sigma_1::BigSigma1I0Columns::to_ids(None));
        ids.extend(big_sigma_1::BigSigma1I1Columns::to_ids(None));
        ids.extend(big_sigma_1::BigSigma1O2Columns::to_ids(None));
        ids.extend(ch_left::ChLeftI0Columns::to_ids(None));
        ids.extend(ch_left::ChLeftI1Columns::to_ids(None));
        ids.extend(ch_right::ChRightI0Columns::to_ids(None));
        ids.extend(ch_right::ChRightI1Columns::to_ids(None));
        ids.extend(maj::MajI0LI1HColumns::to_ids(None));
        ids.extend(maj::MajI0H0I1L0Columns::to_ids(None));
        ids.extend(maj::MajI0H1I1L1Columns::to_ids(None));
        ids.extend(range_check_add::RangeCheckAddColumns::to_ids(None));
        ids.extend(sigma_0::Sigma0I0I1Columns::to_ids(None));
        ids.extend(sigma_0::Sigma0O2Columns::to_ids(None));
        ids.extend(sigma_1::Sigma1I0I1Columns::to_ids(None));
        ids.extend(sigma_1::Sigma1O2Columns::to_ids(None));

        ids
    }
}
