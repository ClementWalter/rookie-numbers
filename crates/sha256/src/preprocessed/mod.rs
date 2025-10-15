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

pub trait PreProcessedColumn {
    fn log_size(&self) -> Vec<u32>;
    fn id(&self) -> Vec<PreProcessedColumnId>;
    fn gen_column_simd(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>;
}

pub struct PreProcessedTrace {
    columns: Vec<Box<dyn PreProcessedColumn>>,
}
impl PreProcessedTrace {
    pub fn log_sizes(&self) -> Vec<u32> {
        self.columns.iter().flat_map(|c| c.log_size()).collect()
    }

    pub fn gen_trace(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        self.columns
            .iter()
            .flat_map(|c| c.gen_column_simd())
            .collect()
    }

    pub fn ids(&self) -> Vec<PreProcessedColumnId> {
        self.columns.iter().flat_map(|c| c.id()).collect()
    }
}

impl Default for PreProcessedTrace {
    fn default() -> Self {
        let columns: Vec<Box<dyn PreProcessedColumn>> = vec![
            Box::new(big_sigma_0::Columns),
            Box::new(big_sigma_1::Columns),
            Box::new(ch_left::Columns),
            Box::new(ch_right::Columns),
            Box::new(maj::Columns),
            Box::new(range_check_add::Columns),
            Box::new(sigma_0::Columns),
            Box::new(sigma_1::Columns),
        ];
        Self { columns }
    }
}
