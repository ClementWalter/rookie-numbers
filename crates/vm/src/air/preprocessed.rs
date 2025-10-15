//! A collection of preprocessed columns, whose values are publicly acknowledged, and independent of
//! the proof.
use stwo::{
    core::fields::m31::BaseField,
    prover::{
        backend::simd::SimdBackend,
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId;

pub trait PreProcessedColumn {
    fn log_size(&self) -> u32;
    fn id(&self) -> PreProcessedColumnId;
    fn gen_column_simd(&self) -> CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>;
}

pub struct PreProcessedTrace {
    columns: Vec<Box<dyn PreProcessedColumn>>,
}
impl PreProcessedTrace {
    pub fn new(columns: Vec<Box<dyn PreProcessedColumn>>) -> Self {
        Self { columns }
    }

    pub fn log_sizes(&self) -> Vec<u32> {
        self.columns.iter().map(|c| c.log_size()).collect()
    }

    pub fn gen_trace(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        self.columns.iter().map(|c| c.gen_column_simd()).collect()
    }

    pub fn ids(&self) -> Vec<PreProcessedColumnId> {
        self.columns.iter().map(|c| c.id()).collect()
    }
}

impl Default for PreProcessedTrace {
    fn default() -> Self {
        let columns = vec![];
        Self::new(columns)
    }
}
