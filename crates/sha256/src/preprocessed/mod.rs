//! A collection of preprocessed columns, whose values are publicly acknowledged, and independent of
//! the proof.
//!
//! They are similar to regular components but are entirely known by the verifier.
use itertools::Itertools;
use stwo::{
    core::fields::m31::BaseField,
    prover::{
        backend::simd::{m31::LOG_N_LANES, SimdBackend},
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

pub struct PreProcessedTrace {
    pub trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    pub ids: Vec<PreProcessedColumnId>,
}

impl PreProcessedTrace {
    #[allow(clippy::cognitive_complexity)]
    pub fn new(log_size: u32) -> Self {
        let mut trace = Vec::new();
        let mut ids = Vec::new();
        debug_assert!(log_size >= LOG_N_LANES);
        let chunk_size = 1 << (log_size - LOG_N_LANES);

        // Helper macro to process each module once
        macro_rules! collect_columns {
            ($mod:ident, $($id_ty:ident),+) => {{
                let module_cols = $mod::gen_column_simd();

                // Collect all ids from all provided id types
                let mut module_ids = Vec::new();
                $(
                    module_ids.extend($mod::$id_ty::to_ids(None));
                )+

                for (id_base, col) in module_ids.into_iter().zip_eq(module_cols.into_iter()) {
                    for (suffix, chunk) in col.chunks(chunk_size).enumerate() {
                        trace.push(circle_evaluation_u32x16!(chunk));
                        ids.push(PreProcessedColumnId {
                            id: format!("{}_{}", id_base.id, suffix),
                        });
                    }
                }
            }};
        }

        // Generate everything in lockstep
        collect_columns!(big_sigma_0, BigSigma0I0I1Columns, BigSigma0O2Columns);
        collect_columns!(
            big_sigma_1,
            BigSigma1I0Columns,
            BigSigma1I1Columns,
            BigSigma1O2Columns
        );
        collect_columns!(ch_left, ChLeftI0Columns, ChLeftI1Columns);
        collect_columns!(ch_right, ChRightI0Columns, ChRightI1Columns);
        collect_columns!(
            maj,
            MajI0LI1HColumns,
            MajI0H0I1L0Columns,
            MajI0H1I1L1Columns
        );
        collect_columns!(range_check_add, RangeCheckAddColumns);
        collect_columns!(sigma_0, Sigma0I0I1Columns, Sigma0O2Columns);
        collect_columns!(sigma_1, Sigma1I0I1Columns, Sigma1O2Columns);

        Self { trace, ids }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let trace = PreProcessedTrace::new(8);
        assert!(trace.trace.iter().map(|t| t.data.len()).max().unwrap() <= 1 << (8 - LOG_N_LANES));
    }
}
