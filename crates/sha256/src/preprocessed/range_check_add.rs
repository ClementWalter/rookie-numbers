use crate::preprocessed::PreProcessedColumn;

use itertools::Itertools;
use stwo_prover::core::backend::simd::column::BaseColumn;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::fields::m31::BaseField;
use stwo_prover::core::poly::circle::{CanonicCoset, CircleEvaluation};
use stwo_prover::core::poly::BitReversedOrder;
use stwo_prover::relation;

use stwo_prover::constraint_framework::preprocessed_columns::PreProcessedColumnId;

// [value, carry_4, carry_6, carry_7]
const N_COLUMNS: usize = 4;

relation!(Relation, N_COLUMNS);

pub struct Columns;

impl PreProcessedColumn for Columns {
    /// For columns, one for RangeCheck16, and one for each carry. Unused carries are 0.
    fn log_size(&self) -> Vec<u32> {
        vec![19, 19, 19, 19]
    }

    fn id(&self) -> Vec<PreProcessedColumnId> {
        ["LIMB", "CARRY_4", "CARRY_6", "CARRY_7"]
            .map(|i| PreProcessedColumnId {
                id: format!("range_check_add_{}", i),
            })
            .to_vec()
    }

    fn gen_column_simd(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        let mut all_columns = Vec::with_capacity(N_COLUMNS);

        // I0 lookup
        let domain = CanonicCoset::new(19).circle_domain();

        let columns = (0..1 << 16).flat_map(move |limb| {
            let carry_4 = [0, 1, 2, 3, 0, 0, 0, 0].into_iter();
            let carry_6 = [0, 1, 2, 3, 4, 5, 0, 0].into_iter();
            let carry_7 = [0, 1, 2, 3, 4, 5, 6, 0].into_iter();
            carry_4
                .zip(carry_6)
                .zip(carry_7)
                .map(move |((c_4, c_6), c_7)| {
                    (
                        BaseField::from_u32_unchecked(limb),
                        BaseField::from_u32_unchecked(c_4),
                        BaseField::from_u32_unchecked(c_6),
                        BaseField::from_u32_unchecked(c_7),
                    )
                })
        });

        let (limb, rest) = columns.tee();
        let (carry_4, rest) = rest.tee();
        let (carry_6, carry_7) = rest.tee();

        all_columns.extend(vec![
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(limb.map(|t| t.0)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(carry_4.map(|t| t.1)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(carry_6.map(|t| t.2)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(carry_7.map(|t| t.3)),
            ),
        ]);

        all_columns
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use stwo_prover::core::backend::Column;

    #[test]
    fn test_ids() {
        assert_eq!(Columns.id().len(), N_COLUMNS);

        assert_eq!(
            Columns.id(),
            vec![
                PreProcessedColumnId {
                    id: "range_check_add_LIMB".to_string(),
                },
                PreProcessedColumnId {
                    id: "range_check_add_CARRY_4".to_string(),
                },
                PreProcessedColumnId {
                    id: "range_check_add_CARRY_6".to_string(),
                },
                PreProcessedColumnId {
                    id: "range_check_add_CARRY_7".to_string(),
                },
            ]
        );
    }

    #[test]
    fn test_gen_column_simd() {
        let columns = Columns.gen_column_simd();
        assert_eq!(columns.len(), N_COLUMNS);
        assert_eq!(columns[0].values.len().ilog2(), 19);
        assert_eq!(columns[1].values.len().ilog2(), 19);
        assert_eq!(columns[2].values.len().ilog2(), 19);
        assert_eq!(columns[3].values.len().ilog2(), 19);
    }
}
