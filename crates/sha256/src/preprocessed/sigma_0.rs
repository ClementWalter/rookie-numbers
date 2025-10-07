use crate::partitions::{Sigma0 as Sigma0Partitions, SubsetIterator};
use crate::preprocessed::PreProcessedColumn;
use crate::sha256::small_sigma0;

use itertools::Itertools;
use stwo_prover::core::backend::simd::column::BaseColumn;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::fields::m31::BaseField;
use stwo_prover::core::poly::circle::{CanonicCoset, CircleEvaluation};
use stwo_prover::core::poly::BitReversedOrder;
use stwo_prover::relation;

use stwo_prover::constraint_framework::preprocessed_columns::PreProcessedColumnId;

const N_IO_COLUMNS: usize = 5;
const N_I1_COLUMNS: usize = 5;
const N_O2_COLUMNS: usize = 4;

relation!(Sigma0, 6);
/// Lookup data for the Sigma0 function.
/// The sigma0 function is emulated with 3 lookups, one for each partition I0, I1,
/// and a final lookup for O2 xor.
pub struct Sigma0LookupData {
    pub i0: [BaseColumn; N_IO_COLUMNS], // [i0_l, i0_h, o0_l, o0_h, o20]
    pub i1: [BaseColumn; N_I1_COLUMNS], // [i1_l, i1_h, o1_l, o1_h, o21]
    pub o2: [BaseColumn; N_O2_COLUMNS], // [o20, o21, o2_l, o2_h]
}

pub struct Sigma0Columns {}

impl PreProcessedColumn for Sigma0Columns {
    fn log_size(&self) -> Vec<u32> {
        vec![
            // IO lookup
            Sigma0Partitions::I0.count_ones(),
            Sigma0Partitions::I0.count_ones(),
            Sigma0Partitions::I0.count_ones(),
            Sigma0Partitions::I0.count_ones(),
            Sigma0Partitions::I0.count_ones(),
            // I1 lookup
            Sigma0Partitions::I1.count_ones(),
            Sigma0Partitions::I1.count_ones(),
            Sigma0Partitions::I1.count_ones(),
            Sigma0Partitions::I1.count_ones(),
            Sigma0Partitions::I1.count_ones(),
            // O2 lookup
            Sigma0Partitions::O2.count_ones(),
            Sigma0Partitions::O2.count_ones(),
            Sigma0Partitions::O2.count_ones(),
            Sigma0Partitions::O2.count_ones(),
        ]
    }

    fn id(&self) -> Vec<PreProcessedColumnId> {
        [
            "I0_L", "I0_H", "O0_L", "O0_H", "O20", // IO lookup
            "I1_L", "I1_H", "O1_L", "O1_H", "O21", // I1 lookup
            "O2_0", "O2_1", "O2_L", "O2_H", // O2 lookup
        ]
        .map(|i| PreProcessedColumnId {
            id: format!("Sigma0_{}", i),
        })
        .to_vec()
    }

    fn gen_column_simd(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        // I0 lookup
        let domain_i0 = CanonicCoset::new(Sigma0Partitions::I0.count_ones()).circle_domain();
        let i0_columns = SubsetIterator::new(Sigma0Partitions::I0)
            .map(|x| (x, small_sigma0(x)))
            .map(|(x, y)| {
                (
                    BaseField::from_u32_unchecked(x & Sigma0Partitions::I0_L),
                    BaseField::from_u32_unchecked(x & Sigma0Partitions::I0_H),
                    BaseField::from_u32_unchecked(y & Sigma0Partitions::O0_L),
                    BaseField::from_u32_unchecked((y >> 16) & Sigma0Partitions::O0_H),
                    BaseField::from_u32_unchecked(y & Sigma0Partitions::O2),
                )
            });

        let (i0_l, rest) = i0_columns.tee();
        let (i0_h, rest) = rest.tee();
        let (o0_l, rest) = rest.tee();
        let (o0_h, o20) = rest.tee();

        let i0_columns = vec![
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i0,
                BaseColumn::from_iter(i0_l.map(|t| t.0)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i0,
                BaseColumn::from_iter(i0_h.map(|t| t.1)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i0,
                BaseColumn::from_iter(o0_l.map(|t| t.2)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i0,
                BaseColumn::from_iter(o0_h.map(|t| t.3)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i0,
                BaseColumn::from_iter(o20.map(|t| t.4)),
            ),
        ];

        // I1 lookup
        let domain_i1 = CanonicCoset::new(Sigma0Partitions::I1.count_ones()).circle_domain();
        let i1_columns = SubsetIterator::new(Sigma0Partitions::I1)
            .map(|x| (x, small_sigma0(x)))
            .map(|(x, y)| {
                (
                    BaseField::from_u32_unchecked(x & Sigma0Partitions::I1_L),
                    BaseField::from_u32_unchecked(x & Sigma0Partitions::I1_H),
                    BaseField::from_u32_unchecked(y & Sigma0Partitions::O1_L),
                    BaseField::from_u32_unchecked((y >> 16) & Sigma0Partitions::O1_H),
                    BaseField::from_u32_unchecked(y & Sigma0Partitions::O2),
                )
            });

        let (i1_l, rest) = i1_columns.tee();
        let (i1_h, rest) = rest.tee();
        let (o1_l, rest) = rest.tee();
        let (o1_h, o21) = rest.tee();

        let i1_columns = vec![
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i1,
                BaseColumn::from_iter(i1_l.map(|t| t.0)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i1,
                BaseColumn::from_iter(i1_h.map(|t| t.1)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i1,
                BaseColumn::from_iter(o1_l.map(|t| t.2)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i1,
                BaseColumn::from_iter(o1_h.map(|t| t.3)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_i1,
                BaseColumn::from_iter(o21.map(|t| t.4)),
            ),
        ];

        // O2 lookup
        let domain_o2 = CanonicCoset::new(Sigma0Partitions::O2.count_ones() * 2).circle_domain();
        let o2_columns = SubsetIterator::new(Sigma0Partitions::O2)
            .flat_map(move |x| SubsetIterator::new(Sigma0Partitions::O2).map(move |y| (x, y)))
            .map(|(x, y)| {
                (
                    BaseField::from_u32_unchecked(x),
                    BaseField::from_u32_unchecked(y),
                    BaseField::from_u32_unchecked((x ^ y) & 0xffff),
                    BaseField::from_u32_unchecked((x ^ y) >> 16),
                )
            });

        let (o2_0, rest) = o2_columns.tee();
        let (o2_1, rest) = rest.tee();
        let (o2_l, o2_h) = rest.tee();

        let o2_columns = vec![
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_o2,
                BaseColumn::from_iter(o2_0.map(|t| t.0)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_o2,
                BaseColumn::from_iter(o2_1.map(|t| t.1)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_o2,
                BaseColumn::from_iter(o2_l.map(|t| t.2)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain_o2,
                BaseColumn::from_iter(o2_h.map(|t| t.3)),
            ),
        ];

        i0_columns
            .into_iter()
            .chain(i1_columns)
            .chain(o2_columns)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use stwo_prover::core::backend::Column;

    const N_COLUMNS: usize = N_IO_COLUMNS + N_I1_COLUMNS + N_O2_COLUMNS;
    #[test]
    fn test_ids() {
        let sigma0 = Sigma0Columns {};
        assert_eq!(sigma0.id().len(), N_COLUMNS);
        assert_eq!(
            sigma0.id(),
            vec![
                PreProcessedColumnId {
                    id: "Sigma0_I0_L".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_I0_H".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O0_L".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O0_H".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O20".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_I1_L".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_I1_H".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O1_L".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O1_H".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O21".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O2_0".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O2_1".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O2_L".to_string(),
                },
                PreProcessedColumnId {
                    id: "Sigma0_O2_H".to_string(),
                },
            ]
        );
    }

    #[test]
    fn test_gen_column_simd() {
        let sigma0 = Sigma0Columns {};
        let columns = sigma0.gen_column_simd();
        assert_eq!(columns.len(), N_COLUMNS);
        assert_eq!(
            columns[0].values.len().ilog2(),
            Sigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[1].values.len().ilog2(),
            Sigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[2].values.len().ilog2(),
            Sigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[3].values.len().ilog2(),
            Sigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[4].values.len().ilog2(),
            Sigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[5].values.len().ilog2(),
            Sigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[6].values.len().ilog2(),
            Sigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[7].values.len().ilog2(),
            Sigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[8].values.len().ilog2(),
            Sigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[9].values.len().ilog2(),
            Sigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[10].values.len().ilog2(),
            Sigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[11].values.len().ilog2(),
            Sigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[12].values.len().ilog2(),
            Sigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[13].values.len().ilog2(),
            Sigma0Partitions::O2.count_ones() * 2
        );
    }
}
