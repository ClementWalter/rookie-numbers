use crate::partitions::{BigSigma1 as BigSigma1Partitions, SubsetIterator};
use crate::preprocessed::PreProcessedColumn;
use crate::sha256::ch_right;

use itertools::Itertools;
use stwo_prover::core::backend::simd::column::BaseColumn;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::channel::Channel;
use stwo_prover::core::fields::m31::BaseField;
use stwo_prover::core::poly::circle::{CanonicCoset, CircleEvaluation};
use stwo_prover::core::poly::BitReversedOrder;
use stwo_prover::relation;

use stwo_prover::constraint_framework::preprocessed_columns::PreProcessedColumnId;

// [e, f, val]
const N_COLUMNS: usize = 3;

relation!(I0_L, N_COLUMNS);
relation!(I0_H, N_COLUMNS);
relation!(I1_L, N_COLUMNS);
relation!(I1_H, N_COLUMNS);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0_low: I0_L,
    pub i0_high: I0_H,
    pub i1_low: I1_L,
    pub i1_high: I1_H,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0_low: I0_L::dummy(),
            i0_high: I0_H::dummy(),
            i1_low: I1_L::dummy(),
            i1_high: I1_H::dummy(),
        }
    }

    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0_low: I0_L::draw(channel),
            i0_high: I0_H::draw(channel),
            i1_low: I1_L::draw(channel),
            i1_high: I1_H::draw(channel),
        }
    }
}

pub struct Columns;

impl PreProcessedColumn for Columns {
    fn log_size(&self) -> Vec<u32> {
        vec![
            // I0_L lookup
            BigSigma1Partitions::I0_L.count_ones() * 2,
            BigSigma1Partitions::I0_L.count_ones() * 2,
            BigSigma1Partitions::I0_L.count_ones() * 2,
            // I0_H lookup
            BigSigma1Partitions::I0_H.count_ones() * 2,
            BigSigma1Partitions::I0_H.count_ones() * 2,
            BigSigma1Partitions::I0_H.count_ones() * 2,
            // I1_L lookup
            BigSigma1Partitions::I1_L.count_ones() * 2,
            BigSigma1Partitions::I1_L.count_ones() * 2,
            BigSigma1Partitions::I1_L.count_ones() * 2,
            // I1_H lookup
            BigSigma1Partitions::I1_H.count_ones() * 2,
            BigSigma1Partitions::I1_H.count_ones() * 2,
            BigSigma1Partitions::I1_H.count_ones() * 2,
        ]
    }

    fn id(&self) -> Vec<PreProcessedColumnId> {
        [
            "I0_L_E", "I0_L_F", "I0_L_RES", "I0_H_E", "I0_H_F", "I0_H_RES", "I1_L_E", "I1_L_F",
            "I1_L_RES", "I1_H_E", "I1_H_F", "I1_H_RES",
        ]
        .map(|i| PreProcessedColumnId {
            id: format!("ch_right_{}", i),
        })
        .to_vec()
    }

    fn gen_column_simd(&self) -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
        let mut all_columns = Vec::with_capacity(4 * N_COLUMNS);

        for partition in [
            BigSigma1Partitions::I0_L,
            BigSigma1Partitions::I0_H,
            BigSigma1Partitions::I1_L,
            BigSigma1Partitions::I1_H,
        ] {
            // I0 lookup
            let domain = CanonicCoset::new(partition.count_ones() * 2).circle_domain();
            let columns = SubsetIterator::new(partition).flat_map(move |e| {
                SubsetIterator::new(partition).map(move |f| {
                    (
                        BaseField::from_u32_unchecked(e),
                        BaseField::from_u32_unchecked(f),
                        BaseField::from_u32_unchecked(ch_right(e, f)),
                    )
                })
            });

            let (e, rest) = columns.tee();
            let (f, res) = rest.tee();

            all_columns.extend(vec![
                CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                    domain,
                    BaseColumn::from_iter(e.map(|t| t.0)),
                ),
                CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                    domain,
                    BaseColumn::from_iter(f.map(|t| t.1)),
                ),
                CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                    domain,
                    BaseColumn::from_iter(res.map(|t| t.2)),
                ),
            ]);
        }

        all_columns
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use stwo_prover::core::backend::Column;

    #[test]
    fn test_ids() {
        assert_eq!(Columns.id().len(), N_COLUMNS * 4);

        assert_eq!(
            Columns.id(),
            vec![
                PreProcessedColumnId {
                    id: "ch_right_I0_L_E".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I0_L_F".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I0_L_RES".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I0_H_E".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I0_H_F".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I0_H_RES".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_L_E".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_L_F".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_L_RES".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_H_E".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_H_F".to_string(),
                },
                PreProcessedColumnId {
                    id: "ch_right_I1_H_RES".to_string(),
                },
            ]
        );
    }

    #[test]
    fn test_gen_column_simd() {
        let columns = Columns.gen_column_simd();
        assert_eq!(columns.len(), N_COLUMNS * 4);
        assert_eq!(
            columns[0].values.len().ilog2(),
            BigSigma1Partitions::I0_L.count_ones() * 2
        );
        assert_eq!(
            columns[1].values.len().ilog2(),
            BigSigma1Partitions::I0_L.count_ones() * 2
        );
        assert_eq!(
            columns[2].values.len().ilog2(),
            BigSigma1Partitions::I0_L.count_ones() * 2
        );
        assert_eq!(
            columns[3].values.len().ilog2(),
            BigSigma1Partitions::I0_H.count_ones() * 2
        );
        assert_eq!(
            columns[4].values.len().ilog2(),
            BigSigma1Partitions::I0_H.count_ones() * 2
        );
        assert_eq!(
            columns[5].values.len().ilog2(),
            BigSigma1Partitions::I0_H.count_ones() * 2
        );
        assert_eq!(
            columns[6].values.len().ilog2(),
            BigSigma1Partitions::I1_L.count_ones() * 2
        );
        assert_eq!(
            columns[7].values.len().ilog2(),
            BigSigma1Partitions::I1_L.count_ones() * 2
        );
        assert_eq!(
            columns[8].values.len().ilog2(),
            BigSigma1Partitions::I1_L.count_ones() * 2
        );
        assert_eq!(
            columns[9].values.len().ilog2(),
            BigSigma1Partitions::I1_H.count_ones() * 2
        );
        assert_eq!(
            columns[10].values.len().ilog2(),
            BigSigma1Partitions::I1_H.count_ones() * 2
        );
        assert_eq!(
            columns[11].values.len().ilog2(),
            BigSigma1Partitions::I1_H.count_ones() * 2
        );
    }

    #[test]
    fn test_random_input() {
        let columns = Columns.gen_column_simd();

        let mut lookup_i0_l: HashMap<(u32, u32), u32> = HashMap::new();
        columns[0]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[1].values.to_cpu())
            .zip(columns[2].values.to_cpu())
            .map(|((a, b), c)| ((a.0, b.0), c.0))
            .for_each(|(key, value)| {
                lookup_i0_l.insert(key, value);
            });

        let mut lookup_i0_h: HashMap<(u32, u32), u32> = HashMap::new();
        columns[3]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[4].values.to_cpu())
            .zip(columns[5].values.to_cpu())
            .map(|((a, b), c)| ((a.0, b.0), c.0))
            .for_each(|(key, value)| {
                lookup_i0_h.insert(key, value);
            });

        let mut lookup_i1_l: HashMap<(u32, u32), u32> = HashMap::new();
        columns[6]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[7].values.to_cpu())
            .zip(columns[8].values.to_cpu())
            .map(|((a, b), c)| ((a.0, b.0), c.0))
            .for_each(|(key, value)| {
                lookup_i1_l.insert(key, value);
            });

        let mut lookup_i1_h: HashMap<(u32, u32), u32> = HashMap::new();
        columns[9]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[10].values.to_cpu())
            .zip(columns[11].values.to_cpu())
            .map(|((a, b), c)| ((a.0, b.0), c.0))
            .for_each(|(key, value)| {
                lookup_i1_h.insert(key, value);
            });

        let (e_low, e_high) = (123456789_u32 & 0xffff, 123456789_u32 >> 16);
        let (f_low, f_high) = (987654321_u32 & 0xffff, 987654321_u32 >> 16);

        // I0_L lookup
        let lookup_i0_l_value = lookup_i0_l.get(&(
            e_low & BigSigma1Partitions::I0_L,
            f_low & BigSigma1Partitions::I0_L,
        ));
        assert!(lookup_i0_l_value.is_some());
        let res_i0_l = lookup_i0_l_value.unwrap();

        // I0_H lookup
        let lookup_i0_h_value = lookup_i0_h.get(&(
            e_high & BigSigma1Partitions::I0_H,
            f_high & BigSigma1Partitions::I0_H,
        ));
        assert!(lookup_i0_h_value.is_some());
        let res_i0_h = lookup_i0_h_value.unwrap();

        // I1_L lookup
        let lookup_i1_l_value = lookup_i1_l.get(&(
            e_low & BigSigma1Partitions::I1_L,
            f_low & BigSigma1Partitions::I1_L,
        ));
        assert!(lookup_i1_l_value.is_some());
        let res_i1_l = lookup_i1_l_value.unwrap();

        // I1_H lookup
        let lookup_i1_h_value = lookup_i1_h.get(&(
            e_high & BigSigma1Partitions::I1_H,
            f_high & BigSigma1Partitions::I1_H,
        ));
        assert!(lookup_i1_h_value.is_some());
        let res_i1_h = lookup_i1_h_value.unwrap();

        // Check the result
        let expected = ch_right(e_low + (e_high << 16), f_low + (f_high << 16));
        assert_eq!(res_i0_l + res_i1_l, expected & 0xffff);
        assert_eq!(res_i0_h + res_i1_h, expected >> 16);
    }
}
