use itertools::Itertools;
use stwo::{
    core::{channel::Channel, fields::m31::BaseField, poly::circle::CanonicCoset},
    prover::{
        backend::simd::{column::BaseColumn, SimdBackend},
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::relation;

use crate::{
    partitions::{pext_u32, BigSigma0 as BigSigma0Partitions, SubsetIterator},
    sha256::big_sigma_0,
    trace_columns,
};

const N_IO_COLUMNS: usize = 6;
const N_I1_COLUMNS: usize = 6;
const N_O2_COLUMNS: usize = 4;

relation!(I0, N_IO_COLUMNS);
relation!(I1, N_I1_COLUMNS);
relation!(O2, N_O2_COLUMNS);

trace_columns!(
    Columns,
    big_sigma_0_i0_low,
    big_sigma_0_i0_high_0,
    big_sigma_0_i0_high_1,
    big_sigma_0_o0_low,
    big_sigma_0_o0_high,
    big_sigma_0_o20,
    big_sigma_0_i1_low_0,
    big_sigma_0_i1_low_1,
    big_sigma_0_i1_high,
    big_sigma_0_o1_low,
    big_sigma_0_o1_high,
    big_sigma_0_o21,
    big_sigma_0_o2_0,
    big_sigma_0_o2_1,
    big_sigma_0_o2_low,
    big_sigma_0_o2_high,
);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0: I0,
    pub i1: I1,
    pub o2: O2,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0: I0::dummy(),
            i1: I1::dummy(),
            o2: O2::dummy(),
        }
    }
}

impl Relation {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0: I0::draw(channel),
            i1: I1::draw(channel),
            o2: O2::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    // I0 lookup
    let domain_i0 = CanonicCoset::new(BigSigma0Partitions::I0.count_ones()).circle_domain();
    let i0_columns = SubsetIterator::new(BigSigma0Partitions::I0)
        .map(|x| (x, big_sigma_0(x)))
        .map(|(x, y)| {
            (
                BaseField::from_u32_unchecked(x & BigSigma0Partitions::I0_L),
                BaseField::from_u32_unchecked((x >> 16) & BigSigma0Partitions::I0_H0),
                BaseField::from_u32_unchecked((x >> 24) & BigSigma0Partitions::I0_H1),
                BaseField::from_u32_unchecked(y & BigSigma0Partitions::O0_L),
                BaseField::from_u32_unchecked((y >> 16) & BigSigma0Partitions::O0_H),
                BaseField::from_u32_unchecked(pext_u32(y, BigSigma0Partitions::O2)),
            )
        });

    let (i0_l, rest) = i0_columns.tee();
    let (i0_h_0, rest) = rest.tee();
    let (i0_h_1, rest) = rest.tee();
    let (o0_l, rest) = rest.tee();
    let (o0_h, o20) = rest.tee();

    let i0_columns = vec![
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(i0_l.map(|t| t.0)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(i0_h_0.map(|t| t.1)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(i0_h_1.map(|t| t.2)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(o0_l.map(|t| t.3)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(o0_h.map(|t| t.4)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i0,
            BaseColumn::from_iter(o20.map(|t| t.5)),
        ),
    ];

    // I1 lookup
    let domain_i1 = CanonicCoset::new(BigSigma0Partitions::I1.count_ones()).circle_domain();
    let i1_columns = SubsetIterator::new(BigSigma0Partitions::I1)
        .map(|x| (x, big_sigma_0(x)))
        .map(|(x, y)| {
            (
                BaseField::from_u32_unchecked(x & BigSigma0Partitions::I1_L0),
                BaseField::from_u32_unchecked((x >> 8) & BigSigma0Partitions::I1_L1),
                BaseField::from_u32_unchecked((x >> 16) & BigSigma0Partitions::I1_H),
                BaseField::from_u32_unchecked(y & BigSigma0Partitions::O1_L),
                BaseField::from_u32_unchecked((y >> 16) & BigSigma0Partitions::O1_H),
                BaseField::from_u32_unchecked(pext_u32(y, BigSigma0Partitions::O2)),
            )
        });

    let (i1_l_0, rest) = i1_columns.tee();
    let (i1_l_1, rest) = rest.tee();
    let (i1_h, rest) = rest.tee();
    let (o1_l, rest) = rest.tee();
    let (o1_h, o21) = rest.tee();

    let i1_columns = vec![
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(i1_l_0.map(|t| t.0)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(i1_l_1.map(|t| t.1)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(i1_h.map(|t| t.2)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(o1_l.map(|t| t.3)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(o1_h.map(|t| t.4)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain_i1,
            BaseColumn::from_iter(o21.map(|t| t.5)),
        ),
    ];

    // O2 lookup
    let domain_o2 = CanonicCoset::new(BigSigma0Partitions::O2.count_ones() * 2).circle_domain();
    let o2_columns = SubsetIterator::new(BigSigma0Partitions::O2)
        .flat_map(move |x| SubsetIterator::new(BigSigma0Partitions::O2).map(move |y| (x, y)))
        .map(|(x, y)| {
            (
                BaseField::from_u32_unchecked(pext_u32(x, BigSigma0Partitions::O2)),
                BaseField::from_u32_unchecked(pext_u32(y, BigSigma0Partitions::O2)),
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

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use stwo::prover::backend::Column;

    use super::*;

    const N_COLUMNS: usize = N_IO_COLUMNS + N_I1_COLUMNS + N_O2_COLUMNS;
    #[test]
    fn test_ids() {
        assert_eq!(Columns::to_ids().len(), N_COLUMNS);
    }

    #[test]
    fn test_gen_column_simd() {
        let columns = gen_column_simd();
        assert_eq!(columns.len(), N_COLUMNS);
        // IO lookup
        assert_eq!(
            columns[0].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[1].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[2].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[3].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[4].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        assert_eq!(
            columns[5].values.len().ilog2(),
            BigSigma0Partitions::I0.count_ones()
        );
        // I1 lookup
        assert_eq!(
            columns[6].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[7].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[8].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[9].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[10].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        assert_eq!(
            columns[11].values.len().ilog2(),
            BigSigma0Partitions::I1.count_ones()
        );
        // O2 lookup
        assert_eq!(
            columns[12].values.len().ilog2(),
            BigSigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[13].values.len().ilog2(),
            BigSigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[14].values.len().ilog2(),
            BigSigma0Partitions::O2.count_ones() * 2
        );
        assert_eq!(
            columns[15].values.len().ilog2(),
            BigSigma0Partitions::O2.count_ones() * 2
        );
    }

    #[test]
    fn test_random_input() {
        let columns = gen_column_simd();

        let mut lookup_i0: HashMap<(u32, u32, u32), (u32, u32, u32)> = HashMap::new();
        columns[0]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[1].values.to_cpu())
            .zip(columns[2].values.to_cpu())
            .zip(columns[3].values.to_cpu())
            .zip(columns[4].values.to_cpu())
            .zip(columns[5].values.to_cpu())
            .map(|(((((a, b), c), d), e), f)| ((a.0, b.0, c.0), (d.0, e.0, f.0)))
            .for_each(|(key, value)| {
                lookup_i0.insert(key, value);
            });

        let mut lookup_i1: HashMap<(u32, u32, u32), (u32, u32, u32)> = HashMap::new();
        columns[6]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[7].values.to_cpu())
            .zip(columns[8].values.to_cpu())
            .zip(columns[9].values.to_cpu())
            .zip(columns[10].values.to_cpu())
            .zip(columns[11].values.to_cpu())
            .map(|(((((a, b), c), d), e), f)| ((a.0, b.0, c.0), (d.0, e.0, f.0)))
            .for_each(|(key, value)| {
                lookup_i1.insert(key, value);
            });

        let mut lookup_o2: HashMap<(u32, u32), (u32, u32)> = HashMap::new();
        columns[12]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[13].values.to_cpu())
            .zip(columns[14].values.to_cpu())
            .zip(columns[15].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0), (c.0, d.0)))
            .for_each(|(key, value)| {
                lookup_o2.insert(key, value);
            });

        let (x_low, x_high) = (123456789_u32 & 0xffff, 123456789_u32 >> 16);

        // I0 lookup
        let x_0_low = x_low & BigSigma0Partitions::I0_L;
        let x_0_high_0 = x_high & BigSigma0Partitions::I0_H0;
        let x_0_high_1 = (x_high >> 8) & BigSigma0Partitions::I0_H1;
        let lookup_i0_value = lookup_i0.get(&(x_0_low, x_0_high_0, x_0_high_1));
        assert!(lookup_i0_value.is_some());
        let (o0_l, o0_h, o20) = lookup_i0_value.unwrap();

        // I1 lookup
        let x_1_low_0 = x_low & BigSigma0Partitions::I1_L0;
        let x_1_low_1 = (x_low >> 8) & BigSigma0Partitions::I1_L1;
        let x_1_high = x_high & BigSigma0Partitions::I1_H;
        let lookup_i1_value = lookup_i1.get(&(x_1_low_0, x_1_low_1, x_1_high));
        assert!(lookup_i1_value.is_some());
        let (o1_l, o1_h, o21) = lookup_i1_value.unwrap();

        // 02 lookup
        let lookup_o2_value = lookup_o2.get(&(*o20, *o21));
        assert!(lookup_o2_value.is_some());
        let (o2_l, o2_h) = lookup_o2_value.unwrap();

        // Check the result
        let expected = big_sigma_0(x_low + (x_high << 16));
        assert_eq!(o0_l + o1_l + o2_l, expected & 0xffff);
        assert_eq!(o0_h + o1_h + o2_h, expected >> 16);
    }
}
