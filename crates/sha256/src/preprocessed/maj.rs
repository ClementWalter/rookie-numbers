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
    partitions::{BigSigma0 as BigSigma0Partitions, SubsetIterator},
    sha256::maj,
    trace_columns,
};

// [a, b, c, val]
const N_COLUMNS: usize = 4;

relation!(I0_L, N_COLUMNS);
relation!(I0_H0, N_COLUMNS);
relation!(I0_H1, N_COLUMNS);
relation!(I1_L0, N_COLUMNS);
relation!(I1_L1, N_COLUMNS);
relation!(I1_H, N_COLUMNS);

trace_columns!(
    Columns,
    maj_i0_low_a,
    maj_i0_low_b,
    maj_i0_low_c,
    maj_i0_low_res,
    maj_i0_high_0_a,
    maj_i0_high_0_b,
    maj_i0_high_0_c,
    maj_i0_high_0_res,
    maj_i0_high_1_a,
    maj_i0_high_1_b,
    maj_i0_high_1_c,
    maj_i0_high_1_res,
    maj_i1_low_0_a,
    maj_i1_low_0_b,
    maj_i1_low_0_c,
    maj_i1_low_0_res,
    maj_i1_low_1_a,
    maj_i1_low_1_b,
    maj_i1_low_1_c,
    maj_i1_low_1_res,
    maj_i1_high_a,
    maj_i1_high_b,
    maj_i1_high_c,
    maj_i1_high_res,
);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0_low: I0_L,
    pub i0_high_0: I0_H0,
    pub i0_high_1: I0_H1,
    pub i1_low_0: I1_L0,
    pub i1_low_1: I1_L1,
    pub i1_high: I1_H,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0_low: I0_L::dummy(),
            i0_high_0: I0_H0::dummy(),
            i0_high_1: I0_H1::dummy(),
            i1_low_0: I1_L0::dummy(),
            i1_low_1: I1_L1::dummy(),
            i1_high: I1_H::dummy(),
        }
    }

    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0_low: I0_L::draw(channel),
            i0_high_0: I0_H0::draw(channel),
            i0_high_1: I0_H1::draw(channel),
            i1_low_0: I1_L0::draw(channel),
            i1_low_1: I1_L1::draw(channel),
            i1_high: I1_H::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let mut all_columns = Vec::with_capacity(6 * N_COLUMNS);

    for partition in [
        BigSigma0Partitions::I0_L,
        BigSigma0Partitions::I0_H0,
        BigSigma0Partitions::I0_H1,
        BigSigma0Partitions::I1_L0,
        BigSigma0Partitions::I1_L1,
        BigSigma0Partitions::I1_H,
    ] {
        let domain = CanonicCoset::new(partition.count_ones() * 3).circle_domain();
        let columns = SubsetIterator::new(partition)
            .flat_map(move |x| SubsetIterator::new(partition).map(move |y| (x, y)))
            .flat_map(move |x| {
                SubsetIterator::new(partition).map(move |y| {
                    (
                        BaseField::from_u32_unchecked(x.0),
                        BaseField::from_u32_unchecked(x.1),
                        BaseField::from_u32_unchecked(y),
                        BaseField::from_u32_unchecked(maj(x.0, x.1, y)),
                    )
                })
            });

        let (a, rest) = columns.tee();
        let (b, rest) = rest.tee();
        let (c, res) = rest.tee();

        all_columns.extend(vec![
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(a.map(|t| t.0)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(b.map(|t| t.1)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(c.map(|t| t.2)),
            ),
            CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                domain,
                BaseColumn::from_iter(res.map(|t| t.3)),
            ),
        ]);
    }

    all_columns
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use stwo::prover::backend::Column;

    use super::*;

    #[test]
    fn test_ids() {
        assert_eq!(Columns::to_ids().len(), N_COLUMNS * 6);
    }

    #[allow(clippy::cognitive_complexity)]
    #[test]
    fn test_gen_column_simd() {
        let columns = gen_column_simd();
        assert_eq!(columns.len(), N_COLUMNS * 6);
        assert_eq!(
            columns[0].values.len().ilog2(),
            BigSigma0Partitions::I0_L.count_ones() * 3
        );
        assert_eq!(
            columns[1].values.len().ilog2(),
            BigSigma0Partitions::I0_L.count_ones() * 3
        );
        assert_eq!(
            columns[2].values.len().ilog2(),
            BigSigma0Partitions::I0_L.count_ones() * 3
        );
        assert_eq!(
            columns[3].values.len().ilog2(),
            BigSigma0Partitions::I0_L.count_ones() * 3
        );
        assert_eq!(
            columns[4].values.len().ilog2(),
            BigSigma0Partitions::I0_H0.count_ones() * 3
        );
        assert_eq!(
            columns[5].values.len().ilog2(),
            BigSigma0Partitions::I0_H0.count_ones() * 3
        );
        assert_eq!(
            columns[6].values.len().ilog2(),
            BigSigma0Partitions::I0_H0.count_ones() * 3
        );
        assert_eq!(
            columns[7].values.len().ilog2(),
            BigSigma0Partitions::I0_H0.count_ones() * 3
        );
        assert_eq!(
            columns[8].values.len().ilog2(),
            BigSigma0Partitions::I0_H1.count_ones() * 3
        );
        assert_eq!(
            columns[9].values.len().ilog2(),
            BigSigma0Partitions::I0_H1.count_ones() * 3
        );
        assert_eq!(
            columns[10].values.len().ilog2(),
            BigSigma0Partitions::I0_H1.count_ones() * 3
        );
        assert_eq!(
            columns[11].values.len().ilog2(),
            BigSigma0Partitions::I0_H1.count_ones() * 3
        );
        assert_eq!(
            columns[12].values.len().ilog2(),
            BigSigma0Partitions::I1_L0.count_ones() * 3
        );
        assert_eq!(
            columns[13].values.len().ilog2(),
            BigSigma0Partitions::I1_L0.count_ones() * 3
        );
        assert_eq!(
            columns[14].values.len().ilog2(),
            BigSigma0Partitions::I1_L0.count_ones() * 3
        );
        assert_eq!(
            columns[15].values.len().ilog2(),
            BigSigma0Partitions::I1_L0.count_ones() * 3
        );
        assert_eq!(
            columns[16].values.len().ilog2(),
            BigSigma0Partitions::I1_L1.count_ones() * 3
        );
        assert_eq!(
            columns[17].values.len().ilog2(),
            BigSigma0Partitions::I1_L1.count_ones() * 3
        );
        assert_eq!(
            columns[18].values.len().ilog2(),
            BigSigma0Partitions::I1_L1.count_ones() * 3
        );
        assert_eq!(
            columns[19].values.len().ilog2(),
            BigSigma0Partitions::I1_L1.count_ones() * 3
        );
        assert_eq!(
            columns[20].values.len().ilog2(),
            BigSigma0Partitions::I1_H.count_ones() * 3
        );
        assert_eq!(
            columns[21].values.len().ilog2(),
            BigSigma0Partitions::I1_H.count_ones() * 3
        );
        assert_eq!(
            columns[22].values.len().ilog2(),
            BigSigma0Partitions::I1_H.count_ones() * 3
        );
        assert_eq!(
            columns[23].values.len().ilog2(),
            BigSigma0Partitions::I1_H.count_ones() * 3
        );
    }

    #[test]
    fn test_random_input() {
        let columns = gen_column_simd();

        let mut lookup_i0_l: HashMap<(u32, u32, u32), u32> = HashMap::new();
        columns[0]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[1].values.to_cpu())
            .zip(columns[2].values.to_cpu())
            .zip(columns[3].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i0_l.insert(key, value);
            });

        let mut lookup_i0_h_0: HashMap<(u32, u32, u32), u32> = HashMap::new();
        columns[4]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[5].values.to_cpu())
            .zip(columns[6].values.to_cpu())
            .zip(columns[7].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i0_h_0.insert(key, value);
            });

        let mut lookup_i0_h_1: HashMap<(u32, u32, u32), u32> = HashMap::new();
        columns[8]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[9].values.to_cpu())
            .zip(columns[10].values.to_cpu())
            .zip(columns[11].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i0_h_1.insert(key, value);
            });

        let mut lookup_i1_l_0: HashMap<(u32, u32, u32), u32> = HashMap::new();
        columns[12]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[13].values.to_cpu())
            .zip(columns[14].values.to_cpu())
            .zip(columns[15].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i1_l_0.insert(key, value);
            });

        let mut lookup_i1_l_1: HashMap<(u32, u32, u32), u32> = HashMap::new();

        columns[16]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[17].values.to_cpu())
            .zip(columns[18].values.to_cpu())
            .zip(columns[19].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i1_l_1.insert(key, value);
            });

        let mut lookup_i1_h: HashMap<(u32, u32, u32), u32> = HashMap::new();
        columns[20]
            .values
            .to_cpu()
            .into_iter()
            .zip(columns[21].values.to_cpu())
            .zip(columns[22].values.to_cpu())
            .zip(columns[23].values.to_cpu())
            .map(|(((a, b), c), d)| ((a.0, b.0, c.0), d.0))
            .for_each(|(key, value)| {
                lookup_i1_h.insert(key, value);
            });

        let (a_low, a_high) = (123456789_u32 & 0xffff, 123456789_u32 >> 16);
        let (b_low, b_high) = (987654321_u32 & 0xffff, 987654321_u32 >> 16);
        let (c_low, c_high) = (543219876_u32 & 0xffff, 543219876_u32 >> 16);

        // I0_L lookup
        let lookup_i0_l_value = lookup_i0_l.get(&(
            a_low & BigSigma0Partitions::I0_L,
            b_low & BigSigma0Partitions::I0_L,
            c_low & BigSigma0Partitions::I0_L,
        ));
        assert!(lookup_i0_l_value.is_some());
        let res_i0_l = lookup_i0_l_value.unwrap();

        // I0_H0 lookup
        let lookup_i0_h_0_value = lookup_i0_h_0.get(&(
            a_high & BigSigma0Partitions::I0_H0,
            b_high & BigSigma0Partitions::I0_H0,
            c_high & BigSigma0Partitions::I0_H0,
        ));
        assert!(lookup_i0_h_0_value.is_some());
        let res_i0_h_0 = lookup_i0_h_0_value.unwrap();

        // I0_H1 lookup
        let lookup_i0_h_1_value = lookup_i0_h_1.get(&(
            (a_high >> 8) & BigSigma0Partitions::I0_H1,
            (b_high >> 8) & BigSigma0Partitions::I0_H1,
            (c_high >> 8) & BigSigma0Partitions::I0_H1,
        ));
        assert!(lookup_i0_h_1_value.is_some());
        let res_i0_h_1 = lookup_i0_h_1_value.unwrap();

        // I1_L0 lookup
        let lookup_i1_l_0_value = lookup_i1_l_0.get(&(
            a_low & BigSigma0Partitions::I1_L0,
            b_low & BigSigma0Partitions::I1_L0,
            c_low & BigSigma0Partitions::I1_L0,
        ));
        assert!(lookup_i1_l_0_value.is_some());
        let res_i1_l_0 = lookup_i1_l_0_value.unwrap();

        // I1_L1 lookup
        let lookup_i1_l_1_value = lookup_i1_l_1.get(&(
            (a_low >> 8) & BigSigma0Partitions::I1_L1,
            (b_low >> 8) & BigSigma0Partitions::I1_L1,
            (c_low >> 8) & BigSigma0Partitions::I1_L1,
        ));
        assert!(lookup_i1_l_1_value.is_some());
        let res_i1_l_1 = lookup_i1_l_1_value.unwrap();

        // I1_H lookup
        let lookup_i1_h_value = lookup_i1_h.get(&(
            a_high & BigSigma0Partitions::I1_H,
            b_high & BigSigma0Partitions::I1_H,
            c_high & BigSigma0Partitions::I1_H,
        ));
        assert!(lookup_i1_h_value.is_some());
        let res_i1_h = lookup_i1_h_value.unwrap();

        // Check the result
        let expected = maj(
            a_low + (a_high << 16),
            b_low + (b_high << 16),
            c_low + (c_high << 16),
        );
        assert_eq!(res_i0_l + res_i1_l_0 + (res_i1_l_1 << 8), expected & 0xffff);
        assert_eq!(res_i0_h_0 + (res_i0_h_1 << 8) + res_i1_h, expected >> 16);
    }
}
