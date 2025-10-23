use std::simd::u32x16;

use itertools::{iproduct, Itertools};
use stwo::core::channel::Channel;
use stwo_constraint_framework::relation;
use utils::trace_columns;

use crate::{
    partitions::{BigSigma1 as BigSigma1Partitions, SubsetIterator},
    sha256::ch_right_u32x16,
};

// [e, f, val]
const N_COLUMNS: usize = 3;

relation!(CH_RIGHT_I0_L, N_COLUMNS);
relation!(CH_RIGHT_I0_H, N_COLUMNS);
relation!(CH_RIGHT_I1_L, N_COLUMNS);
relation!(CH_RIGHT_I1_H, N_COLUMNS);

trace_columns!(
    ChRightColumns,
    i0_low_e,
    i0_low_g,
    i0_low_res,
    i0_high_e,
    i0_high_g,
    i0_high_res,
    i1_low_e,
    i1_low_g,
    i1_low_res,
    i1_high_e,
    i1_high_g,
    i1_high_res,
);
trace_columns!(
    ChRightI0Columns,
    i0_low_e,
    i0_low_g,
    i0_low_res,
    i0_high_e,
    i0_high_g,
    i0_high_res,
);
trace_columns!(
    ChRightI1Columns,
    i1_low_e,
    i1_low_g,
    i1_low_res,
    i1_high_e,
    i1_high_g,
    i1_high_res,
);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0_low: CH_RIGHT_I0_L,
    pub i0_high: CH_RIGHT_I0_H,
    pub i1_low: CH_RIGHT_I1_L,
    pub i1_high: CH_RIGHT_I1_H,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0_low: CH_RIGHT_I0_L::dummy(),
            i0_high: CH_RIGHT_I0_H::dummy(),
            i1_low: CH_RIGHT_I1_L::dummy(),
            i1_high: CH_RIGHT_I1_H::dummy(),
        }
    }
}

impl Relation {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0_low: CH_RIGHT_I0_L::draw(channel),
            i0_high: CH_RIGHT_I0_H::draw(channel),
            i1_low: CH_RIGHT_I1_L::draw(channel),
            i1_high: CH_RIGHT_I1_H::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<Vec<u32x16>> {
    let mut all_columns: Vec<Vec<u32x16>> = vec![Vec::new(); ChRightColumns::SIZE];

    for (i, partition) in [
        BigSigma1Partitions::I0_L,
        BigSigma1Partitions::I0_H,
        BigSigma1Partitions::I1_L,
        BigSigma1Partitions::I1_H,
    ]
    .iter()
    .enumerate()
    {
        // lookup
        let tuples: Vec<(u32x16, u32x16, u32x16)> = iproduct!(
            SubsetIterator::new(*partition),
            SubsetIterator::new(*partition)
        )
        .chunks(16)
        .into_iter()
        .map(|chunk| {
            let (xs, ys): (Vec<_>, Vec<_>) = chunk.collect::<Vec<_>>().into_iter().unzip();
            (u32x16::from_slice(&xs), u32x16::from_slice(&ys))
        })
        .map(|(e, f)| (e, f, ch_right_u32x16(e, f)))
        .collect();

        for (e, f, res) in tuples {
            all_columns[3 * i].push(e);
            all_columns[3 * i + 1].push(f);
            all_columns[3 * i + 2].push(res);
        }
    }

    all_columns
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::sha256::ch_right;

    #[test]
    fn test_gen_column_simd() {
        let columns = gen_column_simd();
        assert_eq!(columns.len(), ChRightColumns::SIZE);
    }

    #[test]
    fn test_random_input() {
        let columns = gen_column_simd();

        // Helper to flatten Vec<u32x16> into Vec<u32>
        fn flatten_simd_column(col: &[u32x16]) -> Vec<u32> {
            col.iter()
                .flat_map(|v| v.as_array().iter().copied())
                .collect()
        }

        let mut lookup_i0_low: HashMap<(u32, u32), u32> = HashMap::new();
        let i0_low_e = flatten_simd_column(&columns[0]);
        let i0_low_g = flatten_simd_column(&columns[1]);
        let i0_low_res = flatten_simd_column(&columns[2]);
        for ((e, f), res) in i0_low_e.iter().zip(&i0_low_g).zip(&i0_low_res) {
            lookup_i0_low.insert((*e, *f), *res);
        }

        let mut lookup_i0_high: HashMap<(u32, u32), u32> = HashMap::new();
        let i0_high_e = flatten_simd_column(&columns[3]);
        let i0_high_g = flatten_simd_column(&columns[4]);
        let i0_high_res = flatten_simd_column(&columns[5]);
        for ((e, f), res) in i0_high_e.iter().zip(&i0_high_g).zip(&i0_high_res) {
            lookup_i0_high.insert((*e, *f), *res);
        }

        let mut lookup_i1_low: HashMap<(u32, u32), u32> = HashMap::new();
        let i1_low_e = flatten_simd_column(&columns[6]);
        let i1_low_g = flatten_simd_column(&columns[7]);
        let i1_low_res = flatten_simd_column(&columns[8]);
        for ((e, f), res) in i1_low_e.iter().zip(&i1_low_g).zip(&i1_low_res) {
            lookup_i1_low.insert((*e, *f), *res);
        }

        let mut lookup_i1_high: HashMap<(u32, u32), u32> = HashMap::new();
        let i1_high_e = flatten_simd_column(&columns[9]);
        let i1_high_g = flatten_simd_column(&columns[10]);
        let i1_high_res = flatten_simd_column(&columns[11]);
        for ((e, f), res) in i1_high_e.iter().zip(&i1_high_g).zip(&i1_high_res) {
            lookup_i1_high.insert((*e, *f), *res);
        }

        let (e_low, e_high) = (123456789_u32 & 0xffff, 123456789_u32 >> 16);
        let (f_low, f_high) = (987654321_u32 & 0xffff, 987654321_u32 >> 16);

        // I0_L lookup
        let lookup_i0_low_value = lookup_i0_low.get(&(
            e_low & BigSigma1Partitions::I0_L,
            f_low & BigSigma1Partitions::I0_L,
        ));
        assert!(lookup_i0_low_value.is_some());
        let res_i0_l = lookup_i0_low_value.unwrap();

        // I0_H lookup
        let lookup_i0_high_value = lookup_i0_high.get(&(
            e_high & BigSigma1Partitions::I0_H,
            f_high & BigSigma1Partitions::I0_H,
        ));
        assert!(lookup_i0_high_value.is_some());
        let res_i0_h = lookup_i0_high_value.unwrap();

        // I1_L lookup
        let lookup_i1_low_value = lookup_i1_low.get(&(
            e_low & BigSigma1Partitions::I1_L,
            f_low & BigSigma1Partitions::I1_L,
        ));
        assert!(lookup_i1_low_value.is_some());
        let res_i1_l = lookup_i1_low_value.unwrap();

        // I1_H lookup
        let lookup_i1_high_value = lookup_i1_high.get(&(
            e_high & BigSigma1Partitions::I1_H,
            f_high & BigSigma1Partitions::I1_H,
        ));
        assert!(lookup_i1_high_value.is_some());
        let res_i1_h = lookup_i1_high_value.unwrap();

        // Check the result
        let expected = ch_right(e_low + (e_high << 16), f_low + (f_high << 16));
        assert_eq!(res_i0_l + res_i1_l, expected & 0xffff);
        assert_eq!(res_i0_h + res_i1_h, expected >> 16);
    }
}
