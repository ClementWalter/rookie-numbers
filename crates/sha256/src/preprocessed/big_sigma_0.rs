use std::simd::u32x16;

use itertools::{iproduct, Itertools};
use stwo::core::channel::Channel;
use stwo_constraint_framework::relation;
use utils::trace_columns;

use crate::{
    partitions::{pext_u32x16, BigSigma0 as BigSigma0Partitions, SubsetIterator},
    sha256::big_sigma_0_u32x16,
};

const N_IO_COLUMNS: usize = 6;
const N_I1_COLUMNS: usize = 6;
const N_O2_COLUMNS: usize = 4;

relation!(BIG_SIGMA_0_I0, N_IO_COLUMNS);
relation!(BIG_SIGMA_0_I1, N_I1_COLUMNS);
relation!(BIG_SIGMA_0_O2, N_O2_COLUMNS);

trace_columns!(
    BigSigma0Columns,
    i0_low,
    i0_high_0,
    i0_high_1,
    o0_low,
    o0_high,
    o20_pext,
    i1_low_0,
    i1_low_1,
    i1_high,
    o1_low,
    o1_high,
    o21_pext,
    o2_0,
    o2_1,
    o2_low,
    o2_high,
);

trace_columns!(
    BigSigma0I0I1Columns,
    i0_low,
    i0_high_0,
    i0_high_1,
    o0_low,
    o0_high,
    o20_pext,
    i1_low_0,
    i1_low_1,
    i1_high,
    o1_low,
    o1_high,
    o21_pext,
);

trace_columns!(BigSigma0O2Columns, o2_0, o2_1, o2_low, o2_high);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0: BIG_SIGMA_0_I0,
    pub i1: BIG_SIGMA_0_I1,
    pub o2: BIG_SIGMA_0_O2,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0: BIG_SIGMA_0_I0::dummy(),
            i1: BIG_SIGMA_0_I1::dummy(),
            o2: BIG_SIGMA_0_O2::dummy(),
        }
    }
}

impl Relation {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0: BIG_SIGMA_0_I0::draw(channel),
            i1: BIG_SIGMA_0_I1::draw(channel),
            o2: BIG_SIGMA_0_O2::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<Vec<u32x16>> {
    // I0 lookup
    let i0_tuples: Vec<(u32x16, u32x16, u32x16, u32x16, u32x16, u32x16)> =
        SubsetIterator::new(BigSigma0Partitions::I0)
            .collect::<Vec<u32>>()
            .chunks(16)
            .map(u32x16::from_slice)
            .map(|x| (x, big_sigma_0_u32x16(x)))
            .map(|(x, y)| {
                (
                    x & u32x16::splat(BigSigma0Partitions::I0_L),
                    (x >> 16) & u32x16::splat(BigSigma0Partitions::I0_H0),
                    (x >> 24) & u32x16::splat(BigSigma0Partitions::I0_H1),
                    y & u32x16::splat(BigSigma0Partitions::O0_L),
                    (y >> 16) & u32x16::splat(BigSigma0Partitions::O0_H),
                    pext_u32x16(y, BigSigma0Partitions::O2),
                )
            })
            .collect();

    let mut i0_columns: Vec<Vec<u32x16>> = vec![Vec::new(); N_IO_COLUMNS];
    for (a, b, c, d, e, f) in i0_tuples.into_iter() {
        i0_columns[0].push(a);
        i0_columns[1].push(b);
        i0_columns[2].push(c);
        i0_columns[3].push(d);
        i0_columns[4].push(e);
        i0_columns[5].push(f);
    }

    // I1 lookup
    let i1_tuples: Vec<(u32x16, u32x16, u32x16, u32x16, u32x16, u32x16)> =
        SubsetIterator::new(BigSigma0Partitions::I1)
            .collect::<Vec<u32>>()
            .chunks(16)
            .map(u32x16::from_slice)
            .map(|x| (x, big_sigma_0_u32x16(x)))
            .map(|(x, y)| {
                (
                    x & u32x16::splat(BigSigma0Partitions::I1_L0),
                    (x >> 8) & u32x16::splat(BigSigma0Partitions::I1_L1),
                    (x >> 16) & u32x16::splat(BigSigma0Partitions::I1_H),
                    y & u32x16::splat(BigSigma0Partitions::O1_L),
                    (y >> 16) & u32x16::splat(BigSigma0Partitions::O1_H),
                    pext_u32x16(y, BigSigma0Partitions::O2),
                )
            })
            .collect();

    let mut i1_columns: Vec<Vec<u32x16>> = vec![Vec::new(); N_I1_COLUMNS];
    for (a, b, c, d, e, f) in i1_tuples.into_iter() {
        i1_columns[0].push(a);
        i1_columns[1].push(b);
        i1_columns[2].push(c);
        i1_columns[3].push(d);
        i1_columns[4].push(e);
        i1_columns[5].push(f);
    }

    // O2 lookup
    let o2_tuples: Vec<(u32x16, u32x16, u32x16, u32x16)> = iproduct!(
        SubsetIterator::new(BigSigma0Partitions::O2),
        SubsetIterator::new(BigSigma0Partitions::O2)
    )
    .chunks(16)
    .into_iter()
    .map(|chunk| {
        let (xs, ys): (Vec<_>, Vec<_>) = chunk.collect::<Vec<_>>().into_iter().unzip();
        (u32x16::from_slice(&xs), u32x16::from_slice(&ys))
    })
    .map(|(x, y)| {
        (
            pext_u32x16(x, BigSigma0Partitions::O2),
            pext_u32x16(y, BigSigma0Partitions::O2),
            (x ^ y) & u32x16::splat(0xffff),
            (x ^ y) >> u32x16::splat(16),
        )
    })
    .collect();

    let mut o2_columns: Vec<Vec<u32x16>> = vec![Vec::new(); N_O2_COLUMNS];
    for (a, b, c, d) in o2_tuples.into_iter() {
        o2_columns[0].push(a);
        o2_columns[1].push(b);
        o2_columns[2].push(c);
        o2_columns[3].push(d);
    }

    i0_columns
        .into_iter()
        .chain(i1_columns)
        .chain(o2_columns)
        .collect()
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::sha256::big_sigma_0;

    #[test]
    fn test_ids() {
        assert_eq!(
            BigSigma0Columns::to_ids(None).len(),
            BigSigma0I0I1Columns::SIZE + BigSigma0O2Columns::SIZE
        );
    }

    #[test]
    fn test_random_input() {
        let columns = gen_column_simd();

        // Helper to flatten Vec<u32x16> into Vec<u32>
        fn flatten_simd(col: &[u32x16]) -> Vec<u32> {
            col.iter()
                .flat_map(|v| v.as_array().iter().copied())
                .collect()
        }

        let mut lookup_i0: HashMap<(u32, u32, u32), (u32, u32, u32)> = HashMap::new();
        let i0_0 = flatten_simd(&columns[0]);
        let i0_1 = flatten_simd(&columns[1]);
        let i0_2 = flatten_simd(&columns[2]);
        let i0_3 = flatten_simd(&columns[3]);
        let i0_4 = flatten_simd(&columns[4]);
        let i0_5 = flatten_simd(&columns[5]);
        for (((((&a, &b), &c), &d), &e), &f) in i0_0
            .iter()
            .zip(&i0_1)
            .zip(&i0_2)
            .zip(&i0_3)
            .zip(&i0_4)
            .zip(&i0_5)
        {
            lookup_i0.insert((a, b, c), (d, e, f));
        }

        let mut lookup_i1: HashMap<(u32, u32, u32), (u32, u32, u32)> = HashMap::new();
        let i1_0 = flatten_simd(&columns[6]);
        let i1_1 = flatten_simd(&columns[7]);
        let i1_2 = flatten_simd(&columns[8]);
        let i1_3 = flatten_simd(&columns[9]);
        let i1_4 = flatten_simd(&columns[10]);
        let i1_5 = flatten_simd(&columns[11]);
        for (((((&a, &b), &c), &d), &e), &f) in i1_0
            .iter()
            .zip(&i1_1)
            .zip(&i1_2)
            .zip(&i1_3)
            .zip(&i1_4)
            .zip(&i1_5)
        {
            lookup_i1.insert((a, b, c), (d, e, f));
        }

        let mut lookup_o2: HashMap<(u32, u32), (u32, u32)> = HashMap::new();
        let o2_0 = flatten_simd(&columns[12]);
        let o2_1 = flatten_simd(&columns[13]);
        let o2_2 = flatten_simd(&columns[14]);
        let o2_3 = flatten_simd(&columns[15]);
        for (((&a, &b), &c), &d) in o2_0.iter().zip(&o2_1).zip(&o2_2).zip(&o2_3) {
            lookup_o2.insert((a, b), (c, d));
        }

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
