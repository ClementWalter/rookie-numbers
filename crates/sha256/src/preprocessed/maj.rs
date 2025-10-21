use std::simd::u32x16;

use itertools::{iproduct, Itertools};
use stwo::core::channel::Channel;
use stwo_constraint_framework::relation;

use crate::{
    partitions::{BigSigma0 as BigSigma0Partitions, SubsetIterator},
    sha256::maj_u32x16,
    trace_columns,
};

// [a, b, c, val]
const N_COLUMNS: usize = 4;

relation!(MAJ_I0_L, N_COLUMNS);
relation!(MAJ_I0_H0, N_COLUMNS);
relation!(MAJ_I0_H1, N_COLUMNS);
relation!(MAJ_I1_L0, N_COLUMNS);
relation!(MAJ_I1_L1, N_COLUMNS);
relation!(MAJ_I1_H, N_COLUMNS);

trace_columns!(
    MajColumns,
    i0_low_a,
    i0_low_b,
    i0_low_c,
    i0_low_res,
    i1_high_a,
    i1_high_b,
    i1_high_c,
    i1_high_res,
    i0_high_0_a,
    i0_high_0_b,
    i0_high_0_c,
    i0_high_0_res,
    i1_low_0_a,
    i1_low_0_b,
    i1_low_0_c,
    i1_low_0_res,
    i0_high_1_a,
    i0_high_1_b,
    i0_high_1_c,
    i0_high_1_res,
    i1_low_1_a,
    i1_low_1_b,
    i1_low_1_c,
    i1_low_1_res,
);

trace_columns!(
    MajI0LI1HColumns,
    i0_low_a,
    i0_low_b,
    i0_low_c,
    i0_low_res,
    i1_high_a,
    i1_high_b,
    i1_high_c,
    i1_high_res,
);

trace_columns!(
    MajI0H0I1L0Columns,
    i0_high_0_a,
    i0_high_0_b,
    i0_high_0_c,
    i0_high_0_res,
    i1_low_0_a,
    i1_low_0_b,
    i1_low_0_c,
    i1_low_0_res,
);

trace_columns!(
    MajI0H1I1L1Columns,
    i0_high_1_a,
    i0_high_1_b,
    i0_high_1_c,
    i0_high_1_res,
    i1_low_1_a,
    i1_low_1_b,
    i1_low_1_c,
    i1_low_1_res,
);

#[derive(Debug, Clone)]
pub struct Relation {
    pub i0_low: MAJ_I0_L,
    pub i0_high_0: MAJ_I0_H0,
    pub i0_high_1: MAJ_I0_H1,
    pub i1_low_0: MAJ_I1_L0,
    pub i1_low_1: MAJ_I1_L1,
    pub i1_high: MAJ_I1_H,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            i0_low: MAJ_I0_L::dummy(),
            i0_high_0: MAJ_I0_H0::dummy(),
            i0_high_1: MAJ_I0_H1::dummy(),
            i1_low_0: MAJ_I1_L0::dummy(),
            i1_low_1: MAJ_I1_L1::dummy(),
            i1_high: MAJ_I1_H::dummy(),
        }
    }

    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            i0_low: MAJ_I0_L::draw(channel),
            i0_high_0: MAJ_I0_H0::draw(channel),
            i0_high_1: MAJ_I0_H1::draw(channel),
            i1_low_0: MAJ_I1_L0::draw(channel),
            i1_low_1: MAJ_I1_L1::draw(channel),
            i1_high: MAJ_I1_H::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<Vec<u32x16>> {
    let mut all_columns: Vec<Vec<u32x16>> = vec![Vec::new(); MajColumns::SIZE];

    for (i, partition) in [
        BigSigma0Partitions::I0_L,
        BigSigma0Partitions::I1_H,
        BigSigma0Partitions::I0_H0,
        BigSigma0Partitions::I1_L0,
        BigSigma0Partitions::I0_H1,
        BigSigma0Partitions::I1_L1,
    ]
    .iter()
    .enumerate()
    {
        // lookup
        let tuples: Vec<(u32x16, u32x16, u32x16, u32x16)> = iproduct!(
            SubsetIterator::new(*partition),
            SubsetIterator::new(*partition),
            SubsetIterator::new(*partition)
        )
        .chunks(16)
        .into_iter()
        .map(|chunk| {
            let mut xs = [0u32; 16];
            let mut ys = [0u32; 16];
            let mut zs = [0u32; 16];
            for (i, (x, y, z)) in chunk.enumerate() {
                xs[i] = x;
                ys[i] = y;
                zs[i] = z;
            }
            let a = u32x16::from_array(xs);
            let b = u32x16::from_array(ys);
            let c = u32x16::from_array(zs);
            (a, b, c, maj_u32x16(a, b, c))
        })
        .collect();

        for (a, b, c, res) in tuples {
            all_columns[4 * i].push(a);
            all_columns[4 * i + 1].push(b);
            all_columns[4 * i + 2].push(c);
            all_columns[4 * i + 3].push(res);
        }
    }

    all_columns
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::sha256::maj;

    #[test]
    fn test_ids() {
        assert_eq!(
            MajColumns::SIZE,
            MajI0LI1HColumns::SIZE + MajI0H0I1L0Columns::SIZE + MajI0H1I1L1Columns::SIZE
        );
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

        let mut lookup_i0_l: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i0_l_0 = flatten_simd_column(&columns[0]);
        let i0_l_1 = flatten_simd_column(&columns[1]);
        let i0_l_2 = flatten_simd_column(&columns[2]);
        let i0_l_3 = flatten_simd_column(&columns[3]);
        for (((&a, &b), &c), &d) in i0_l_0.iter().zip(&i0_l_1).zip(&i0_l_2).zip(&i0_l_3) {
            lookup_i0_l.insert((a, b, c), d);
        }

        let mut lookup_i1_h: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i1_h_0 = flatten_simd_column(&columns[4]);
        let i1_h_1 = flatten_simd_column(&columns[5]);
        let i1_h_2 = flatten_simd_column(&columns[6]);
        let i1_h_3 = flatten_simd_column(&columns[7]);
        for (((&a, &b), &c), &d) in i1_h_0.iter().zip(&i1_h_1).zip(&i1_h_2).zip(&i1_h_3) {
            lookup_i1_h.insert((a, b, c), d);
        }

        let mut lookup_i0_h_0: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i0_h_0_0 = flatten_simd_column(&columns[8]);
        let i0_h_0_1 = flatten_simd_column(&columns[9]);
        let i0_h_0_2 = flatten_simd_column(&columns[10]);
        let i0_h_0_3 = flatten_simd_column(&columns[11]);
        for (((&a, &b), &c), &d) in i0_h_0_0.iter().zip(&i0_h_0_1).zip(&i0_h_0_2).zip(&i0_h_0_3) {
            lookup_i0_h_0.insert((a, b, c), d);
        }

        let mut lookup_i1_l_0: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i1_l_0_0 = flatten_simd_column(&columns[12]);
        let i1_l_0_1 = flatten_simd_column(&columns[13]);
        let i1_l_0_2 = flatten_simd_column(&columns[14]);
        let i1_l_0_3 = flatten_simd_column(&columns[15]);
        for (((&a, &b), &c), &d) in i1_l_0_0.iter().zip(&i1_l_0_1).zip(&i1_l_0_2).zip(&i1_l_0_3) {
            lookup_i1_l_0.insert((a, b, c), d);
        }

        let mut lookup_i0_h_1: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i0_h_1_0 = flatten_simd_column(&columns[16]);
        let i0_h_1_1 = flatten_simd_column(&columns[17]);
        let i0_h_1_2 = flatten_simd_column(&columns[18]);
        let i0_h_1_3 = flatten_simd_column(&columns[19]);
        for (((&a, &b), &c), &d) in i0_h_1_0.iter().zip(&i0_h_1_1).zip(&i0_h_1_2).zip(&i0_h_1_3) {
            lookup_i0_h_1.insert((a, b, c), d);
        }

        let mut lookup_i1_l_1: HashMap<(u32, u32, u32), u32> = HashMap::new();
        let i1_l_1_0 = flatten_simd_column(&columns[20]);
        let i1_l_1_1 = flatten_simd_column(&columns[21]);
        let i1_l_1_2 = flatten_simd_column(&columns[22]);
        let i1_l_1_3 = flatten_simd_column(&columns[23]);
        for (((&a, &b), &c), &d) in i1_l_1_0.iter().zip(&i1_l_1_1).zip(&i1_l_1_2).zip(&i1_l_1_3) {
            lookup_i1_l_1.insert((a, b, c), d);
        }

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
