use std::simd::u32x16;

use itertools::izip;
use stwo::{
    core::{
        fields::{m31::BaseField, qm31::QM31},
        ColumnVec,
    },
    prover::{
        backend::simd::{
            m31::{PackedM31, LOG_N_LANES},
            qm31::PackedQM31,
            SimdBackend,
        },
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::{LogupTraceGenerator, Relation};

use crate::{
    combine,
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns,
        preprocessed::big_sigma_0::i0_i1::columns::ComponentColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma0},
    preprocessed::big_sigma_0::{self, BigSigma0Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
    simd_vec, write_pair,
};

pub fn gen_trace(
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i0_mult = vec![0u32; 1 << BigSigma0::I0.count_ones()];
    let mut i1_mult = vec![0u32; 1 << BigSigma0::I1.count_ones()];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.a_i0_low, cols.a_i0_high_0, cols.a_i0_high_1).for_each(
            |(a_i0_low, a_i0_high_0, a_i0_high_1)| {
                let idx_i0 = pext_u32x16(
                    a_i0_low + (a_i0_high_0 << 16) + (a_i0_high_1 << 24),
                    BigSigma0::I0,
                );
                idx_i0
                    .to_array()
                    .iter()
                    .for_each(|x| i0_mult[*x as usize] += 1);
            },
        );
        izip!(cols.a_i1_low_0, cols.a_i1_low_1, cols.a_i1_high).for_each(
            |(a_i1_low_0, a_i1_low_1, a_i1_high)| {
                let idx_i1 = pext_u32x16(
                    a_i1_low_0 + (a_i1_low_1 << 8) + (a_i1_high << 16),
                    BigSigma0::I1,
                );
                idx_i1
                    .to_array()
                    .iter()
                    .for_each(|x| i1_mult[*x as usize] += 1);
            },
        );
    }

    simd_vec!(i0_mult, i1_mult)
}

pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let preprocessed_columns = big_sigma_0::gen_column_simd();
    let BigSigma0Columns {
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
        ..
    } = BigSigma0Columns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let i0 = combine!(
        relations.big_sigma_0.i0,
        i0_low,
        i0_high_0,
        i0_high_1,
        o0_low,
        o0_high,
        o20_pext
    );
    let i1 = combine!(
        relations.big_sigma_0.i1,
        i1_low_0,
        i1_low_1,
        i1_high,
        o1_low,
        o1_high,
        o21_pext
    );

    write_pair!(
        cols.i0_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i0,
        cols.i1_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i1,
        interaction_trace
    );
    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{partitions::SubsetIterator, preprocessed::big_sigma_0};

    #[test]
    fn test_value_at_index() {
        let preprocessed_cols: Vec<Vec<u32x16>> = big_sigma_0::gen_column_simd();

        let mut iterator = SubsetIterator::new(BigSigma0::I0);

        let x: [u32; 16] = std::array::from_fn(|_| iterator.next().unwrap());
        let x_low = u32x16::from_slice(&x.map(|x| x & BigSigma0::I0_L));
        let x_high_0 = u32x16::from_slice(&x.map(|x| (x >> 16) & BigSigma0::I0_H0));
        let x_high_1 = u32x16::from_slice(&x.map(|x| (x >> 24) & BigSigma0::I0_H1));
        let idx_i0 = pext_u32x16(x_low + (x_high_0 << 16) + (x_high_1 << 32), BigSigma0::I0);

        assert_eq!(
            idx_i0
                .to_array()
                .iter()
                .map(|x| {
                    preprocessed_cols[0][(*x >> LOG_N_LANES) as usize].to_array()
                        [(*x % (1 << LOG_N_LANES)) as usize]
                })
                .collect::<Vec<u32>>(),
            x_low.to_array().to_vec()
        );
        assert_eq!(
            idx_i0
                .to_array()
                .iter()
                .map(|x| {
                    preprocessed_cols[1][(*x >> LOG_N_LANES) as usize].to_array()
                        [(*x % (1 << LOG_N_LANES)) as usize]
                })
                .collect::<Vec<u32>>(),
            x_high_0.to_array().to_vec()
        );
        assert_eq!(
            idx_i0
                .to_array()
                .iter()
                .map(|x| {
                    preprocessed_cols[2][(*x >> LOG_N_LANES) as usize].to_array()
                        [(*x % (1 << LOG_N_LANES)) as usize]
                })
                .collect::<Vec<u32>>(),
            x_high_1.to_array().to_vec()
        );
    }
}
