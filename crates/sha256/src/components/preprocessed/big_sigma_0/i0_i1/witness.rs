use std::simd::u32x16;

use itertools::{izip, Itertools};
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
use utils::{aligned_vec, combine, simd::into_simd, write_pair};

use crate::{
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma0},
    preprocessed::big_sigma_0::{self, BigSigma0Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i0_mult = aligned_vec![0u32; 1 << BigSigma0::I0.count_ones()];
    let mut i1_mult = aligned_vec![0u32; 1 << BigSigma0::I1.count_ones()];

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

    // Split into chunks of size log_size - LOG_N_LANES
    into_simd(&i0_mult)
        .chunks((1 << (log_size - LOG_N_LANES)) as usize)
        .zip_eq(into_simd(&i1_mult).chunks((1 << (log_size - LOG_N_LANES)) as usize))
        .flat_map(|(i0, i1)| [i0.to_vec(), i1.to_vec()])
        .collect()
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
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    for (i, [i0_mult, i1_mult]) in trace.array_chunks::<2>().enumerate() {
        let start = i * simd_size;
        let end = start + simd_size;

        let i0 = combine!(
            relations.big_sigma_0.i0,
            [
                &i0_low[start..end],
                &i0_high_0[start..end],
                &i0_high_1[start..end],
                &o0_low[start..end],
                &o0_high[start..end],
                &o20_pext[start..end]
            ]
        );
        let i1 = combine!(
            relations.big_sigma_0.i1,
            [
                &i1_low_0[start..end],
                &i1_low_1[start..end],
                &i1_high[start..end],
                &o1_low[start..end],
                &o1_high[start..end],
                &o21_pext[start..end]
            ]
        );
        write_pair!(
            i0_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0,
            i1_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1,
            interaction_trace
        );
    }
    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        components::{
            compression::witness::gen_trace as gen_compression_trace,
            preprocessed::big_sigma_0::i0_i1::witness::gen_trace,
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        partitions::SubsetIterator,
        preprocessed::big_sigma_0,
    };

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

    #[test_log::test]
    fn test_trace() {
        const LOG_N_SHA256: u32 = 8;
        let log_size = 21;
        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_SHA256);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(log_size, &scheduling_lookup_data, &compression_lookup_data);
        assert!(trace.iter().all(|t| t.len() == trace[0].len()));
    }
}
