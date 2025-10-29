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
    partitions::{pext_u32x16, BigSigma1},
    preprocessed::ch_right::{self, ChRightColumns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i1_low_mult = aligned_vec![0u32; 1 << (BigSigma1::I1_L.count_ones() * 2)];
    let mut i1_high_mult = aligned_vec![0u32; 1 << (BigSigma1::I1_H.count_ones() * 2)];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.e_i1_low, cols.g_i1_low).for_each(|(e_i1_low, g_i1_low)| {
            let idx_i1_low = (pext_u32x16(*e_i1_low, BigSigma1::I1_L)
                << BigSigma1::I1_L.count_ones())
                + pext_u32x16(*g_i1_low, BigSigma1::I1_L);
            idx_i1_low
                .to_array()
                .iter()
                .for_each(|x| i1_low_mult[*x as usize] += 1);
        });
        izip!(cols.e_i1_high, cols.g_i1_high).for_each(|(e_i1_high, g_i1_high)| {
            let idx_i1_high = (pext_u32x16(*e_i1_high, BigSigma1::I1_H)
                << BigSigma1::I1_H.count_ones())
                + pext_u32x16(*g_i1_high, BigSigma1::I1_H);
            idx_i1_high
                .to_array()
                .iter()
                .for_each(|x| i1_high_mult[*x as usize] += 1);
        });
    }

    into_simd(&i1_low_mult)
        .chunks((1 << (log_size - LOG_N_LANES)) as usize)
        .zip_eq(into_simd(&i1_high_mult).chunks((1 << (log_size - LOG_N_LANES)) as usize))
        .flat_map(|(i1_low, i1_high)| [i1_low.to_vec(), i1_high.to_vec()])
        .collect()
}

pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let preprocessed_columns = ch_right::gen_column_simd();
    let ChRightColumns {
        i1_low_e,
        i1_low_g,
        i1_low_res,
        i1_high_e,
        i1_high_g,
        i1_high_res,
        ..
    } = ChRightColumns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    for (i, [i1_low_mult, i1_high_mult]) in trace.array_chunks::<2>().enumerate() {
        let start = i * simd_size;
        let end = start + simd_size;

        let i1_low = combine!(
            relations.ch_right.i1_low,
            [
                &i1_low_e[start..end],
                &i1_low_g[start..end],
                &i1_low_res[start..end]
            ]
        );
        let i1_high = combine!(
            relations.ch_right.i1_high,
            [
                &i1_high_e[start..end],
                &i1_high_g[start..end],
                &i1_high_res[start..end]
            ]
        );
        write_pair!(
            i1_low_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_low,
            i1_high_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_high,
            interaction_trace
        );
    }
    interaction_trace.finalize_last()
}
