use std::simd::u32x16;

use itertools::izip;
#[cfg(feature = "dynamic-preprocessed-shape")]
use itertools::Itertools;
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
#[cfg(feature = "dynamic-preprocessed-shape")]
use utils::{aligned_vec, combine, simd::into_simd, write_col, write_pair};
#[cfg(not(feature = "dynamic-preprocessed-shape"))]
use utils::{combine, simd_vec, write_col};

use crate::{
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma1},
    preprocessed::big_sigma_1::{self, BigSigma1Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};
#[cfg(feature = "dynamic-preprocessed-shape")]
use crate::preprocessed_log_size;
#[cfg(not(feature = "dynamic-preprocessed-shape"))]
use crate::components::preprocessed::big_sigma_1::i0::columns::ComponentColumns;

#[cfg(feature = "dynamic-preprocessed-shape")]
pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    let mut i0_mult = aligned_vec![0u32; 1 << BigSigma1::I0.count_ones()];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.e_i0_low, cols.e_i0_high).for_each(|(e_i0_low, e_i0_high)| {
            let idx_i0 = pext_u32x16(e_i0_low + (e_i0_high << 16), BigSigma1::I0);
            idx_i0
                .to_array()
                .iter()
                .for_each(|x| i0_mult[*x as usize] += 1);
        });
    }

    let effective_log_size = preprocessed_log_size(log_size);
    into_simd(&i0_mult)
        .chunks((1 << (effective_log_size - LOG_N_LANES)) as usize)
        .map(|chunk| chunk.to_vec())
        .collect()
}

#[cfg(not(feature = "dynamic-preprocessed-shape"))]
pub fn gen_trace(
    _log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    let mut i0_mult = vec![0u32; 1 << BigSigma1::I0.count_ones()];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.e_i0_low, cols.e_i0_high).for_each(|(e_i0_low, e_i0_high)| {
            let idx_i0 = pext_u32x16(e_i0_low + (e_i0_high << 16), BigSigma1::I0);
            idx_i0
                .to_array()
                .iter()
                .for_each(|x| i0_mult[*x as usize] += 1);
        });
    }

    simd_vec!(i0_mult)
}

#[cfg(feature = "dynamic-preprocessed-shape")]
pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let big_sigma_1_cols = big_sigma_1::gen_column_simd();
    let BigSigma1Columns {
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext,
        ..
    } = BigSigma1Columns::from_slice(&big_sigma_1_cols[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    let i0_den = combine!(
        relations.big_sigma_1.i0,
        [&i0_low, &i0_high, &o0_low, &o0_high, &o20_pext]
    );

    for ([i0_mult_0, i0_mult_1], (i0_den_0, i0_den_1)) in trace
        .array_chunks::<2>()
        .zip(i0_den.chunks(simd_size).tuples())
    {
        write_pair!(
            i0_mult_0
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0_den_0.to_vec(),
            i0_mult_1
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0_den_1.to_vec(),
            interaction_trace
        );
    }

    if trace.len() % 2 == 1 {
        let i0_mult = trace.last().unwrap();
        let i0_den_chunk = i0_den.chunks(simd_size).last().unwrap();
        write_col!(
            i0_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0_den_chunk.to_vec(),
            interaction_trace
        );
    }

    interaction_trace.finalize_last()
}

#[cfg(not(feature = "dynamic-preprocessed-shape"))]
pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let big_sigma_1_cols = big_sigma_1::gen_column_simd();
    let BigSigma1Columns {
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext,
        ..
    } = BigSigma1Columns::from_slice(&big_sigma_1_cols[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let i0 = combine!(
        relations.big_sigma_1.i0,
        [i0_low, i0_high, o0_low, o0_high, o20_pext]
    );

    write_col!(
        cols.i0_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i0,
        interaction_trace
    );
    interaction_trace.finalize_last()
}
