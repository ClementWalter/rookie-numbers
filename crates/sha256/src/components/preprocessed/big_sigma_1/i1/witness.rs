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
use utils::{aligned_vec, combine, simd::into_simd, write_col, write_pair};

use crate::{
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma1},
    preprocessed::big_sigma_1::{self, BigSigma1Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    let mut i1_mult = aligned_vec![0u32; 1 << BigSigma1::I1.count_ones()];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.e_i1_low, cols.e_i1_high).for_each(|(e_i1_low, e_i1_high)| {
            let idx_i1 = pext_u32x16(e_i1_low + (e_i1_high << 16), BigSigma1::I1);
            idx_i1
                .to_array()
                .iter()
                .for_each(|x| i1_mult[*x as usize] += 1);
        });
    }

    into_simd(&i1_mult)
        .chunks((1 << (log_size - LOG_N_LANES)) as usize)
        .map(|chunk| chunk.to_vec())
        .collect()
}

pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let big_sigma_1_cols = big_sigma_1::gen_column_simd();
    let BigSigma1Columns {
        i1_low,
        i1_high,
        o1_low,
        o1_high,
        o21_pext,
        ..
    } = BigSigma1Columns::from_slice(&big_sigma_1_cols[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    let i1_den = combine!(
        relations.big_sigma_1.i1,
        [&i1_low, &i1_high, &o1_low, &o1_high, &o21_pext]
    );

    for ([i1_mult_0, i1_mult_1], (i1_den_0, i1_den_1)) in trace
        .array_chunks::<2>()
        .zip(i1_den.chunks(simd_size).tuples())
    {
        write_pair!(
            i1_mult_0
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_den_0.to_vec(),
            i1_mult_1
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_den_1.to_vec(),
            interaction_trace
        );
    }

    if trace.len() % 2 == 1 {
        let i1_mult = trace.last().unwrap();
        let i1_den_chunk = i1_den.chunks(simd_size).last().unwrap();
        write_col!(
            i1_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_den_chunk.to_vec(),
            interaction_trace
        );
    }

    interaction_trace.finalize_last()
}
