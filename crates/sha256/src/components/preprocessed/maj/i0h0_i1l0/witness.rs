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
    preprocessed::maj::{self, MajColumns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i0_high_0_mult = aligned_vec![0u32; 1 << (BigSigma0::I0_H0.count_ones() * 3)];
    let mut i1_low_0_mult = aligned_vec![0u32; 1 << (BigSigma0::I1_L0.count_ones() * 3)];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.a_i0_high_0, cols.b_i0_high_0, cols.c_i0_high_0).for_each(
            |(a_i0_high_0, b_i0_high_0, c_i0_high_0)| {
                let a_pos = pext_u32x16(*a_i0_high_0, BigSigma0::I0_H0)
                    << (BigSigma0::I0_H0.count_ones() * 2);
                let b_pos =
                    pext_u32x16(*b_i0_high_0, BigSigma0::I0_H0) << (BigSigma0::I0_H0.count_ones());
                let c_pos = pext_u32x16(*c_i0_high_0, BigSigma0::I0_H0);
                let idx_i0_high_0 = a_pos + b_pos + c_pos;
                idx_i0_high_0
                    .to_array()
                    .iter()
                    .for_each(|x| i0_high_0_mult[*x as usize] += 1);
            },
        );

        izip!(cols.a_i1_low_0, cols.b_i1_low_0, cols.c_i1_low_0).for_each(
            |(a_i1_low_0, b_i1_low_0, c_i1_low_0)| {
                let a_pos = pext_u32x16(*a_i1_low_0, BigSigma0::I1_L0)
                    << (BigSigma0::I1_L0.count_ones() * 2);
                let b_pos =
                    pext_u32x16(*b_i1_low_0, BigSigma0::I1_L0) << (BigSigma0::I1_L0.count_ones());
                let c_pos = pext_u32x16(*c_i1_low_0, BigSigma0::I1_L0);
                let idx_i1_low_0 = a_pos + b_pos + c_pos;
                idx_i1_low_0
                    .to_array()
                    .iter()
                    .for_each(|x| i1_low_0_mult[*x as usize] += 1);
            },
        );
    }

    into_simd(&i0_high_0_mult)
        .chunks((1 << (log_size - LOG_N_LANES)) as usize)
        .zip_eq(into_simd(&i1_low_0_mult).chunks((1 << (log_size - LOG_N_LANES)) as usize))
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
    let preprocessed_columns = maj::gen_column_simd();
    let MajColumns {
        i0_high_0_a,
        i0_high_0_b,
        i0_high_0_c,
        i0_high_0_res,
        i1_low_0_a,
        i1_low_0_b,
        i1_low_0_c,
        i1_low_0_res,
        ..
    } = MajColumns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    let i0_high_0 = combine!(
        relations.maj.i0_high_0,
        [&i0_high_0_a, &i0_high_0_b, &i0_high_0_c, &i0_high_0_res]
    );
    let i1_low_0 = combine!(
        relations.maj.i1_low_0,
        [&i1_low_0_a, &i1_low_0_b, &i1_low_0_c, &i1_low_0_res]
    );

    for ([i0_high_0_mult, i1_low_0_mult], (i0_high_0_den, i1_low_0_den)) in trace
        .array_chunks::<2>()
        .zip(i0_high_0.chunks(simd_size).zip(i1_low_0.chunks(simd_size)))
    {
        write_pair!(
            i0_high_0_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0_high_0_den.to_vec(),
            i1_low_0_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_low_0_den.to_vec(),
            interaction_trace
        );
    }

    interaction_trace.finalize_last()
}
