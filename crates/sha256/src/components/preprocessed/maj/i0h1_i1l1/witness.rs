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
        preprocessed::maj::i0h1_i1l1::columns::ComponentColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma0},
    preprocessed::maj::{self, MajColumns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
    simd_vec, write_pair,
};

pub fn gen_trace(
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i0_high_1_mult = vec![0u32; 1 << (BigSigma0::I0_H1.count_ones() * 3)];
    let mut i1_low_1_mult = vec![0u32; 1 << (BigSigma0::I1_L1.count_ones() * 3)];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.a_i0_high_1, cols.b_i0_high_1, cols.c_i0_high_1).for_each(
            |(a_i0_high_1, b_i0_high_1, c_i0_high_1)| {
                let a_pos = pext_u32x16(*a_i0_high_1, BigSigma0::I0_H1)
                    << (BigSigma0::I0_H1.count_ones() * 2);
                let b_pos =
                    pext_u32x16(*b_i0_high_1, BigSigma0::I0_H1) << (BigSigma0::I0_H1.count_ones());
                let c_pos = pext_u32x16(*c_i0_high_1, BigSigma0::I0_H1);
                let idx_i0_high_1 = a_pos + b_pos + c_pos;
                idx_i0_high_1
                    .to_array()
                    .iter()
                    .for_each(|x| i0_high_1_mult[*x as usize] += 1);
            },
        );

        izip!(cols.a_i1_low_1, cols.b_i1_low_1, cols.c_i1_low_1).for_each(
            |(a_i1_low_1, b_i1_low_1, c_i1_low_1)| {
                let a_pos = pext_u32x16(*a_i1_low_1, BigSigma0::I1_L1)
                    << (BigSigma0::I1_L1.count_ones() * 2);
                let b_pos =
                    pext_u32x16(*b_i1_low_1, BigSigma0::I1_L1) << (BigSigma0::I1_L1.count_ones());
                let c_pos = pext_u32x16(*c_i1_low_1, BigSigma0::I1_L1);
                let idx_i1_low_1 = a_pos + b_pos + c_pos;
                idx_i1_low_1
                    .to_array()
                    .iter()
                    .for_each(|x| i1_low_1_mult[*x as usize] += 1);
            },
        );
    }

    simd_vec!(i0_high_1_mult, i1_low_1_mult)
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
        i0_high_1_a,
        i0_high_1_b,
        i0_high_1_c,
        i0_high_1_res,
        i1_low_1_a,
        i1_low_1_b,
        i1_low_1_c,
        i1_low_1_res,
        ..
    } = MajColumns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let i0_high_1 = combine!(
        relations.maj.i0_high_1,
        i0_high_1_a,
        i0_high_1_b,
        i0_high_1_c,
        i0_high_1_res,
    );
    let i1_low_1 = combine!(
        relations.maj.i1_low_1,
        i1_low_1_a,
        i1_low_1_b,
        i1_low_1_c,
        i1_low_1_res
    );

    write_pair!(
        cols.i0_high_1_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i0_high_1,
        cols.i1_low_1_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i1_low_1,
        interaction_trace
    );
    interaction_trace.finalize_last()
}
