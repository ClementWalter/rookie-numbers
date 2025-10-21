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
        preprocessed::big_sigma_1::i0::columns::ComponentColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma1},
    preprocessed::big_sigma_1::{self, BigSigma1Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
    simd_vec, write_col,
};

pub fn gen_trace(
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
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext
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
