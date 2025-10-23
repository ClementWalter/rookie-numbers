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
use utils::{combine, simd_vec, write_pair};

use crate::{
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns,
        preprocessed::ch_left::i0::columns::ComponentColumns, W_SIZE,
    },
    partitions::{pext_u32x16, BigSigma1},
    preprocessed::ch_left::{self, ChLeftColumns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut i0_low_mult = vec![0u32; 1 << (BigSigma1::I0_L.count_ones() * 2)];
    let mut i0_high_mult = vec![0u32; 1 << (BigSigma1::I0_H.count_ones() * 2)];

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.e_i0_low, cols.f_i0_low).for_each(|(e_i0_low, f_i0_low)| {
            let idx_i0_low = (pext_u32x16(*e_i0_low, BigSigma1::I0_L)
                << BigSigma1::I0_L.count_ones())
                + pext_u32x16(*f_i0_low, BigSigma1::I0_L);
            idx_i0_low
                .to_array()
                .iter()
                .for_each(|x| i0_low_mult[*x as usize] += 1);
        });
        izip!(cols.e_i0_high, cols.f_i0_high).for_each(|(e_i0_high, f_i0_high)| {
            let idx_i0_high = (pext_u32x16(*e_i0_high, BigSigma1::I0_H)
                << BigSigma1::I0_H.count_ones())
                + pext_u32x16(*f_i0_high, BigSigma1::I0_H);
            idx_i0_high
                .to_array()
                .iter()
                .for_each(|x| i0_high_mult[*x as usize] += 1);
        });
    }

    simd_vec!(i0_low_mult, i0_high_mult)
}

pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let preprocessed_columns = ch_left::gen_column_simd();
    let ChLeftColumns {
        i0_low_e,
        i0_low_f,
        i0_low_res,
        i0_high_e,
        i0_high_f,
        i0_high_res,
        ..
    } = ChLeftColumns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let i0_low = combine!(relations.ch_left.i0_low, [i0_low_e, i0_low_f, i0_low_res]);
    let i0_high = combine!(
        relations.ch_left.i0_high,
        [i0_high_e, i0_high_f, i0_high_res]
    );

    write_pair!(
        cols.i0_low_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i0_low,
        cols.i0_high_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        i0_high,
        interaction_trace
    );
    interaction_trace.finalize_last()
}
