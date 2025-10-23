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
            SimdBackend,
        },
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::{LogupTraceGenerator, Relation};
use utils::{combine, simd_vec, write_col};

use crate::{
    components::{
        compression::columns::RoundInteractionColumns as CompressionInteractionColumns,
        preprocessed::big_sigma_1::o2::columns::ComponentColumns, W_SIZE,
    },
    partitions::BigSigma1,
    preprocessed::big_sigma_1::{self, BigSigma1Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    let mut o2_mult = vec![0u32; 1 << (BigSigma1::O2.count_ones() * 2)];

    // Aggregate over all scheduling lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.sigma_1_o20_pext, cols.sigma_1_o21_pext).for_each(
            |(sigma_1_o20_pext, sigma_1_o21_pext)| {
                let idx_o2 = (sigma_1_o20_pext << BigSigma1::O2.count_ones()) + sigma_1_o21_pext;
                idx_o2
                    .to_array()
                    .iter()
                    .for_each(|x| o2_mult[*x as usize] += 1);
            },
        );
    }

    simd_vec!(o2_mult)
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
        o2_0,
        o2_1,
        o2_low,
        o2_high,
        ..
    } = BigSigma1Columns::from_slice(&big_sigma_1_cols[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let big_sigma_1_o2 = combine!(relations.big_sigma_1.o2, [o2_0, o2_1, o2_low, o2_high]);

    write_col!(
        cols.o2_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(|v| v.into()),
        big_sigma_1_o2,
        interaction_trace
    );
    interaction_trace.finalize_last()
}
