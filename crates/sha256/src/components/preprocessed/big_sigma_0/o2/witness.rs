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
    partitions::BigSigma0,
    preprocessed::big_sigma_0::{self, BigSigma0Columns},
    relations::Relations,
    sha256::N_COMPRESSION_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    _scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    let mut o2_mult = aligned_vec![0u32; 1 << (BigSigma0::O2.count_ones() * 2)];

    // Aggregate over all scheduling lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.sigma_0_o20_pext, cols.sigma_0_o21_pext).for_each(
            |(sigma_0_o20_pext, sigma_0_o21_pext)| {
                let idx_o2 = (sigma_0_o20_pext << BigSigma0::O2.count_ones()) + sigma_0_o21_pext;
                idx_o2
                    .to_array()
                    .iter()
                    .for_each(|x| o2_mult[*x as usize] += 1);
            },
        );
    }

    into_simd(&o2_mult)
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
    let preprocessed_columns = big_sigma_0::gen_column_simd();
    let BigSigma0Columns {
        o2_0,
        o2_1,
        o2_low,
        o2_high,
        ..
    } = BigSigma0Columns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    let o2_den = combine!(relations.big_sigma_0.o2, [&o2_0, &o2_1, &o2_low, &o2_high]);

    for ([o2_mult_0, o2_mult_1], (o2_den_0, o2_den_1)) in trace
        .array_chunks::<2>()
        .zip(o2_den.chunks(simd_size).tuple_windows())
    {
        write_pair!(
            o2_mult_0
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            o2_den_0.to_vec(),
            o2_mult_1
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            o2_den_1.to_vec(),
            interaction_trace
        );
    }

    if trace.len() % 2 == 1 {
        let o2_mult = trace.last().unwrap();
        let o2_den_chunk = o2_den.chunks(simd_size).last().unwrap();
        write_col!(
            o2_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            o2_den_chunk.to_vec(),
            interaction_trace
        );
    }

    interaction_trace.finalize_last()
}
