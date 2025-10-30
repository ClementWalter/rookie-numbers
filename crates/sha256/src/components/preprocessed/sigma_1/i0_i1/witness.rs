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
        scheduling::columns::RoundInteractionColumns as SchedulingInteractionColumns, W_SIZE,
    },
    partitions::{pext_u32x16, Sigma1},
    preprocessed::sigma_1::{self, Sigma1Columns},
    relations::Relations,
    sha256::N_SCHEDULING_ROUNDS,
};

pub fn gen_trace(
    log_size: u32,
    scheduling_lookup_data: &[Vec<u32x16>],
    _compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    debug_assert_eq!(Sigma1::I0.count_ones(), Sigma1::I1.count_ones());

    // Dense counters for each relation
    let mut sigma_1_i0_mult = aligned_vec![0u32; 1 << Sigma1::I0.count_ones()];
    let mut sigma_1_i1_mult = aligned_vec![0u32; 1 << Sigma1::I1.count_ones()];

    // Aggregate over all scheduling lookups
    for round in 0..N_SCHEDULING_ROUNDS {
        let start = W_SIZE + round * SchedulingInteractionColumns::SIZE;
        let end = start + SchedulingInteractionColumns::SIZE;

        let cols = SchedulingInteractionColumns::from_slice(&scheduling_lookup_data[start..end]);

        izip!(cols.w_2_i0_low, cols.w_2_i0_high).for_each(|(w_2_i0_low, w_2_i0_high)| {
            let idx_i0 = pext_u32x16(w_2_i0_low + (w_2_i0_high << 16), Sigma1::I0);
            idx_i0
                .to_array()
                .iter()
                .for_each(|x| sigma_1_i0_mult[*x as usize] += 1);
        });
        izip!(cols.w_2_i1_low, cols.w_2_i1_high).for_each(|(w_2_i1_low, w_2_i1_high)| {
            let idx_i1 = pext_u32x16(w_2_i1_low + (w_2_i1_high << 16), Sigma1::I1);
            idx_i1
                .to_array()
                .iter()
                .for_each(|x| sigma_1_i1_mult[*x as usize] += 1);
        });
    }

    into_simd(&sigma_1_i0_mult)
        .chunks((1 << (log_size - LOG_N_LANES)) as usize)
        .zip_eq(into_simd(&sigma_1_i1_mult).chunks((1 << (log_size - LOG_N_LANES)) as usize))
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
    let preprocessed_columns = sigma_1::gen_column_simd();
    let Sigma1Columns {
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext,
        i1_low,
        i1_high,
        o1_low,
        o1_high,
        o21_pext,
        ..
    } = Sigma1Columns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let log_size = simd_size.ilog2() + LOG_N_LANES;
    let mut interaction_trace = LogupTraceGenerator::new(log_size);

    let i0 = combine!(
        relations.sigma_1.i0,
        [&i0_low, &i0_high, &o0_low, &o0_high, &o20_pext]
    );
    let i1 = combine!(
        relations.sigma_1.i1,
        [&i1_low, &i1_high, &o1_low, &o1_high, &o21_pext]
    );

    for ([sigma_1_i0_mult, sigma_1_i1_mult], (i0_den, i1_den)) in trace
        .array_chunks::<2>()
        .zip(i0.chunks(simd_size).zip(i1.chunks(simd_size)))
    {
        write_pair!(
            sigma_1_i0_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i0_den.to_vec(),
            sigma_1_i1_mult
                .iter()
                .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
                .map(PackedQM31::from),
            i1_den.to_vec(),
            interaction_trace
        );
    }

    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{partitions::SubsetIterator, preprocessed::sigma_1};

    #[test]
    fn test_value_at_index() {
        let preprocessed_cols: Vec<Vec<u32x16>> = sigma_1::gen_column_simd();

        let mut iterator = SubsetIterator::new(Sigma1::I0);

        let x: [u32; 16] = std::array::from_fn(|_| iterator.next().unwrap());
        let x_low = u32x16::from_slice(&x.map(|x| x & Sigma1::I0_L));
        let x_high = u32x16::from_slice(&x.map(|x| (x >> 16) & Sigma1::I0_H));
        let idx_i0 = pext_u32x16(x_low + (x_high << 16), Sigma1::I0);

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
            x_high.to_array().to_vec()
        );
    }
}
