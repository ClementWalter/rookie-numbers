use std::simd::u32x16;

use itertools::{izip, Itertools};
use stwo::{
    core::{
        fields::{m31::BaseField, qm31::QM31},
        poly::circle::CanonicCoset,
        ColumnVec,
    },
    prover::{
        backend::simd::{
            column::BaseColumn,
            m31::{PackedM31, LOG_N_LANES},
            SimdBackend,
        },
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::{LogupTraceGenerator, Relation};

use crate::{
    column_vec, combine,
    components::{
        scheduling::{
            columns::InteractionColumns as SchedulingInteractionColumns,
            witness::N_INTERACTION_COLUMNS,
        },
        W_SIZE,
    },
    partitions::{pext_u32x16, Sigma0},
    preprocessed::{
        sigma_0::{Columns, InteractionColumnsIndex},
        PreProcessedColumn,
    },
    relations::Relations,
    write_col, write_pair,
};

#[allow(clippy::type_complexity)]
pub fn gen_trace(
    scheduling: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    debug_assert_eq!(Sigma0::I0.count_ones(), Sigma0::I1.count_ones());

    // Dense counters for each relation
    let mut sigma_0_i0_mult = vec![0u32; 1 << Sigma0::I0.count_ones()];
    let mut sigma_0_i1_mult = vec![0u32; 1 << Sigma0::I1.count_ones()];

    // Aggregate over all scheduling lookups
    for round in scheduling
        .iter()
        .skip(W_SIZE)
        .map(|col| col.data.iter().map(|x| x.into_simd()))
        .chunks(N_INTERACTION_COLUMNS)
        .into_iter()
    {
        let cols = SchedulingInteractionColumns::from_iter(round);

        izip!(cols.w_15_i0_low, cols.w_15_i0_high).for_each(|(w_15_i0_low, w_15_i0_high)| {
            let idx_i0 = pext_u32x16(w_15_i0_low + (w_15_i0_high << 16), Sigma0::I0);
            idx_i0
                .to_array()
                .iter()
                .for_each(|x| sigma_0_i0_mult[*x as usize] += 1);
        });
        izip!(cols.w_15_i1_low, cols.w_15_i1_high).for_each(|(w_15_i1_low, w_15_i1_high)| {
            let idx_i1 = pext_u32x16(w_15_i1_low + (w_15_i1_high << 16), Sigma0::I1);
            idx_i1
                .to_array()
                .iter()
                .for_each(|x| sigma_0_i1_mult[*x as usize] += 1);
        });
    }

    column_vec!(sigma_0_i0_mult, sigma_0_i1_mult)
}

pub fn gen_interaction_trace(
    trace: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let preprocessed_trace: Vec<Vec<u32x16>> = Columns
        .gen_column_simd()
        .into_iter()
        .map(|col| col.data.iter().map(|p| p.into_simd()).collect())
        .collect();
    let trace: Vec<&Vec<PackedM31>> = trace.iter().map(|col| &col.data).collect();

    let sigma_0_i0 = combine!(
        relations.sigma_0.i0,
        preprocessed_trace,
        0,
        sigma_0_i0_low,
        sigma_0_i0_high,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext
    );
    let sigma_0_i1 = combine!(
        relations.sigma_0.i1,
        preprocessed_trace,
        0,
        sigma_0_i1_low,
        sigma_0_i1_high,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext
    );

    debug_assert_eq!(trace[0].len(), trace[1].len());
    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);
    write_pair!(
        trace[0],
        sigma_0_i0,
        trace[1],
        sigma_0_i1,
        interaction_trace
    );
    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partitions::SubsetIterator;

    #[test]
    fn test_value_at_index() {
        let preprocessed_trace: Vec<Vec<u32x16>> = Columns
            .gen_column_simd()
            .into_iter()
            .map(|col| col.data.iter().map(|p| p.into_simd()).collect())
            .collect();

        let mut iterator = SubsetIterator::new(Sigma0::I0);

        let x: [u32; 16] = std::array::from_fn(|_| iterator.next().unwrap());
        let x_low = u32x16::from_slice(&x.map(|x| x & Sigma0::I0_L));
        let x_high = u32x16::from_slice(&x.map(|x| (x >> 16) & Sigma0::I0_H));
        let idx_i0 = pext_u32x16(x_low + (x_high << 16), Sigma0::I0);

        assert_eq!(
            idx_i0
                .to_array()
                .iter()
                .map(|x| {
                    preprocessed_trace[0][(*x >> LOG_N_LANES) as usize].to_array()
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
                    preprocessed_trace[1][(*x >> LOG_N_LANES) as usize].to_array()
                        [(*x % (1 << LOG_N_LANES)) as usize]
                })
                .collect::<Vec<u32>>(),
            x_high.to_array().to_vec()
        );
    }
}
