use std::simd::u32x16;

use itertools::izip;
use stwo::{
    core::{
        fields::{m31::BaseField, qm31::QM31},
        poly::circle::CanonicCoset,
        ColumnVec,
    },
    prover::{
        backend::simd::{column::BaseColumn, m31::LOG_N_LANES, qm31::PackedQM31, SimdBackend},
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::{LogupTraceGenerator, Relation};

use crate::{
    column_vec, combine,
    components::{
        preprocessed::sigma_0_i0_i1::columns::ComponentColumns,
        scheduling::columns::RoundInteractionColumns as SchedulingInteractionColumns, W_SIZE,
    },
    partitions::{pext_u32x16, Sigma0},
    preprocessed::sigma_0::{self, Columns as Sigma0Columns},
    relations::Relations,
    sha256::N_SCHEDULING_ROUNDS,
    write_pair,
};

#[allow(clippy::type_complexity)]
pub fn gen_trace(
    scheduling_lookup_data: &[Vec<u32x16>],
) -> ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    debug_assert_eq!(Sigma0::I0.count_ones(), Sigma0::I1.count_ones());

    // Dense counters for each relation
    let mut sigma_0_i0_mult = vec![0u32; 1 << Sigma0::I0.count_ones()];
    let mut sigma_0_i1_mult = vec![0u32; 1 << Sigma0::I1.count_ones()];

    // Aggregate over all scheduling lookups
    for round in 0..N_SCHEDULING_ROUNDS {
        let start = W_SIZE + round * SchedulingInteractionColumns::SIZE;
        let end = start + SchedulingInteractionColumns::SIZE;

        let cols = SchedulingInteractionColumns::from_slice(&scheduling_lookup_data[start..end]);

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
    let preprocessed_columns = sigma_0::gen_column_simd()
        .into_iter()
        .map(|col| {
            col.data
                .iter()
                .map(|p| p.into_simd())
                .collect::<Vec<u32x16>>()
        })
        .collect::<Vec<Vec<u32x16>>>();
    let Sigma0Columns {
        sigma_0_i0_low,
        sigma_0_i0_high,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext,
        sigma_0_i1_low,
        sigma_0_i1_high,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext,
        ..
    } = Sigma0Columns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].data.len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let sigma_0_i0 = combine!(
        relations.sigma_0.i0,
        sigma_0_i0_low,
        sigma_0_i0_high,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext
    );
    let sigma_0_i1 = combine!(
        relations.sigma_0.i1,
        sigma_0_i1_low,
        sigma_0_i1_high,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext
    );

    write_pair!(
        cols.sigma_0_i0_mult
            .data
            .iter()
            .map(|v| (*v).into())
            .collect::<Vec<PackedQM31>>(),
        sigma_0_i0,
        cols.sigma_0_i1_mult
            .data
            .iter()
            .map(|v| (*v).into())
            .collect::<Vec<PackedQM31>>(),
        sigma_0_i1,
        interaction_trace
    );
    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{partitions::SubsetIterator, preprocessed::sigma_0};

    #[test]
    fn test_value_at_index() {
        let preprocessed_cols: Vec<Vec<u32x16>> = sigma_0::gen_column_simd()
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
