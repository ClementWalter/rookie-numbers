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
        preprocessed::range_check_add::columns::ComponentColumns,
        scheduling::columns::RoundInteractionColumns as SchedulingInteractionColumns, W_SIZE,
    },
    preprocessed::range_check_add::{self, RangeCheckAddColumns},
    relations::Relations,
    sha256::{N_COMPRESSION_ROUNDS, N_SCHEDULING_ROUNDS},
    simd_vec, write_col, write_pair,
};

pub fn gen_trace(
    scheduling_lookup_data: &[Vec<u32x16>],
    compression_lookup_data: &[Vec<u32x16>],
) -> Vec<Vec<u32x16>> {
    // Dense counters for each relation
    let mut carry_4_mult = vec![0u32; 1 << 19];
    let mut carry_7_mult = vec![0u32; 1 << 19];
    let mut carry_8_mult = vec![0u32; 1 << 19];

    // Aggregate over all scheduling lookups
    for round in 0..N_SCHEDULING_ROUNDS {
        let start = W_SIZE + round * SchedulingInteractionColumns::SIZE;
        let end = start + SchedulingInteractionColumns::SIZE;

        let cols = SchedulingInteractionColumns::from_slice(&scheduling_lookup_data[start..end]);

        izip!(cols.new_w_low, cols.carry_low).for_each(|(new_w_low, carry_low)| {
            let idx = (new_w_low << 3) + carry_low;
            idx.to_array()
                .iter()
                .for_each(|x| carry_4_mult[*x as usize] += 1);
        });
        izip!(cols.new_w_high, cols.carry_high).for_each(|(new_w_high, carry_high)| {
            let idx = (new_w_high << 3) + carry_high;
            idx.to_array()
                .iter()
                .for_each(|x| carry_4_mult[*x as usize] += 1);
        });
    }

    // Aggregate over all compression lookups
    for round in 0..N_COMPRESSION_ROUNDS {
        let start = W_SIZE + round * CompressionInteractionColumns::SIZE;
        let end = start + CompressionInteractionColumns::SIZE;

        let cols = CompressionInteractionColumns::from_slice(&compression_lookup_data[start..end]);

        izip!(cols.new_e_low, cols.e_carry_low).for_each(|(new_e_low, e_carry_low)| {
            let idx = (new_e_low << 3) + e_carry_low;
            idx.to_array()
                .iter()
                .for_each(|x| carry_7_mult[*x as usize] += 1);
        });
        izip!(cols.new_e_high, cols.e_carry_high).for_each(|(new_e_high, e_carry_high)| {
            let idx = (new_e_high << 3) + e_carry_high;
            idx.to_array()
                .iter()
                .for_each(|x| carry_7_mult[*x as usize] += 1);
        });
        izip!(cols.new_a_low, cols.a_carry_low).for_each(|(new_a_low, a_carry_low)| {
            let idx = (new_a_low << 3) + a_carry_low;
            idx.to_array()
                .iter()
                .for_each(|x| carry_8_mult[*x as usize] += 1);
        });
        izip!(cols.new_a_high, cols.a_carry_high).for_each(|(new_a_high, a_carry_high)| {
            let idx = (new_a_high << 3) + a_carry_high;
            idx.to_array()
                .iter()
                .for_each(|x| carry_8_mult[*x as usize] += 1);
        });
    }

    simd_vec!(carry_4_mult, carry_7_mult, carry_8_mult)
}

pub fn gen_interaction_trace(
    trace: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let preprocessed_columns = range_check_add::gen_column_simd();
    let RangeCheckAddColumns {
        value,
        carry_4,
        carry_7,
        carry_8,
    } = RangeCheckAddColumns::from_slice(&preprocessed_columns[..]);

    let simd_size = trace[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    let cols = ComponentColumns::from_slice(trace);

    let carry_4 = combine!(relations.range_check_add.add_4, value, carry_4);
    let carry_7 = combine!(relations.range_check_add.add_7, value, carry_7);
    let carry_8 = combine!(relations.range_check_add.add_8, value, carry_8);

    write_pair!(
        cols.carry_4_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        carry_4,
        cols.carry_7_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        carry_7,
        interaction_trace
    );
    write_col!(
        cols.carry_8_mult
            .iter()
            .map(|v| unsafe { PackedM31::from_simd_unchecked(*v) })
            .map(PackedQM31::from),
        carry_8,
        interaction_trace
    );
    interaction_trace.finalize_last()
}
