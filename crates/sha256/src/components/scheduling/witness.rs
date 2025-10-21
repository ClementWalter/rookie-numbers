use std::simd::u32x16;

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
    combine,
    components::{
        combine_w,
        scheduling::columns::{RoundColumns, RoundInteractionColumns},
        W_SIZE,
    },
    consume_pair, emit_col,
    partitions::{pext_u32x16, Sigma0, Sigma1},
    relations::Relations,
    sha256::{small_sigma_0_u32x16, small_sigma_1_u32x16, CHUNK_SIZE, N_SCHEDULING_ROUNDS},
};

const N_COLUMNS: usize = W_SIZE + RoundColumns::SIZE * N_SCHEDULING_ROUNDS;
const N_INTERACTION_COLUMNS: usize = W_SIZE + RoundInteractionColumns::SIZE * N_SCHEDULING_ROUNDS;

#[inline]
fn generate_simd_sequence_bulk(start: usize, len: usize) -> Vec<u32x16> {
    assert!(len.is_multiple_of(16));
    let n = len / 16;
    let base = start as u32;

    (0..n)
        .map(|k| {
            let b = base + (k as u32) * 16;
            u32x16::from_array(core::array::from_fn(|i| (b + i as u32) & 0xffff))
        })
        .collect()
}

#[allow(clippy::type_complexity)]
pub fn gen_trace(
    log_size: u32,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    Vec<Vec<u32x16>>,
) {
    assert!(log_size >= LOG_N_LANES);
    let simd_size = 1 << (log_size - LOG_N_LANES);

    // Initialize vec for all groups of columns
    let mut evals: Vec<Vec<u32x16>> = (0..N_COLUMNS)
        .map(|_| Vec::with_capacity(simd_size))
        .collect::<Vec<_>>();
    let mut lookup_data: Vec<Vec<u32x16>> = (0..N_INTERACTION_COLUMNS)
        .map(|_| Vec::with_capacity(simd_size))
        .collect::<Vec<_>>();

    // Generate random inputs
    evals
        .iter_mut()
        .enumerate()
        .take(CHUNK_SIZE)
        .for_each(|(i, eval)| {
            *eval = generate_simd_sequence_bulk(i, 1 << log_size);
        });
    lookup_data
        .iter_mut()
        .enumerate()
        .take(CHUNK_SIZE)
        .for_each(|(i, eval)| {
            *eval = generate_simd_sequence_bulk(i, 1 << log_size);
        });

    for t in 16..(16 + N_SCHEDULING_ROUNDS) {
        let index = W_SIZE + (t - 16) * RoundColumns::SIZE;
        let interaction_index = W_SIZE + (t - 16) * RoundInteractionColumns::SIZE;

        for simd_row in 0..simd_size {
            // Load the W values
            let w_16_low = evals[2 * (t - 16)][simd_row];
            let w_16_high = evals[2 * (t - 16) + 1][simd_row];
            let w_15_low = evals[2 * (t - 15)][simd_row];
            let w_15_high = evals[2 * (t - 15) + 1][simd_row];
            let w_7_low = evals[2 * (t - 7)][simd_row];
            let w_7_high = evals[2 * (t - 7) + 1][simd_row];
            let w_2_low = evals[2 * (t - 2)][simd_row];
            let w_2_high = evals[2 * (t - 2) + 1][simd_row];

            // SIGMA0
            // Decomposition over I0
            let w_15_i0_low = w_15_low & u32x16::splat(Sigma0::I0_L);
            let w_15_i0_high = w_15_high & u32x16::splat(Sigma0::I0_H);
            let sigma_0 = small_sigma_0_u32x16(w_15_i0_low + (w_15_i0_high << 16));
            let sigma_0_o0_low = sigma_0 & u32x16::splat(Sigma0::O0_L);
            let sigma_0_o0_high = (sigma_0 >> 16) & u32x16::splat(Sigma0::O0_H);
            let sigma_0_o20 = sigma_0 & u32x16::splat(Sigma0::O2);
            let sigma_0_o20_pext = pext_u32x16(sigma_0_o20, Sigma0::O2);

            // Decomposition over I1
            let w_15_i1_low = w_15_low & u32x16::splat(Sigma0::I1_L);
            let w_15_i1_high = w_15_high & u32x16::splat(Sigma0::I1_H);
            let sigma_0 = small_sigma_0_u32x16(w_15_i1_low + (w_15_i1_high << 16));
            let sigma_0_o1_low = sigma_0 & u32x16::splat(Sigma0::O1_L);
            let sigma_0_o1_high = (sigma_0 >> 16) & u32x16::splat(Sigma0::O1_H);
            let sigma_0_o21 = sigma_0 & u32x16::splat(Sigma0::O2);
            let sigma_0_o21_pext = pext_u32x16(sigma_0_o21, Sigma0::O2);

            // XOR the two O2 values
            let sigma_0_o2 = sigma_0_o20 ^ sigma_0_o21;
            let sigma_0_o2_low = sigma_0_o2 & u32x16::splat(0xffff);
            let sigma_0_o2_high = sigma_0_o2 >> 16;

            // Compute sigma_0 output
            let sigma_0_low = sigma_0_o0_low + sigma_0_o1_low + sigma_0_o2_low;
            let sigma_0_high = sigma_0_o0_high + sigma_0_o1_high + sigma_0_o2_high;

            // SIGMA1
            // Decomposition over I0
            let w_2_i0_low = w_2_low & u32x16::splat(Sigma1::I0_L);
            let w_2_i0_high = w_2_high & u32x16::splat(Sigma1::I0_H);
            let sigma_1 = small_sigma_1_u32x16(w_2_i0_low + (w_2_i0_high << 16));
            let sigma_1_o0_low = sigma_1 & u32x16::splat(Sigma1::O0_L);
            let sigma_1_o0_high = (sigma_1 >> 16) & u32x16::splat(Sigma1::O0_H);
            let sigma_1_o20 = sigma_1 & u32x16::splat(Sigma1::O2);
            let sigma_1_o20_pext = pext_u32x16(sigma_1_o20, Sigma1::O2);

            // Decomposition over I1
            let w_2_i1_low = w_2_low & u32x16::splat(Sigma1::I1_L);
            let w_2_i1_high = w_2_high & u32x16::splat(Sigma1::I1_H);
            let sigma_1 = small_sigma_1_u32x16(w_2_i1_low + (w_2_i1_high << 16));
            let sigma_1_o1_low = sigma_1 & u32x16::splat(Sigma1::O1_L);
            let sigma_1_o1_high = (sigma_1 >> 16) & u32x16::splat(Sigma1::O1_H);
            let sigma_1_o21 = sigma_1 & u32x16::splat(Sigma1::O2);
            let sigma_1_o21_pext = pext_u32x16(sigma_1_o21, Sigma1::O2);

            // XOR the two O2 values
            let sigma_1_o2 = sigma_1_o20 ^ sigma_1_o21;
            let sigma_1_o2_low = sigma_1_o2 & u32x16::splat(0xffff);
            let sigma_1_o2_high = sigma_1_o2 >> 16;

            // Compute sigma_1 output
            let sigma_1_low = sigma_1_o0_low + sigma_1_o1_low + sigma_1_o2_low;
            let sigma_1_high = sigma_1_o0_high + sigma_1_o1_high + sigma_1_o2_high;

            // Compute the final output
            let round_low = w_16_low + sigma_0_low + w_7_low + sigma_1_low;
            let round_high = w_16_high + sigma_0_high + w_7_high + sigma_1_high;
            let carry_low = round_low >> 16;
            let carry_high = (round_high + carry_low) >> 16;
            let new_w_low = round_low - (carry_low << 16);
            let new_w_high = round_high + carry_low - (carry_high << 16);

            let trace_values: RoundColumns<u32x16> = RoundColumns {
                w_15_i0_low: &w_15_i0_low,
                w_15_i0_high: &w_15_i0_high,
                sigma_0_o0_low: &sigma_0_o0_low,
                sigma_0_o0_high: &sigma_0_o0_high,
                sigma_0_o20_pext: &sigma_0_o20_pext,
                sigma_0_o1_low: &sigma_0_o1_low,
                sigma_0_o1_high: &sigma_0_o1_high,
                sigma_0_o21_pext: &sigma_0_o21_pext,
                sigma_0_o2_low: &sigma_0_o2_low,
                sigma_0_o2_high: &sigma_0_o2_high,
                w_2_i0_low: &w_2_i0_low,
                w_2_i0_high: &w_2_i0_high,
                sigma_1_o0_low: &sigma_1_o0_low,
                sigma_1_o0_high: &sigma_1_o0_high,
                sigma_1_o20_pext: &sigma_1_o20_pext,
                sigma_1_o1_low: &sigma_1_o1_low,
                sigma_1_o1_high: &sigma_1_o1_high,
                sigma_1_o21_pext: &sigma_1_o21_pext,
                sigma_1_o2_low: &sigma_1_o2_low,
                sigma_1_o2_high: &sigma_1_o2_high,
                carry_low: &carry_low,
                carry_high: &carry_high,
            };
            for (i, value) in trace_values.iter().enumerate() {
                evals[index + i].push(*value);
            }

            let interaction_values: RoundInteractionColumns<u32x16> = RoundInteractionColumns {
                w_15_i0_low: &w_15_i0_low,
                w_15_i0_high: &w_15_i0_high,
                sigma_0_o0_low: &sigma_0_o0_low,
                sigma_0_o0_high: &sigma_0_o0_high,
                sigma_0_o20_pext: &sigma_0_o20_pext,
                w_15_i1_low: &w_15_i1_low,
                w_15_i1_high: &w_15_i1_high,
                sigma_0_o1_low: &sigma_0_o1_low,
                sigma_0_o1_high: &sigma_0_o1_high,
                sigma_0_o21_pext: &sigma_0_o21_pext,
                sigma_0_o2_low: &sigma_0_o2_low,
                sigma_0_o2_high: &sigma_0_o2_high,
                w_2_i0_low: &w_2_i0_low,
                w_2_i0_high: &w_2_i0_high,
                sigma_1_o0_low: &sigma_1_o0_low,
                sigma_1_o0_high: &sigma_1_o0_high,
                sigma_1_o20_pext: &sigma_1_o20_pext,
                w_2_i1_low: &w_2_i1_low,
                w_2_i1_high: &w_2_i1_high,
                sigma_1_o1_low: &sigma_1_o1_low,
                sigma_1_o1_high: &sigma_1_o1_high,
                sigma_1_o21_pext: &sigma_1_o21_pext,
                sigma_1_o2_low: &sigma_1_o2_low,
                sigma_1_o2_high: &sigma_1_o2_high,
                new_w_low: &new_w_low,
                new_w_high: &new_w_high,
                carry_low: &carry_low,
                carry_high: &carry_high,
            };
            for (i, value) in interaction_values.iter().enumerate() {
                lookup_data[interaction_index + i].push(*value);
            }

            evals[2 * t].push(new_w_low);
            evals[2 * t + 1].push(new_w_high);
            lookup_data[2 * t].push(new_w_low);
            lookup_data[2 * t + 1].push(new_w_high);
        }
    }

    let domain = CanonicCoset::new(log_size).circle_domain();
    let trace = evals
        .into_iter()
        .map(|values| {
            CircleEvaluation::new(
                domain,
                BaseColumn::from_simd(
                    values
                        .into_iter()
                        .map(|simd_chunk| unsafe { PackedM31::from_simd_unchecked(simd_chunk) })
                        .collect(),
                ),
            )
        })
        .collect();

    (trace, lookup_data)
}

pub fn gen_interaction_trace(
    lookup_data: &[Vec<u32x16>],
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    QM31,
) {
    let simd_size = lookup_data[0].len();
    let mut interaction_trace = LogupTraceGenerator::new(simd_size.ilog2() + LOG_N_LANES);

    for round in lookup_data[W_SIZE..].array_chunks::<{ RoundInteractionColumns::SIZE }>() {
        let RoundInteractionColumns {
            w_15_i0_low,
            w_15_i0_high,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
            w_15_i1_low,
            w_15_i1_high,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
            w_2_i0_low,
            w_2_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext,
            w_2_i1_low,
            w_2_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
            sigma_1_o2_low,
            sigma_1_o2_high,
            new_w_low,
            new_w_high,
            carry_low,
            carry_high,
        } = RoundInteractionColumns::from_slice(round);

        // SIGMA 0
        let sigma_0_i0 = combine!(
            relations.sigma_0.i0,
            w_15_i0_low,
            w_15_i0_high,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext
        );
        let sigma_0_i1 = combine!(
            relations.sigma_0.i1,
            w_15_i1_low,
            w_15_i1_high,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext
        );
        let sigma_0_o2 = combine!(
            relations.sigma_0.o2,
            sigma_0_o20_pext,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high
        );
        // SIGMA 1
        let sigma_1_i0 = combine!(
            relations.sigma_1.i0,
            w_2_i0_low,
            w_2_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext
        );
        let sigma_1_i1 = combine!(
            relations.sigma_1.i1,
            w_2_i1_low,
            w_2_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext
        );
        let sigma_1_o2 = combine!(
            relations.sigma_1.o2,
            sigma_1_o20_pext,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high
        );
        // ADD
        let carry_low = combine!(relations.range_check_add.add_4, new_w_low, carry_low,);
        let carry_high = combine!(relations.range_check_add.add_4, new_w_high, carry_high);

        consume_pair!(
            interaction_trace;
            sigma_0_i0, sigma_0_i1,


            sigma_0_o2, sigma_1_i0, sigma_1_i1, sigma_1_o2, carry_low,
            carry_high,
        );
    }

    // Emit W consumed by compression
    emit_col!(combine_w(relations, lookup_data), interaction_trace);

    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use stwo::prover::backend::Column;

    use super::*;
    use crate::sha256::{small_sigma_0, small_sigma_1};

    #[test]
    fn test_generate_simd_sequence_bulk() {
        let sequence = generate_simd_sequence_bulk(0, 1 << LOG_N_LANES);
        assert_eq!(sequence.len(), 1);
        assert_eq!(
            sequence[0],
            u32x16::from_array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        );
    }

    #[test]
    fn test_gen_trace_columns_count() {
        let (trace, _) = gen_trace(LOG_N_LANES + 3);
        assert_eq!(trace.len(), N_COLUMNS);
    }

    #[test]
    fn test_gen_trace_values() {
        let log_size = LOG_N_LANES;
        let size = 1 << log_size;
        let (mut trace, _) = gen_trace(log_size);
        let w = trace
            .drain(0..W_SIZE)
            .map(|eval| {
                eval.values
                    .to_cpu()
                    .into_iter()
                    .map(|x| x.0)
                    .collect::<Vec<u32>>()
            })
            .collect::<Vec<Vec<u32>>>();
        for row in 0..size {
            for t in 16..64 {
                let w_16_low = w[2 * (t - 16)][row];
                let w_16_high = w[2 * (t - 16) + 1][row];
                let w_16 = w_16_low + (w_16_high << 16);
                let w_15_low = w[2 * (t - 15)][row];
                let w_15_high = w[2 * (t - 15) + 1][row];
                let w_15 = w_15_low + (w_15_high << 16);
                let w_7_low = w[2 * (t - 7)][row];
                let w_7_high = w[2 * (t - 7) + 1][row];
                let w_7 = w_7_low + (w_7_high << 16);
                let w_2_low = w[2 * (t - 2)][row];
                let w_2_high = w[2 * (t - 2) + 1][row];
                let w_2 = w_2_low + (w_2_high << 16);
                let w_current_low = w[2 * t][row];
                let w_current_high = w[2 * t + 1][row];
                let w_current = w_current_low + (w_current_high << 16);
                assert_eq!(
                    w_current,
                    w_16.overflowing_add(small_sigma_0(w_15))
                        .0
                        .overflowing_add(w_7.overflowing_add(small_sigma_1(w_2)).0)
                        .0
                );
            }
        }
    }
}
