//! The compression component is responsible for proving the compression part of sha256.
//!
//! This is, 64 iterations of the main loop, doing Sigma0, Sigma1 and adding the
//! results with the buffer values.

use std::simd::u32x16;

use itertools::izip;
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
        compression::columns::{RoundColumns, RoundInteractionColumns},
        W_SIZE,
    },
    consume_col, consume_pair,
    partitions::{pext_u32x16, BigSigma0, BigSigma1},
    relations::Relations,
    sha256::{
        big_sigma1_u32x16, big_sigma_0_u32x16, ch_left_u32x16, ch_right_u32x16, maj_u32x16, H, K,
        N_COMPRESSION_ROUNDS,
    },
};

const N_COLUMNS: usize = W_SIZE + RoundColumns::SIZE * N_COMPRESSION_ROUNDS;
const N_INTERACTION_COLUMNS: usize = W_SIZE + RoundInteractionColumns::SIZE * N_COMPRESSION_ROUNDS;

#[allow(clippy::type_complexity)]
pub fn gen_trace(
    w: &ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    Vec<Vec<u32x16>>,
) {
    let simd_size = w[0].data.len();

    // Initialize vec for all groups of columns
    let mut evals: Vec<Vec<u32x16>> = (0..N_COLUMNS)
        .map(|_| Vec::with_capacity(simd_size))
        .collect::<Vec<_>>();
    let mut lookup_data: Vec<Vec<u32x16>> = (0..N_INTERACTION_COLUMNS)
        .map(|_| Vec::with_capacity(simd_size))
        .collect::<Vec<_>>();

    // Generate round constants
    let k: [u32x16; K.len() * 2] = K
        .iter()
        .flat_map(|k| [u32x16::splat(k & 0xffff), u32x16::splat(k >> 16)])
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    // Get initial hash value
    let mut hash_buffer: [Vec<u32x16>; H.len() * 2] = H
        .iter()
        .flat_map(|h| {
            [
                vec![u32x16::splat(h & 0xffff); simd_size],
                vec![u32x16::splat(h >> 16); simd_size],
            ]
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    // Fill initial trace and lookup data
    evals
        .iter_mut()
        .enumerate()
        .take(W_SIZE)
        .for_each(|(i, eval)| {
            *eval = w[i]
                .data
                .clone()
                .into_iter()
                .map(|x| x.into_simd())
                .collect();
        });
    lookup_data
        .iter_mut()
        .enumerate()
        .take(W_SIZE)
        .for_each(|(i, eval)| {
            *eval = w[i]
                .data
                .clone()
                .into_iter()
                .map(|x| x.into_simd())
                .collect();
        });

    for round in 0..N_COMPRESSION_ROUNDS {
        let index = W_SIZE + round * RoundColumns::SIZE;
        let interaction_index = W_SIZE + round * RoundInteractionColumns::SIZE;

        let a_low = &hash_buffer[0].clone();
        let a_high = &hash_buffer[1].clone();
        let b_low = &hash_buffer[2].clone();
        let b_high = &hash_buffer[3].clone();
        let c_low = &hash_buffer[4].clone();
        let c_high = &hash_buffer[5].clone();
        let d_low = &hash_buffer[6].clone();
        let d_high = &hash_buffer[7].clone();
        let e_low = &hash_buffer[8].clone();
        let e_high = &hash_buffer[9].clone();
        let f_low = &hash_buffer[10].clone();
        let f_high = &hash_buffer[11].clone();
        let g_low = &hash_buffer[12].clone();
        let g_high = &hash_buffer[13].clone();
        let h_low = &hash_buffer[14].clone();
        let h_high = &hash_buffer[15].clone();

        // Load K value
        let k_low = k[2 * round];
        let k_high = k[2 * round + 1];

        for simd_row in 0..simd_size {
            // Load W value
            let w_low = evals[2 * round][simd_row];
            let w_high = evals[2 * round + 1][simd_row];

            // BIG_SIGMA1
            // Decomposition over I0
            let e_i0_low = e_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let e_i0_high = e_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let sigma_1 = big_sigma1_u32x16(e_i0_low + (e_i0_high << 16));
            let sigma_1_o0_low = sigma_1 & u32x16::splat(BigSigma1::O0_L);
            let sigma_1_o0_high = (sigma_1 >> 16) & u32x16::splat(BigSigma1::O0_H);
            let sigma_1_o20 = sigma_1 & u32x16::splat(BigSigma1::O2);
            let sigma_1_o20_pext = pext_u32x16(sigma_1_o20, BigSigma1::O2);

            // Decomposition over I1
            let e_i1_low = e_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let e_i1_high = e_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let sigma_1 = big_sigma1_u32x16(e_i1_low + (e_i1_high << 16));
            let sigma_1_o1_low = sigma_1 & u32x16::splat(BigSigma1::O1_L);
            let sigma_1_o1_high = (sigma_1 >> 16) & u32x16::splat(BigSigma1::O1_H);
            let sigma_1_o21 = sigma_1 & u32x16::splat(BigSigma1::O2);
            let sigma_1_o21_pext = pext_u32x16(sigma_1_o21, BigSigma1::O2);

            // XOR the two O2 values
            let big_sigma_1_o2 = sigma_1_o20 ^ sigma_1_o21;
            let sigma_1_o2_low = big_sigma_1_o2 & u32x16::splat(0xffff);
            let sigma_1_o2_high = big_sigma_1_o2 >> 16;

            // Output
            let sigma1_low = sigma_1_o0_low + sigma_1_o1_low + sigma_1_o2_low;
            let sigma1_high = sigma_1_o0_high + sigma_1_o1_high + sigma_1_o2_high;

            // CH
            // left side
            let f_i0_low = f_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let f_i0_high = f_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let f_i1_low = f_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let f_i1_high = f_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let ch_left_i0_low = ch_left_u32x16(e_i0_low, f_i0_low);
            let ch_left_i0_high = ch_left_u32x16(e_i0_high, f_i0_high);
            let ch_left_i1_low = ch_left_u32x16(e_i1_low, f_i1_low);
            let ch_left_i1_high = ch_left_u32x16(e_i1_high, f_i1_high);

            // right side
            let g_i0_low = g_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let g_i0_high = g_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let g_i1_low = g_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let g_i1_high = g_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let ch_right_i0_low = ch_right_u32x16(e_i0_low, g_i0_low);
            let ch_right_i0_high = ch_right_u32x16(e_i0_high, g_i0_high);
            let ch_right_i1_low = ch_right_u32x16(e_i1_low, g_i1_low);
            let ch_right_i1_high = ch_right_u32x16(e_i1_high, g_i1_high);

            // Output
            let ch_low = ch_left_i0_low + ch_left_i1_low + ch_right_i0_low + ch_right_i1_low;
            let ch_high = ch_left_i0_high + ch_left_i1_high + ch_right_i0_high + ch_right_i1_high;

            // BIG_SIGMA0
            // Decomposition over I0
            let a_i0_low = a_low[simd_row] & u32x16::splat(BigSigma0::I0_L);
            let a_i0_high_0 = a_high[simd_row] & u32x16::splat(BigSigma0::I0_H0);
            let a_i0_high_1 = (a_high[simd_row] >> 8) & u32x16::splat(BigSigma0::I0_H1);

            let sigma_0 = big_sigma_0_u32x16(a_i0_low + (a_i0_high_0 << 16) + (a_i0_high_1 << 24));
            let sigma_0_o0_low = sigma_0 & u32x16::splat(BigSigma0::O0_L);
            let sigma_0_o0_high = (sigma_0 >> 16) & u32x16::splat(BigSigma0::O0_H);
            let sigma_0_o20 = sigma_0 & u32x16::splat(BigSigma0::O2);
            let sigma_0_o20_pext = pext_u32x16(sigma_0_o20, BigSigma0::O2);

            // Decomposition over I1
            let a_i1_low_0 = a_low[simd_row] & u32x16::splat(BigSigma0::I1_L0);
            let a_i1_low_1 = (a_low[simd_row] >> 8) & u32x16::splat(BigSigma0::I1_L1);
            let a_i1_high = a_high[simd_row] & u32x16::splat(BigSigma0::I1_H);

            let sigma_0 = big_sigma_0_u32x16(a_i1_low_0 + (a_i1_low_1 << 8) + (a_i1_high << 16));
            let sigma_0_o1_low = sigma_0 & u32x16::splat(BigSigma0::O1_L);
            let sigma_0_o1_high = (sigma_0 >> 16) & u32x16::splat(BigSigma0::O1_H);
            let sigma_0_o21 = sigma_0 & u32x16::splat(BigSigma0::O2);
            let sigma_0_o21_pext = pext_u32x16(sigma_0_o21, BigSigma0::O2);

            // XOR the two O2 values
            let big_sigma_0_o2 = sigma_0_o20 ^ sigma_0_o21;
            let sigma_0_o2_low = big_sigma_0_o2 & u32x16::splat(0xffff);
            let sigma_0_o2_high = big_sigma_0_o2 >> 16;

            // Output
            let sigma_0_low = sigma_0_o0_low + sigma_0_o1_low + sigma_0_o2_low;
            let sigma_0_high = sigma_0_o0_high + sigma_0_o1_high + sigma_0_o2_high;

            // MAJ
            let b_i0_low = b_low[simd_row] & u32x16::splat(BigSigma0::I0_L);
            let b_i0_high_0 = b_high[simd_row] & u32x16::splat(BigSigma0::I0_H0);
            let b_i0_high_1 = (b_high[simd_row] >> 8) & u32x16::splat(BigSigma0::I0_H1);
            let b_i1_low_0 = b_low[simd_row] & u32x16::splat(BigSigma0::I1_L0);
            let b_i1_low_1 = (b_low[simd_row] >> 8) & u32x16::splat(BigSigma0::I1_L1);
            let b_i1_high = b_high[simd_row] & u32x16::splat(BigSigma0::I1_H);
            let c_i0_low = c_low[simd_row] & u32x16::splat(BigSigma0::I0_L);
            let c_i0_high_0 = c_high[simd_row] & u32x16::splat(BigSigma0::I0_H0);
            let c_i0_high_1 = (c_high[simd_row] >> 8) & u32x16::splat(BigSigma0::I0_H1);
            let c_i1_low_0 = c_low[simd_row] & u32x16::splat(BigSigma0::I1_L0);
            let c_i1_low_1 = (c_low[simd_row] >> 8) & u32x16::splat(BigSigma0::I1_L1);
            let c_i1_high = c_high[simd_row] & u32x16::splat(BigSigma0::I1_H);
            let maj_i0_low = maj_u32x16(a_i0_low, b_i0_low, c_i0_low);
            let maj_i0_high_0 = maj_u32x16(a_i0_high_0, b_i0_high_0, c_i0_high_0);
            let maj_i0_high_1 = maj_u32x16(a_i0_high_1, b_i0_high_1, c_i0_high_1);
            let maj_i1_low_0 = maj_u32x16(a_i1_low_0, b_i1_low_0, c_i1_low_0);
            let maj_i1_low_1 = maj_u32x16(a_i1_low_1, b_i1_low_1, c_i1_low_1);
            let maj_i1_high = maj_u32x16(a_i1_high, b_i1_high, c_i1_high);

            // Output
            let maj_low = maj_i0_low + maj_i1_low_0 + (maj_i1_low_1 << 8);
            let maj_high = maj_i0_high_0 + (maj_i0_high_1 << 8) + maj_i1_high;

            // TEMP
            let temp1_low = h_low[simd_row] + sigma1_low + ch_low + k_low + w_low;
            let temp1_high = h_high[simd_row] + sigma1_high + ch_high + k_high + w_high;
            let temp2_low = sigma_0_low + maj_low;
            let temp2_high = sigma_0_high + maj_high;

            let e_carry_low = (temp1_low + d_low[simd_row]) >> 16;
            let e_carry_high = (temp1_high + d_high[simd_row] + e_carry_low) >> 16;
            let new_e_low = temp1_low + d_low[simd_row] - (e_carry_low << 16);
            let new_e_high = temp1_high + d_high[simd_row] + e_carry_low - (e_carry_high << 16);
            let a_carry_low = (temp1_low + temp2_low) >> 16;
            let a_carry_high = (temp1_high + temp2_high + a_carry_low) >> 16;
            let new_a_low = temp1_low + temp2_low - (a_carry_low << 16);
            let new_a_high = temp1_high + temp2_high + a_carry_low - (a_carry_high << 16);

            let trace_values: RoundColumns<u32x16> = RoundColumns {
                e_i0_low: &e_i0_low,
                e_i0_high: &e_i0_high,
                sigma_1_o0_low: &sigma_1_o0_low,
                sigma_1_o0_high: &sigma_1_o0_high,
                sigma_1_o20_pext: &sigma_1_o20_pext,
                sigma_1_o1_low: &sigma_1_o1_low,
                sigma_1_o1_high: &sigma_1_o1_high,
                sigma_1_o21_pext: &sigma_1_o21_pext,
                sigma_1_o2_low: &sigma_1_o2_low,
                sigma_1_o2_high: &sigma_1_o2_high,
                f_i0_low: &f_i0_low,
                f_i0_high: &f_i0_high,
                ch_left_i0_low: &ch_left_i0_low,
                ch_left_i0_high: &ch_left_i0_high,
                ch_left_i1_low: &ch_left_i1_low,
                ch_left_i1_high: &ch_left_i1_high,
                g_i0_low: &g_i0_low,
                g_i0_high: &g_i0_high,
                ch_right_i0_low: &ch_right_i0_low,
                ch_right_i0_high: &ch_right_i0_high,
                ch_right_i1_low: &ch_right_i1_low,
                ch_right_i1_high: &ch_right_i1_high,
                a_i0_high_0: &a_i0_high_0,
                a_i0_high_1: &a_i0_high_1,
                a_i1_low_0: &a_i1_low_0,
                a_i1_low_1: &a_i1_low_1,
                sigma_0_o0_low: &sigma_0_o0_low,
                sigma_0_o0_high: &sigma_0_o0_high,
                sigma_0_o20_pext: &sigma_0_o20_pext,
                sigma_0_o1_low: &sigma_0_o1_low,
                sigma_0_o1_high: &sigma_0_o1_high,
                sigma_0_o21_pext: &sigma_0_o21_pext,
                sigma_0_o2_low: &sigma_0_o2_low,
                sigma_0_o2_high: &sigma_0_o2_high,
                b_i0_high_0: &b_i0_high_0,
                b_i0_high_1: &b_i0_high_1,
                b_i1_low_0: &b_i1_low_0,
                b_i1_low_1: &b_i1_low_1,
                c_i0_high_0: &c_i0_high_0,
                c_i0_high_1: &c_i0_high_1,
                c_i1_low_0: &c_i1_low_0,
                c_i1_low_1: &c_i1_low_1,
                maj_i0_low: &maj_i0_low,
                maj_i0_high_0: &maj_i0_high_0,
                maj_i0_high_1: &maj_i0_high_1,
                maj_i1_low_0: &maj_i1_low_0,
                maj_i1_low_1: &maj_i1_low_1,
                maj_i1_high: &maj_i1_high,
                e_carry_low: &e_carry_low,
                e_carry_high: &e_carry_high,
                a_carry_low: &a_carry_low,
                a_carry_high: &a_carry_high,
            };
            for (i, value) in trace_values.iter().enumerate() {
                evals[index + i].push(*value);
            }

            let interaction_values: RoundInteractionColumns<u32x16> = RoundInteractionColumns {
                e_i0_low: &e_i0_low,
                e_i0_high: &e_i0_high,
                sigma_1_o0_low: &sigma_1_o0_low,
                sigma_1_o0_high: &sigma_1_o0_high,
                sigma_1_o20_pext: &sigma_1_o20_pext,
                e_i1_low: &e_i1_low,
                e_i1_high: &e_i1_high,
                sigma_1_o1_low: &sigma_1_o1_low,
                sigma_1_o1_high: &sigma_1_o1_high,
                sigma_1_o21_pext: &sigma_1_o21_pext,
                sigma_1_o2_low: &sigma_1_o2_low,
                sigma_1_o2_high: &sigma_1_o2_high,
                f_i0_low: &f_i0_low,
                f_i0_high: &f_i0_high,
                f_i1_low: &f_i1_low,
                f_i1_high: &f_i1_high,
                ch_left_i0_low: &ch_left_i0_low,
                ch_left_i0_high: &ch_left_i0_high,
                ch_left_i1_low: &ch_left_i1_low,
                ch_left_i1_high: &ch_left_i1_high,
                g_i0_low: &g_i0_low,
                g_i0_high: &g_i0_high,
                g_i1_low: &g_i1_low,
                g_i1_high: &g_i1_high,
                ch_right_i0_low: &ch_right_i0_low,
                ch_right_i0_high: &ch_right_i0_high,
                ch_right_i1_low: &ch_right_i1_low,
                ch_right_i1_high: &ch_right_i1_high,
                a_i0_low: &a_i0_low,
                a_i0_high_0: &a_i0_high_0,
                a_i0_high_1: &a_i0_high_1,
                sigma_0_o0_low: &sigma_0_o0_low,
                sigma_0_o0_high: &sigma_0_o0_high,
                sigma_0_o20_pext: &sigma_0_o20_pext,
                a_i1_low_0: &a_i1_low_0,
                a_i1_low_1: &a_i1_low_1,
                a_i1_high: &a_i1_high,
                sigma_0_o1_low: &sigma_0_o1_low,
                sigma_0_o1_high: &sigma_0_o1_high,
                sigma_0_o21_pext: &sigma_0_o21_pext,
                sigma_0_o2_low: &sigma_0_o2_low,
                sigma_0_o2_high: &sigma_0_o2_high,
                b_i0_low: &b_i0_low,
                b_i0_high_0: &b_i0_high_0,
                b_i0_high_1: &b_i0_high_1,
                b_i1_low_0: &b_i1_low_0,
                b_i1_low_1: &b_i1_low_1,
                b_i1_high: &b_i1_high,
                c_i0_low: &c_i0_low,
                c_i0_high_0: &c_i0_high_0,
                c_i0_high_1: &c_i0_high_1,
                c_i1_low_0: &c_i1_low_0,
                c_i1_low_1: &c_i1_low_1,
                c_i1_high: &c_i1_high,
                maj_i0_low: &maj_i0_low,
                maj_i0_high_0: &maj_i0_high_0,
                maj_i0_high_1: &maj_i0_high_1,
                maj_i1_low_0: &maj_i1_low_0,
                maj_i1_low_1: &maj_i1_low_1,
                maj_i1_high: &maj_i1_high,
                e_carry_low: &e_carry_low,
                e_carry_high: &e_carry_high,
                new_e_low: &new_e_low,
                new_e_high: &new_e_high,
                a_carry_low: &a_carry_low,
                a_carry_high: &a_carry_high,
                new_a_low: &new_a_low,
                new_a_high: &new_a_high,
            };
            for (i, value) in interaction_values.iter().enumerate() {
                lookup_data[interaction_index + i].push(*value);
            }
        }

        update_hash_buffer(&mut hash_buffer, &evals, round);
    }

    let domain = CanonicCoset::new(simd_size.ilog2() + LOG_N_LANES).circle_domain();
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

/// Update the hash buffer with the values from the trace
fn update_hash_buffer(hash_buffer: &mut [Vec<u32x16>], evals: &[Vec<u32x16>], round: usize) {
    let d_low = &hash_buffer[6];
    let d_high = &hash_buffer[7];
    let h_low = &hash_buffer[14];
    let h_high = &hash_buffer[15];

    let w_low = evals[2 * round].clone();
    let w_high = evals[2 * round + 1].clone();

    let k_low = u32x16::splat(K[round] & 0xffff);
    let k_high = u32x16::splat(K[round] >> 16);

    let index = W_SIZE + RoundColumns::SIZE * round;
    let RoundColumns {
        e_i0_low: _,
        e_i0_high: _,
        sigma_1_o0_low,
        sigma_1_o0_high,
        sigma_1_o20_pext: _,
        sigma_1_o1_low,
        sigma_1_o1_high,
        sigma_1_o21_pext: _,
        sigma_1_o2_low,
        sigma_1_o2_high,
        f_i0_low: _,
        f_i0_high: _,
        ch_left_i0_low,
        ch_left_i0_high,
        ch_left_i1_low,
        ch_left_i1_high,
        g_i0_low: _,
        g_i0_high: _,
        ch_right_i0_low,
        ch_right_i0_high,
        ch_right_i1_low,
        ch_right_i1_high,
        a_i0_high_0: _,
        a_i0_high_1: _,
        a_i1_low_0: _,
        a_i1_low_1: _,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext: _,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext: _,
        sigma_0_o2_low,
        sigma_0_o2_high,
        b_i0_high_0: _,
        b_i0_high_1: _,
        b_i1_low_0: _,
        b_i1_low_1: _,
        c_i0_high_0: _,
        c_i0_high_1: _,
        c_i1_low_0: _,
        c_i1_low_1: _,
        maj_i0_low,
        maj_i0_high_0,
        maj_i0_high_1,
        maj_i1_low_0,
        maj_i1_low_1,
        maj_i1_high,
        e_carry_low,
        e_carry_high,
        a_carry_low,
        a_carry_high,
    } = RoundColumns::from_slice(&evals[index..(index + RoundColumns::SIZE)]);

    let sigma1_low: Vec<u32x16> = izip!(sigma_1_o0_low, sigma_1_o1_low, sigma_1_o2_low)
        .map(|(a, b, c)| a + b + c)
        .collect();
    let sigma1_high: Vec<u32x16> = izip!(sigma_1_o0_high, sigma_1_o1_high, sigma_1_o2_high)
        .map(|(a, b, c)| a + b + c)
        .collect();

    let ch_low: Vec<u32x16> = izip!(
        ch_left_i0_low,
        ch_left_i1_low,
        ch_right_i0_low,
        ch_right_i1_low
    )
    .map(|(a, b, c, d)| a + b + c + d)
    .collect();
    let ch_high: Vec<u32x16> = izip!(
        ch_left_i0_high,
        ch_left_i1_high,
        ch_right_i0_high,
        ch_right_i1_high
    )
    .map(|(a, b, c, d)| a + b + c + d)
    .collect();

    let sigma_0_high: Vec<u32x16> = izip!(sigma_0_o0_high, sigma_0_o1_high, sigma_0_o2_high)
        .map(|(a, b, c)| a + b + c)
        .collect();
    let sigma_0_low: Vec<u32x16> = izip!(sigma_0_o0_low, sigma_0_o1_low, sigma_0_o2_low)
        .map(|(a, b, c)| a + b + c)
        .collect();

    let maj_high: Vec<u32x16> = izip!(maj_i0_high_0, maj_i0_high_1, maj_i1_high)
        .map(|(a, b, c)| a + (b << 8) + c)
        .collect();
    let maj_low: Vec<u32x16> = izip!(maj_i0_low, maj_i1_low_0, maj_i1_low_1)
        .map(|(a, b, c)| a + b + (c << 8))
        .collect();

    let temp1_high: Vec<u32x16> = izip!(h_high, sigma1_high, ch_high, w_high)
        .map(|(a, b, c, d)| a + b + c + d + k_high)
        .collect();
    let temp1_low: Vec<u32x16> = izip!(h_low, sigma1_low, ch_low, w_low)
        .map(|(a, b, c, d)| a + b + c + d + k_low)
        .collect();

    let temp2_high: Vec<u32x16> = izip!(sigma_0_high, maj_high).map(|(a, b)| a + b).collect();
    let temp2_low: Vec<u32x16> = izip!(sigma_0_low, maj_low).map(|(a, b)| a + b).collect();

    let e_low: Vec<u32x16> = izip!(d_low.clone(), temp1_low.clone(), e_carry_low)
        .map(|(a, b, carry_low)| a + b - (carry_low << u32x16::splat(16)))
        .collect();
    let e_high: Vec<u32x16> = izip!(d_high, temp1_high.clone(), e_carry_low, e_carry_high)
        .map(|(a, b, carry_low, carry_high)| a + b + carry_low - (carry_high << u32x16::splat(16)))
        .collect();

    let a_low: Vec<u32x16> = izip!(temp1_low, temp2_low, a_carry_low.clone())
        .map(|(a, b, carry_low)| a + b - (carry_low << u32x16::splat(16)))
        .collect();
    let a_high = izip!(temp1_high, temp2_high, a_carry_low, a_carry_high)
        .map(|(a, b, carry_low, carry_high)| a + b + carry_low - (carry_high << u32x16::splat(16)))
        .collect();

    hash_buffer[15] = hash_buffer[13].clone(); // h_high = g_high
    hash_buffer[14] = hash_buffer[12].clone(); // h_low = g_low
    hash_buffer[13] = hash_buffer[11].clone(); // g_high = f_high
    hash_buffer[12] = hash_buffer[10].clone(); // g_low = f_low
    hash_buffer[11] = hash_buffer[9].clone(); // f_high = e_high
    hash_buffer[10] = hash_buffer[8].clone(); // f_low = e_low
    hash_buffer[9] = e_high; // e_high = d_high + temp1_high
    hash_buffer[8] = e_low; // e_low = d_low + temp1_low
    hash_buffer[7] = hash_buffer[5].clone(); // d_high = c_high
    hash_buffer[6] = hash_buffer[4].clone(); // d_low = c_low
    hash_buffer[5] = hash_buffer[3].clone(); // c_high = b_high
    hash_buffer[4] = hash_buffer[2].clone(); // c_low = b_low
    hash_buffer[3] = hash_buffer[1].clone(); // b_high = a_high
    hash_buffer[2] = hash_buffer[0].clone(); // b_low = a_low
    hash_buffer[1] = a_high; // a_high = temp1_high + temp2_high
    hash_buffer[0] = a_low; // a_low = temp1_low + temp2_low
}

#[allow(clippy::cognitive_complexity)]
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
            e_i0_low,
            e_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext,
            e_i1_low,
            e_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high,
            f_i0_low,
            f_i0_high,
            f_i1_low,
            f_i1_high,
            ch_left_i0_low,
            ch_left_i0_high,
            ch_left_i1_low,
            ch_left_i1_high,
            g_i0_low,
            g_i0_high,
            g_i1_low,
            g_i1_high,
            ch_right_i0_low,
            ch_right_i0_high,
            ch_right_i1_low,
            ch_right_i1_high,
            a_i0_low,
            a_i0_high_0,
            a_i0_high_1,
            a_i1_low_0,
            a_i1_low_1,
            a_i1_high,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
            b_i0_low,
            b_i0_high_0,
            b_i0_high_1,
            b_i1_low_0,
            b_i1_low_1,
            b_i1_high,
            c_i0_low,
            c_i0_high_0,
            c_i0_high_1,
            c_i1_low_0,
            c_i1_low_1,
            c_i1_high,
            maj_i0_low,
            maj_i0_high_0,
            maj_i0_high_1,
            maj_i1_low_0,
            maj_i1_low_1,
            maj_i1_high,
            e_carry_low,
            e_carry_high,
            a_carry_low,
            a_carry_high,
            new_e_low,
            new_e_high,
            new_a_low,
            new_a_high,
        } = RoundInteractionColumns::from_slice(round);

        // BIG_SIGMA1
        let big_sigma_1_i0 = combine!(
            relations.big_sigma_1.i0,
            e_i0_low,
            e_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext
        );
        let big_sigma_1_i1 = combine!(
            relations.big_sigma_1.i1,
            e_i1_low,
            e_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext
        );
        let big_sigma_1_o2 = combine!(
            relations.big_sigma_1.o2,
            sigma_1_o20_pext,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high
        );
        // CH_LEFT
        let ch_left_i0_low = combine!(relations.ch_left.i0_low, e_i0_low, f_i0_low, ch_left_i0_low);
        let ch_left_i0_high = combine!(
            relations.ch_left.i0_high,
            e_i0_high,
            f_i0_high,
            ch_left_i0_high,
        );
        let ch_left_i1_low =
            combine!(relations.ch_left.i1_low, e_i1_low, f_i1_low, ch_left_i1_low,);
        let ch_left_i1_high = combine!(
            relations.ch_left.i1_high,
            e_i1_high,
            f_i1_high,
            ch_left_i1_high,
        );
        // CH_RIGHT
        let ch_right_i0_low = combine!(
            relations.ch_right.i0_low,
            e_i0_low,
            g_i0_low,
            ch_right_i0_low,
        );
        let ch_right_i0_high = combine!(
            relations.ch_right.i0_high,
            e_i0_high,
            g_i0_high,
            ch_right_i0_high,
        );
        let ch_right_i1_low = combine!(
            relations.ch_right.i1_low,
            e_i1_low,
            g_i1_low,
            ch_right_i1_low,
        );
        let ch_right_i1_high = combine!(
            relations.ch_right.i1_high,
            e_i1_high,
            g_i1_high,
            ch_right_i1_high,
        );
        // BIG SIGMA0
        let big_sigma_0_i0 = combine!(
            relations.big_sigma_0.i0,
            a_i0_low,
            a_i0_high_0,
            a_i0_high_1,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
        );
        let big_sigma_0_i1 = combine!(
            relations.big_sigma_0.i1,
            a_i1_low_0,
            a_i1_low_1,
            a_i1_high,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
        );
        let big_sigma_0_o2 = combine!(
            relations.big_sigma_0.o2,
            sigma_0_o20_pext,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
        );
        // MAJ
        let maj_i0_low = combine!(
            relations.maj.i0_low,
            a_i0_low,
            b_i0_low,
            c_i0_low,
            maj_i0_low,
        );
        let maj_i0_high_0 = combine!(
            relations.maj.i0_high_0,
            a_i0_high_0,
            b_i0_high_0,
            c_i0_high_0,
            maj_i0_high_0,
        );
        let maj_i0_high_1 = combine!(
            relations.maj.i0_high_1,
            a_i0_high_1,
            b_i0_high_1,
            c_i0_high_1,
            maj_i0_high_1,
        );
        let maj_i1_low_0 = combine!(
            relations.maj.i1_low_0,
            a_i1_low_0,
            b_i1_low_0,
            c_i1_low_0,
            maj_i1_low_0,
        );
        let maj_i1_low_1 = combine!(
            relations.maj.i1_low_1,
            a_i1_low_1,
            b_i1_low_1,
            c_i1_low_1,
            maj_i1_low_1,
        );
        let maj_i1_high = combine!(
            relations.maj.i1_high,
            a_i1_high,
            b_i1_high,
            c_i1_high,
            maj_i1_high,
        );
        // ADD
        let e_carry_low = combine!(relations.range_check_add.add_7, new_e_low, e_carry_low,);
        let e_carry_high = combine!(relations.range_check_add.add_7, new_e_high, e_carry_high,);
        let a_carry_low = combine!(relations.range_check_add.add_8, new_a_low, a_carry_low,);
        let a_carry_high = combine!(relations.range_check_add.add_8, new_a_high, a_carry_high,);

        let secure_columns = [
            big_sigma_1_i0,
            big_sigma_1_i1,
            big_sigma_1_o2,
            ch_left_i0_low,
            ch_left_i0_high,
            ch_left_i1_low,
            ch_left_i1_high,
            ch_right_i0_low,
            ch_right_i0_high,
            ch_right_i1_low,
            ch_right_i1_high,
            big_sigma_0_i0,
            big_sigma_0_i1,
            big_sigma_0_o2,
            maj_i0_low,
            maj_i0_high_0,
            maj_i0_high_1,
            maj_i1_low_0,
            maj_i1_low_1,
            maj_i1_high,
            e_carry_low,
            e_carry_high,
            a_carry_low,
            a_carry_high,
        ];
        for i in 0..(secure_columns.len() / 2) {
            consume_pair!(
                secure_columns[2 * i],
                secure_columns[2 * i + 1],
                interaction_trace
            );
        }
    }

    // Consume W emitted by scheduling
    consume_col!(combine_w(relations, lookup_data), interaction_trace);

    interaction_trace.finalize_last()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        components::scheduling::witness::gen_trace as gen_schedule,
        sha256::{process_chunk_u32x16, CHUNK_SIZE},
    };

    #[test]
    fn test_gen_trace_columns_count() {
        let (schedule, _) = gen_schedule(LOG_N_LANES);
        let (trace, _) = gen_trace(&schedule);
        assert_eq!(trace.len(), N_COLUMNS);
    }

    #[test]
    fn test_gen_trace_values() {
        let log_size = LOG_N_LANES;
        let (schedule, _) = gen_schedule(log_size);
        let (trace, _) = gen_trace(&schedule);
        let chunk = trace[0..CHUNK_SIZE]
            .iter()
            .map(|eval| {
                eval.data
                    .clone()
                    .into_iter()
                    .map(|x| x.into_simd())
                    .collect::<Vec<u32x16>>()
            })
            .collect::<Vec<_>>();
        let chunk = std::array::from_fn(|i| chunk[2 * i][0] + (chunk[2 * i + 1][0] << 16));
        let h = std::array::from_fn(|i| u32x16::splat(H[i]));
        let expected = process_chunk_u32x16(chunk, h);

        let mut hash_buffer: [Vec<u32x16>; H.len() * 2] = H
            .iter()
            .flat_map(|h| {
                [
                    vec![u32x16::splat(h & 0xffff); 1],
                    vec![u32x16::splat(h >> 16); 1],
                ]
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        let evals: [Vec<u32x16>; N_COLUMNS] = std::array::from_fn(|i| {
            trace[i]
                .data
                .clone()
                .into_iter()
                .map(|x| x.into_simd())
                .collect()
        });
        for round in 0..64 {
            update_hash_buffer(&mut hash_buffer, &evals, round);
        }

        let result: [u32x16; 8] = std::array::from_fn(|i| {
            hash_buffer[2 * i][0] + (hash_buffer[2 * i + 1][0] << 16) + u32x16::splat(H[i])
        });

        assert_eq!(result, expected);
    }
}
