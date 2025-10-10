//! The compression component is responsible for proving the compression part of sha256.
//!
//! This is, 64 iterations of the main loop, doing Sigma0, Sigma1 and adding the
//! results with the buffer values.

use crate::partitions::{pext_u32x16, BigSigma0, BigSigma1};
use crate::sha256::{
    big_sigma0_u32x16, big_sigma1_u32x16, ch_left_u32x16, ch_right_u32x16, maj_u32x16, H, K,
};
use itertools::izip;
use stwo_prover::core::backend::simd::column::BaseColumn;
use stwo_prover::core::backend::simd::m31::{PackedM31, LOG_N_LANES};
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::fields::m31::BaseField;
use stwo_prover::core::poly::circle::CanonicCoset;
use stwo_prover::core::poly::circle::CircleEvaluation;
use stwo_prover::core::poly::BitReversedOrder;
use stwo_prover::core::ColumnVec;
use tracing::span;
use tracing::Level;

use crate::components::scheduling::W_COLUMNS;
use std::simd::u32x16;

const N_ROUNDS: usize = 64;

const BIG_SIGMA1_COLUMNS: usize = 10;
const CH_COLUMNS: usize = 12;
const BIG_SIGMA0_COLUMNS: usize = 12;
const MAJ_COLUMNS: usize = 14;
const CARRIES_COLUMNS: usize = 4;

const COL_PER_ROUND: usize =
    BIG_SIGMA1_COLUMNS + CH_COLUMNS + BIG_SIGMA0_COLUMNS + MAJ_COLUMNS + CARRIES_COLUMNS;
const INTERACTION_COL_PER_ROUND: usize = COL_PER_ROUND + 2 * 6;
const N_COLUMNS: usize = W_COLUMNS + COL_PER_ROUND * N_ROUNDS;
const N_INTERACTION_COLUMNS: usize = W_COLUMNS + INTERACTION_COL_PER_ROUND * N_ROUNDS;

#[allow(dead_code)]
#[repr(usize)]
enum ColumnIndex {
    e_i0_low,
    e_i0_high,
    sigma1_o0_low,
    sigma1_o0_high,
    sigma1_o20_pext,
    sigma1_o1_low,
    sigma1_o1_high,
    sigma1_o21_pext,
    sigma1_o2_low,
    sigma1_o2_high,
    f_i0_low,
    f_i0_high,
    ch_left_i0_l,
    ch_left_i0_h,
    ch_left_i1_l,
    ch_left_i1_h,
    g_i0_low,
    g_i0_high,
    ch_right_i0_l,
    ch_right_i0_h,
    ch_right_i1_l,
    ch_right_i1_h,
    a_i0_high_0,
    a_i0_high_1,
    a_i1_low_0,
    a_i1_low_1,
    sigma0_o0_low,
    sigma0_o0_high,
    sigma0_o20_pext,
    sigma0_o1_low,
    sigma0_o1_high,
    sigma0_o21_pext,
    sigma0_o2_low,
    sigma0_o2_high,
    b_i0_high_0,
    b_i0_high_1,
    b_i1_low_0,
    b_i1_low_1,
    c_i0_high_0,
    c_i0_high_1,
    c_i1_low_0,
    c_i1_low_1,
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
}

#[inline(always)]
const fn trace_index(round: usize, column: ColumnIndex) -> usize {
    W_COLUMNS + round * COL_PER_ROUND + column as usize
}

#[allow(clippy::type_complexity)]
pub fn gen_trace(
    w: &[Vec<u32x16>],
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    Vec<Vec<u32x16>>,
) {
    let _span = span!(Level::INFO, "Generation").entered();
    let simd_size = w[0].len();

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
        .take(W_COLUMNS)
        .for_each(|(i, eval)| {
            *eval = w[i].clone();
        });
    lookup_data
        .iter_mut()
        .enumerate()
        .take(W_COLUMNS)
        .for_each(|(i, eval)| {
            *eval = w[i].clone();
        });

    for t in 0..64 {
        let index = W_COLUMNS + t * COL_PER_ROUND;
        let interaction_index = W_COLUMNS + t * INTERACTION_COL_PER_ROUND;

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
        let k_low = k[2 * t];
        let k_high = k[2 * t + 1];

        for simd_row in 0..simd_size {
            // Load W value
            let w_low = evals[2 * t][simd_row];
            let w_high = evals[2 * t + 1][simd_row];

            // BIG_SIGMA1
            // Decomposition over I0
            let e_i0_low = e_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let e_i0_high = e_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let sigma1 = big_sigma1_u32x16(e_i0_low + (e_i0_high << 16));
            let sigma1_o0_low = sigma1 & u32x16::splat(BigSigma1::O0_L);
            let sigma1_o0_high = (sigma1 >> 16) & u32x16::splat(BigSigma1::O0_H);
            let sigma1_o20 = sigma1 & u32x16::splat(BigSigma1::O2);
            let sigma1_o20_pext = pext_u32x16(sigma1_o20, BigSigma1::O2);

            // Decomposition over I1
            let e_i1_low = e_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let e_i1_high = e_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let sigma1 = big_sigma1_u32x16(e_i1_low + (e_i1_high << 16));
            let sigma1_o1_low = sigma1 & u32x16::splat(BigSigma1::O1_L);
            let sigma1_o1_high = (sigma1 >> 16) & u32x16::splat(BigSigma1::O1_H);
            let sigma1_o21 = sigma1 & u32x16::splat(BigSigma1::O2);
            let sigma1_o21_pext = pext_u32x16(sigma1_o21, BigSigma1::O2);

            // XOR the two O2 values
            let sigma1_o2 = sigma1_o20 ^ sigma1_o21;
            let sigma1_o2_low = sigma1_o2 & u32x16::splat(0xffff);
            let sigma1_o2_high = sigma1_o2 >> 16;

            // Compute sigma output
            let sigma1_low = sigma1_o0_low + sigma1_o1_low + sigma1_o2_low;
            let sigma1_high = sigma1_o0_high + sigma1_o1_high + sigma1_o2_high;

            // CH
            // left side
            let f_i0_low = f_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let f_i0_high = f_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let f_i1_low = f_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let f_i1_high = f_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let ch_left_i0_l = ch_left_u32x16(e_i0_low, f_i0_low);
            let ch_left_i0_h = ch_left_u32x16(e_i0_high, f_i0_high);
            let ch_left_i1_l = ch_left_u32x16(e_i1_low, f_i1_low);
            let ch_left_i1_h = ch_left_u32x16(e_i1_high, f_i1_high);

            // right side
            let g_i0_low = g_low[simd_row] & u32x16::splat(BigSigma1::I0_L);
            let g_i0_high = g_high[simd_row] & u32x16::splat(BigSigma1::I0_H);
            let g_i1_low = g_low[simd_row] & u32x16::splat(BigSigma1::I1_L);
            let g_i1_high = g_high[simd_row] & u32x16::splat(BigSigma1::I1_H);
            let ch_right_i0_l = ch_right_u32x16(e_i0_low, g_i0_low);
            let ch_right_i0_h = ch_right_u32x16(e_i0_high, g_i0_high);
            let ch_right_i1_l = ch_right_u32x16(e_i1_low, g_i1_low);
            let ch_right_i1_h = ch_right_u32x16(e_i1_high, g_i1_high);

            let ch_low = ch_left_i0_l + ch_left_i1_l + ch_right_i0_l + ch_right_i1_l;
            let ch_high = ch_left_i0_h + ch_left_i1_h + ch_right_i0_h + ch_right_i1_h;

            // BIG_SIGMA0
            // Decomposition over I0
            let a_i0_low = a_low[simd_row] & u32x16::splat(BigSigma0::I0_L);
            let a_i0_high_0 = a_high[simd_row] & u32x16::splat(BigSigma0::I0_H0);
            let a_i0_high_1 = (a_high[simd_row] >> 8) & u32x16::splat(BigSigma0::I0_H1);

            let sigma0 = big_sigma0_u32x16(a_i0_low + (a_i0_high_0 << 16) + (a_i0_high_1 << 24));
            let sigma0_o0_low = sigma0 & u32x16::splat(BigSigma0::O0_L);
            let sigma0_o0_high = (sigma0 >> 16) & u32x16::splat(BigSigma0::O0_H);
            let sigma0_o20 = sigma0 & u32x16::splat(BigSigma0::O2);
            let sigma0_o20_pext = pext_u32x16(sigma0_o20, BigSigma0::O2);

            // Decomposition over I1
            let a_i1_low_0 = a_low[simd_row] & u32x16::splat(BigSigma0::I1_L0);
            let a_i1_low_1 = (a_low[simd_row] >> 8) & u32x16::splat(BigSigma0::I1_L1);
            let a_i1_high = a_high[simd_row] & u32x16::splat(BigSigma0::I1_H);

            let sigma0 = big_sigma0_u32x16(a_i1_low_0 + (a_i1_low_1 << 8) + (a_i1_high << 16));
            let sigma0_o1_low = sigma0 & u32x16::splat(BigSigma0::O1_L);
            let sigma0_o1_high = (sigma0 >> 16) & u32x16::splat(BigSigma0::O1_H);
            let sigma0_o21 = sigma0 & u32x16::splat(BigSigma0::O2);
            let sigma0_o21_pext = pext_u32x16(sigma0_o21, BigSigma0::O2);

            // XOR the two O2 values
            let sigma0_o2 = sigma0_o20 ^ sigma0_o21;
            let sigma0_o2_low = sigma0_o2 & u32x16::splat(0xffff);
            let sigma0_o2_high = sigma0_o2 >> 16;

            // Compute sigma0 output
            let sigma0_low = sigma0_o0_low + sigma0_o1_low + sigma0_o2_low;
            let sigma0_high = sigma0_o0_high + sigma0_o1_high + sigma0_o2_high;

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
            let maj_low = maj_i0_low + maj_i1_low_0 + maj_i1_low_1;
            let maj_high = maj_i0_high_0 + maj_i0_high_1 + maj_i1_high;

            let temp1_low = h_low[simd_row] + sigma1_low + ch_low + k_low + w_low;
            let temp1_high = h_high[simd_row] + sigma1_high + ch_high + k_high + w_high;
            let temp2_low = sigma0_low + maj_low;
            let temp2_high = sigma0_high + maj_high;

            let e_carry_low = (temp1_low + d_low[simd_row]) >> 16;
            let e_carry_high = (temp1_high + d_high[simd_row] + e_carry_low) >> 16;
            let a_carry_low = (temp1_low + temp2_low) >> 16;
            let a_carry_high = (temp2_high + temp2_high + a_carry_low) >> 16;

            let trace_values: [u32x16; COL_PER_ROUND] = [
                // BIG_SIGMA1
                e_i0_low,
                e_i0_high,
                sigma1_o0_low,
                sigma1_o0_high,
                sigma1_o20_pext,
                sigma1_o1_low,
                sigma1_o1_high,
                sigma1_o21_pext,
                sigma1_o2_low,
                sigma1_o2_high,
                // CH
                f_i0_low,
                f_i0_high,
                ch_left_i0_l,
                ch_left_i0_h,
                ch_left_i1_l,
                ch_left_i1_h,
                g_i0_low,
                g_i0_high,
                ch_right_i0_l,
                ch_right_i0_h,
                ch_right_i1_l,
                ch_right_i1_h,
                // BIG_SIGMA0
                a_i0_high_0,
                a_i0_high_1,
                a_i1_low_0,
                a_i1_low_1,
                sigma0_o0_low,
                sigma0_o0_high,
                sigma0_o20_pext,
                sigma0_o1_low,
                sigma0_o1_high,
                sigma0_o21_pext,
                sigma0_o2_low,
                sigma0_o2_high,
                // MAJ
                b_i0_high_0,
                b_i0_high_1,
                b_i1_low_0,
                b_i1_low_1,
                c_i0_high_0,
                c_i0_high_1,
                c_i1_low_0,
                c_i1_low_1,
                maj_i0_low,
                maj_i0_high_0,
                maj_i0_high_1,
                maj_i1_low_0,
                maj_i1_low_1,
                maj_i1_high,
                // ADD
                e_carry_low,
                e_carry_high,
                a_carry_low,
                a_carry_high,
            ];
            for (i, value) in trace_values.iter().enumerate() {
                evals[index + i].push(*value);
            }

            let interaction_values: [u32x16; INTERACTION_COL_PER_ROUND] = [
                // BIG_SIGMA1
                e_i0_low,
                e_i0_high,
                sigma1_o0_low,
                sigma1_o0_high,
                sigma1_o20_pext,
                e_i1_low,
                e_i1_high,
                sigma1_o1_low,
                sigma1_o1_high,
                sigma1_o21_pext,
                sigma1_o2_low,
                sigma1_o2_high,
                // CH
                f_i0_low,
                f_i0_high,
                f_i1_low,
                f_i1_high,
                ch_left_i0_l,
                ch_left_i0_h,
                ch_left_i1_l,
                ch_left_i1_h,
                g_i0_low,
                g_i0_high,
                g_i1_low,
                g_i1_high,
                ch_right_i0_l,
                ch_right_i0_h,
                ch_right_i1_l,
                ch_right_i1_h,
                // BIG_SIGMA0
                a_i0_low,
                a_i0_high_0,
                a_i0_high_1,
                a_i1_low_0,
                a_i1_low_1,
                a_i1_high,
                sigma0_o0_low,
                sigma0_o0_high,
                sigma0_o20_pext,
                sigma0_o1_low,
                sigma0_o1_high,
                sigma0_o21_pext,
                sigma0_o2_low,
                sigma0_o2_high,
                // MAJ
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
                // ADD
                e_carry_low,
                e_carry_high,
                a_carry_low,
                a_carry_high,
            ];
            for (i, value) in interaction_values.iter().enumerate() {
                lookup_data[interaction_index + i].push(*value);
            }
        }

        let w_low = evals[2 * t].clone();
        let w_high = evals[2 * t + 1].clone();
        let sigma1_o0_low = evals[trace_index(t, ColumnIndex::sigma1_o0_low)].clone();
        let sigma1_o0_high = evals[trace_index(t, ColumnIndex::sigma1_o0_high)].clone();
        let sigma1_o1_low = evals[trace_index(t, ColumnIndex::sigma1_o1_low)].clone();
        let sigma1_o1_high = evals[trace_index(t, ColumnIndex::sigma1_o1_high)].clone();
        let sigma1_o2_low = evals[trace_index(t, ColumnIndex::sigma1_o2_low)].clone();
        let sigma1_o2_high = evals[trace_index(t, ColumnIndex::sigma1_o2_high)].clone();
        let ch_left_i0_l = evals[trace_index(t, ColumnIndex::ch_left_i0_l)].clone();
        let ch_left_i0_h = evals[trace_index(t, ColumnIndex::ch_left_i0_h)].clone();
        let ch_left_i1_l = evals[trace_index(t, ColumnIndex::ch_left_i1_l)].clone();
        let ch_left_i1_h = evals[trace_index(t, ColumnIndex::ch_left_i1_h)].clone();
        let ch_right_i0_l = evals[trace_index(t, ColumnIndex::ch_right_i0_l)].clone();
        let ch_right_i0_h = evals[trace_index(t, ColumnIndex::ch_right_i0_h)].clone();
        let ch_right_i1_l = evals[trace_index(t, ColumnIndex::ch_right_i1_l)].clone();
        let ch_right_i1_h = evals[trace_index(t, ColumnIndex::ch_right_i1_h)].clone();
        let sigma0_o0_low = evals[trace_index(t, ColumnIndex::sigma0_o0_low)].clone();
        let sigma0_o0_high = evals[trace_index(t, ColumnIndex::sigma0_o0_high)].clone();
        let sigma0_o1_low = evals[trace_index(t, ColumnIndex::sigma0_o1_low)].clone();
        let sigma0_o1_high = evals[trace_index(t, ColumnIndex::sigma0_o1_high)].clone();
        let sigma0_o2_low = evals[trace_index(t, ColumnIndex::sigma0_o2_low)].clone();
        let sigma0_o2_high = evals[trace_index(t, ColumnIndex::sigma0_o2_high)].clone();
        let maj_i0_low = evals[trace_index(t, ColumnIndex::maj_i0_low)].clone();
        let maj_i0_high_0 = evals[trace_index(t, ColumnIndex::maj_i0_high_0)].clone();
        let maj_i0_high_1 = evals[trace_index(t, ColumnIndex::maj_i0_high_1)].clone();
        let maj_i1_low_0 = evals[trace_index(t, ColumnIndex::maj_i1_low_0)].clone();
        let maj_i1_low_1 = evals[trace_index(t, ColumnIndex::maj_i1_low_1)].clone();
        let maj_i1_high = evals[trace_index(t, ColumnIndex::maj_i1_high)].clone();

        let sigma1_high: Vec<u32x16> = izip!(sigma1_o0_high, sigma1_o1_high, sigma1_o2_high)
            .map(|(a, b, c)| a + b + c)
            .collect();
        let sigma1_low: Vec<u32x16> = izip!(sigma1_o0_low, sigma1_o1_low, sigma1_o2_low)
            .map(|(a, b, c)| a + b + c)
            .collect();

        let ch_low: Vec<u32x16> = izip!(ch_left_i0_l, ch_left_i1_l, ch_right_i0_l, ch_right_i1_l)
            .map(|(a, b, c, d)| a + b + c + d)
            .collect();
        let ch_high: Vec<u32x16> = izip!(ch_left_i0_h, ch_left_i1_h, ch_right_i0_h, ch_right_i1_h)
            .map(|(a, b, c, d)| a + b + c + d)
            .collect();

        let sigma0_high: Vec<u32x16> = izip!(sigma0_o0_high, sigma0_o1_high, sigma0_o2_high)
            .map(|(a, b, c)| a + b + c)
            .collect();
        let sigma0_low: Vec<u32x16> = izip!(sigma0_o0_low, sigma0_o1_low, sigma0_o2_low)
            .map(|(a, b, c)| a + b + c)
            .collect();

        let maj_high: Vec<u32x16> = izip!(maj_i0_high_0, maj_i0_high_1, maj_i1_high)
            .map(|(a, b, c)| a + b + c)
            .collect();
        let maj_low: Vec<u32x16> = izip!(maj_i0_low, maj_i1_low_0, maj_i1_low_1)
            .map(|(a, b, c)| a + b + c)
            .collect();

        let temp1_high: Vec<u32x16> = izip!(h_high, sigma1_high, ch_high, w_high)
            .map(|(a, b, c, d)| a + b + c + d + k_high)
            .collect();
        let temp1_low: Vec<u32x16> = izip!(h_low, sigma1_low, ch_low, w_low)
            .map(|(a, b, c, d)| a + b + c + d + k_low)
            .collect();

        let temp2_high: Vec<u32x16> = izip!(sigma0_high, maj_high).map(|(a, b)| a + b).collect();
        let temp2_low: Vec<u32x16> = izip!(sigma0_low, maj_low).map(|(a, b)| a + b).collect();

        let e_low: Vec<u32x16> = izip!(d_low.clone(), temp1_low.clone())
            .map(|(a, b)| (a + b) & u32x16::splat(0xffff))
            .collect();
        let e_carry_low: Vec<u32x16> = izip!(d_low.clone(), temp1_low.clone())
            .map(|(a, b)| (a + b) >> 16)
            .collect();
        let e_high: Vec<u32x16> = izip!(d_high, temp1_high.clone(), e_carry_low)
            .map(|(a, b, c)| (a + b + c) & u32x16::splat(0xffff))
            .collect();

        let a_low: Vec<u32x16> = izip!(temp1_low.clone(), temp2_low.clone())
            .map(|(a, b)| (a + b) & u32x16::splat(0xffff))
            .collect();
        let a_carry_low: Vec<u32x16> = izip!(temp1_low.clone(), temp2_low.clone())
            .map(|(a, b)| (a + b) >> 16)
            .collect();
        let a_high = izip!(temp1_high.clone(), temp2_high.clone(), a_carry_low.clone())
            .map(|(a, b, c)| (a + b + c) & u32x16::splat(0xffff))
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

#[cfg(test)]
mod tests {

    use crate::components::scheduling::{gen_trace as gen_schedule, CHUNK_COLUMNS};

    use super::*;
    use crate::sha256::{big_sigma0, big_sigma1, ch_left, ch_right, maj, process_chunk};

    #[test]
    fn test_gen_trace_columns_count() {
        let w = vec![vec![u32x16::splat(0); LOG_N_LANES as usize]; W_COLUMNS];
        let (trace, _) = gen_trace(&w);
        assert_eq!(trace.len(), N_COLUMNS);
    }

    #[test]
    fn test_gen_trace_values() {
        let log_size = LOG_N_LANES;
        let (mut schedule, _) = gen_schedule(log_size);
        let schedule: [Vec<u32x16>; W_COLUMNS] = schedule
            .drain(0..W_COLUMNS)
            .map(|eval| {
                eval.data
                    .clone()
                    .into_iter()
                    .map(|x| x.into_simd())
                    .collect()
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let (mut trace, _) = gen_trace(&schedule);
        let chunk = &schedule[0..CHUNK_COLUMNS];

        // for row in 0..size {
        //     let [h, g, f, e, d, c, b, a] = H;
        //     for t in 0..64 {
        //         let w_low = schedule[2 * t][row];
        //         let w_high = schedule[2 * t + 1][row];
        //         let w = w_low + (w_high << 16);
        //         let temp1 = h + big_sigma1(e) + ch_left(e, f) + ch_right(e, g) + w + K[t];
        //         let temp2 = big_sigma0(a) + maj(a, b, c);
        //         let h = g;
        //         let g = f;
        //         let f = e;
        //         let e = d + temp1;
        //         let d = c;
        //         let c = b;
        //         let b = a;
        //         let a = temp1 + temp2;
        //     }
        //     let chunk: [u32; 64] = std::array::from_fn(|i| schedule[i][row]);
        //     let result = [0; 8];
        //     result[0] = H[0] + a;
        //     result[1] = H[1] + b;
        //     result[2] = H[2] + c;
        //     result[3] = H[3] + d;
        //     result[4] = H[4] + e;
        //     result[5] = H[5] + f;
        //     result[6] = H[6] + g;
        //     result[7] = H[7] + h;
        //     assert_eq!(process_chunk(chunk, H), result);
        // }
    }
}
