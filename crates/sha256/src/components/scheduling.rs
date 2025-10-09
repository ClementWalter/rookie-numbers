//! The scheduling component is responsible for proving the scheduling part of sha256.
//!
//! This is, 48 iterations of the main loop, doing sigma0, sigma1 and adding the
//! results with some previous W values to the current W value.
//! The W array is actually write-once in the sense that no value is updated.
//! For the AIR, we just need to define an encoding for a row
//! that would contained all the required data, and access the required values with
//! deterministically computed indices.
//!
//! The AIR works as follows:
//!
//! - row[0: 128] = the resulting W array, decomposed in low and high parts
//!
//! For sigma0, only indexes from 1..=48 are used. Each operation requires
//! to hint the decomposition over I0 and the 3 output values, so 5 columns, or 5*48=240 columns.
//!
//! - row[128:(128 + 240)] = sigma0 values
//!
//! For sigma1, only indexes from 14..=61 are used. Each operation requires
//! to hint the decomposition over I0 and the 3 output values, so 5 columns, or 5*48=240 columns.
//!
//! - row[(128 + 240):(128 + 2*240)] = sigma1 values
//!
//! The final addition at each iteration requires to hint the 2 carries. There are 48 iterations,
//! so 2*48=96 columns.
//!
//! - row[(128 + 2*240):(128 + 2*240 + 96)] = the carries
//!
//! In short, the main trace of the AIR is then defined as follows:
//! - row[0: 128] = W
//! - row[128:368] = sigma0 values
//! - row[368:608] = sigma1 values
//! - row[608:704] = the carries

use crate::partitions::{pext_u32x16, Sigma0, Sigma1};
use crate::sha256::{small_sigma0_u32x16, small_sigma1_u32x16};
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

use std::simd::u32x16;

use crate::relations::LookupData;

const CHUNK_COLUMNS: usize = 16 * 2;
const W_COLUMNS: usize = CHUNK_COLUMNS + 2 * 48;
const SIGMA0_COLUMNS: usize = 48 * 10;
const SIGMA1_COLUMNS: usize = 48 * 10;
const CARRIES_COLUMNS: usize = 48 * 2;
const N_COLUMNS: usize = W_COLUMNS + SIGMA0_COLUMNS + SIGMA1_COLUMNS + CARRIES_COLUMNS;

fn generate_simd_sequence_bulk(start: usize, len: usize) -> Vec<u32x16> {
    let mut buffer: Vec<u32> = Vec::with_capacity(len);
    for i in 0..len {
        buffer.push((start + i) as u32);
    }
    let simd_ptr = buffer.as_ptr() as *const u32x16;
    let simd_vec = unsafe { Vec::from_raw_parts(simd_ptr as *mut u32x16, len / 16, len / 16) };
    std::mem::forget(buffer); // Prevent double-free
    simd_vec
}

pub fn gen_trace(
    log_size: u32,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    LookupData,
) {
    let _span = span!(Level::INFO, "Generation").entered();
    assert!(log_size >= LOG_N_LANES);
    let simd_log_size = log_size - 4;
    let simd_size = 1 << simd_log_size;

    // Initialize vec for all groups of columns
    let mut w: Vec<Vec<u32x16>> = Vec::with_capacity(W_COLUMNS);
    let mut sigma0_columns: Vec<Vec<u32x16>> = Vec::with_capacity(SIGMA0_COLUMNS);
    let mut sigma1_columns: Vec<Vec<u32x16>> = Vec::with_capacity(SIGMA1_COLUMNS);
    let mut carries_columns: Vec<Vec<u32x16>> = Vec::with_capacity(CARRIES_COLUMNS);

    // Generate random inputs
    for i in 0..CHUNK_COLUMNS {
        w.push(generate_simd_sequence_bulk(i, 1 << log_size));
    }

    for t in 16..64 {
        // Generate new w values
        let mut new_w_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut new_w_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        // Generate 10 sigma0 values per round
        let mut sigma0_i0_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_i0_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o0_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o0_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o20_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o1_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o1_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o21_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o2_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma0_o2_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        // Generate 10 sigma1 values per round
        let mut sigma1_i0_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_i0_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o0_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o0_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o20_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o1_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o1_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o21_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o2_low_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut sigma1_o2_high_col: Vec<u32x16> = Vec::with_capacity(simd_size);
        // Generate 2 carries per round
        let mut carry_0: Vec<u32x16> = Vec::with_capacity(simd_size);
        let mut carry_1: Vec<u32x16> = Vec::with_capacity(simd_size);

        for simd_row in 0..simd_size {
            // Load the W values
            let w_16_low = w[2 * (t - 16)][simd_row];
            let w_16_high = w[2 * (t - 16) + 1][simd_row];
            let w_15_low = w[2 * (t - 15)][simd_row];
            let w_15_high = w[2 * (t - 15) + 1][simd_row];
            let w_7_low = w[2 * (t - 7)][simd_row];
            let w_7_high = w[2 * (t - 7) + 1][simd_row];
            let w_2_low = w[2 * (t - 2)][simd_row];
            let w_2_high = w[2 * (t - 2) + 1][simd_row];

            // SIGMA0
            // Decomposition over I0
            let w_15_i0_low = w_15_low & u32x16::splat(Sigma0::I0_L);
            sigma0_i0_low_col.push(w_15_i0_low);
            let w_15_i0_high = w_15_high & u32x16::splat(Sigma0::I0_H);
            sigma0_i0_high_col.push(w_15_i0_high);
            let (sigma0_low, sigma0_high) = small_sigma0_u32x16(w_15_i0_low, w_15_i0_high);
            let sigma0_o0_low = sigma0_low & u32x16::splat(Sigma0::O0_L);
            sigma0_o0_low_col.push(sigma0_o0_low);
            let sigma0_o0_high = sigma0_high & u32x16::splat(Sigma0::O0_H);
            sigma0_o0_high_col.push(sigma0_o0_high);
            let sigma0_o20 = (sigma0_low + (sigma0_high << 16)) & u32x16::splat(Sigma0::O2);
            let sigma0_o20_pext = pext_u32x16(sigma0_o20, Sigma0::O2);
            sigma0_o20_col.push(sigma0_o20_pext);

            // Decomposition over I1
            let w_15_i1_low = w_15_low & u32x16::splat(Sigma0::I1_L);
            let w_15_i1_high = w_15_high & u32x16::splat(Sigma0::I1_H);
            let (sigma0_low, sigma0_high) = small_sigma0_u32x16(w_15_i1_low, w_15_i1_high);
            let sigma0_o1_low = sigma0_low & u32x16::splat(Sigma0::O1_L);
            sigma0_o1_low_col.push(sigma0_o1_low);
            let sigma0_o1_high = sigma0_high & u32x16::splat(Sigma0::O1_H);
            sigma0_o1_high_col.push(sigma0_o1_high);
            let sigma0_o21 = (sigma0_low + (sigma0_high << 16)) & u32x16::splat(Sigma0::O2);
            let sigma0_o21_pext = pext_u32x16(sigma0_o21, Sigma0::O2);
            sigma0_o21_col.push(sigma0_o21_pext);

            // XOR the two O2 values
            let sigma0_o2 = sigma0_o20 ^ sigma0_o21;
            let sigma0_o2_low = sigma0_o2 & u32x16::splat(0xffff);
            let sigma0_o2_high = sigma0_o2 >> 16;
            sigma0_o2_low_col.push(sigma0_o2_low);
            sigma0_o2_high_col.push(sigma0_o2_high);

            // Compute sigma0 output
            let sigma0_low = sigma0_o0_low + sigma0_o1_low + sigma0_o2_low;
            let sigma0_high = sigma0_o0_high + sigma0_o1_high + sigma0_o2_high;

            // SIGMA1
            // Decomposition over I0
            let w_7_i0_low = w_7_low & u32x16::splat(Sigma1::I0_L);
            sigma1_i0_low_col.push(w_7_i0_low);
            let w_7_i0_high = w_7_high & u32x16::splat(Sigma1::I0_H);
            sigma1_i0_high_col.push(w_7_i0_high);
            let (sigma1_low, sigma1_high) = small_sigma1_u32x16(w_7_i0_low, w_7_i0_high);
            let sigma1_o0_low = sigma1_low & u32x16::splat(Sigma1::O0_L);
            sigma1_o0_low_col.push(sigma1_o0_low);
            let sigma1_o0_high = sigma1_high & u32x16::splat(Sigma1::O0_H);
            sigma1_o0_high_col.push(sigma1_o0_high);
            let sigma1_o20 = (sigma1_low + (sigma1_high << 16)) & u32x16::splat(Sigma1::O2);
            let sigma1_o20_pext = pext_u32x16(sigma1_o20, Sigma1::O2);
            sigma1_o20_col.push(sigma1_o20_pext);

            // Decomposition over I1
            let w_7_i1_low = w_7_low & u32x16::splat(Sigma1::I1_L);
            let w_7_i1_high = w_7_high & u32x16::splat(Sigma1::I1_H);
            let (sigma1_low, sigma1_high) = small_sigma1_u32x16(w_7_i1_low, w_7_i1_high);
            let sigma1_o1_low = sigma1_low & u32x16::splat(Sigma1::O1_L);
            sigma1_o1_low_col.push(sigma1_o1_low);
            let sigma1_o1_high = sigma1_high & u32x16::splat(Sigma1::O1_H);
            sigma1_o1_high_col.push(sigma1_o1_high);
            let sigma1_o21 = (sigma1_low + (sigma1_high << 16)) & u32x16::splat(Sigma1::O2);
            let sigma1_o21_pext = pext_u32x16(sigma1_o21, Sigma1::O2);
            sigma1_o21_col.push(sigma1_o21_pext);

            // XOR the two O2 values
            let sigma1_o2 = sigma1_o20 ^ sigma1_o21;
            let sigma1_o2_low = sigma1_o2 & u32x16::splat(0xffff);
            let sigma1_o2_high = sigma1_o2 >> 16;
            sigma1_o2_low_col.push(sigma1_o2_low);
            sigma1_o2_high_col.push(sigma1_o2_high);

            // Compute sigma1 output
            let sigma1_low = sigma1_o0_low + sigma1_o1_low + sigma1_o2_low;
            let sigma1_high = sigma1_o0_high + sigma1_o1_high + sigma1_o2_high;

            // Compute the final output
            let carry_low = (w_16_low + sigma0_low + w_2_low + sigma1_low) >> 16;
            let carry_high = (w_16_high + sigma0_high + w_2_high + sigma1_high + carry_low) >> 16;
            carry_0.push(carry_low);
            carry_1.push(carry_high);
            let new_w_low = w_16_low + sigma0_low + w_2_low + sigma1_low - (carry_low << 16);
            let new_w_high =
                w_16_high + sigma0_high + w_2_high + sigma1_high + carry_low - (carry_high << 16);
            new_w_low_col.push(new_w_low);
            new_w_high_col.push(new_w_high);
        }

        w.push(new_w_low_col);
        w.push(new_w_high_col);
        sigma0_columns.push(sigma0_i0_low_col);
        sigma0_columns.push(sigma0_i0_high_col);
        sigma0_columns.push(sigma0_o0_low_col);
        sigma0_columns.push(sigma0_o0_high_col);
        sigma0_columns.push(sigma0_o20_col);
        sigma0_columns.push(sigma0_o1_low_col);
        sigma0_columns.push(sigma0_o1_high_col);
        sigma0_columns.push(sigma0_o21_col);
        sigma0_columns.push(sigma0_o2_low_col);
        sigma0_columns.push(sigma0_o2_high_col);
        sigma1_columns.push(sigma1_i0_low_col);
        sigma1_columns.push(sigma1_i0_high_col);
        sigma1_columns.push(sigma1_o0_low_col);
        sigma1_columns.push(sigma1_o0_high_col);
        sigma1_columns.push(sigma1_o20_col);
        sigma1_columns.push(sigma1_o1_low_col);
        sigma1_columns.push(sigma1_o1_high_col);
        sigma1_columns.push(sigma1_o21_col);
        sigma1_columns.push(sigma1_o2_low_col);
        sigma1_columns.push(sigma1_o2_high_col);
        carries_columns.push(carry_0);
        carries_columns.push(carry_1);
    }

    let lookup_data = LookupData::new(log_size);

    w.extend(sigma0_columns);
    w.extend(sigma1_columns);
    w.extend(carries_columns);
    let domain = CanonicCoset::new(log_size).circle_domain();
    let trace = w
        .into_iter()
        .map(|col| {
            let eval = BaseColumn::from_simd(
                col.into_iter()
                    .map(|simd_chunk| unsafe { PackedM31::from_simd_unchecked(simd_chunk) })
                    .collect(),
            );
            CircleEvaluation::new(domain, eval)
        })
        .collect();

    (trace, lookup_data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gen_trace() {
        let (trace, _) = gen_trace(LOG_N_LANES * 2);
        assert_eq!(trace.len(), N_COLUMNS);
    }
}
