//! SHA-256 functions

use std::simd::u32x16;

/// Right rotate a 32-bit value expressed as two 16-bit values (low, high)
#[inline(always)]
pub fn rotr_u32x16<const K: u32>(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let s = u32x16::splat(K);
    let s_inv = u32x16::splat(32 - K);
    let lo_r = (lo >> s) | (hi << s_inv);
    let hi_r = (hi >> s) | (lo << s_inv);
    (lo_r, hi_r)
}

/// Right shift a 32-bit value expressed as two 16-bit values (low, high)
#[inline(always)]
pub fn shr_u32x16<const K: u32>(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let s = u32x16::splat(K);
    let s_inv = u32x16::splat(32 - K);
    ((lo >> s) | (hi << s_inv), hi >> s)
}

/// SHA-256 σ₀ function used in message scheduling
#[inline(always)]
pub fn small_sigma0_u32x16(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let (r7_lo, r7_hi) = rotr_u32x16::<7>(lo, hi);
    let (r18_lo, r18_hi) = rotr_u32x16::<18>(lo, hi);
    let (sh3_lo, sh3_hi) = shr_u32x16::<3>(lo, hi);
    (r7_lo ^ r18_lo ^ sh3_lo, r7_hi ^ r18_hi ^ sh3_hi)
}

/// SHA-256 σ₁ function used in message scheduling
#[inline(always)]
pub fn small_sigma1_u32x16(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let (r17_lo, r17_hi) = rotr_u32x16::<17>(lo, hi);
    let (r19_lo, r19_hi) = rotr_u32x16::<19>(lo, hi);
    let (sh10_lo, sh10_hi) = shr_u32x16::<10>(lo, hi);
    (r17_lo ^ r19_lo ^ sh10_lo, r17_hi ^ r19_hi ^ sh10_hi)
}

/// SHA-256 Σ₀ function used in the compression function
#[inline(always)]
pub fn big_sigma0_u32x16(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let (r2_lo, r2_hi) = rotr_u32x16::<2>(lo, hi);
    let (r13_lo, r13_hi) = rotr_u32x16::<13>(lo, hi);
    let (r22_lo, r22_hi) = rotr_u32x16::<22>(lo, hi);
    (r2_lo ^ r13_lo ^ r22_lo, r2_hi ^ r13_hi ^ r22_hi)
}

/// SHA-256 Σ₁ function used in the compression function
#[inline(always)]
pub fn big_sigma1_u32x16(lo: u32x16, hi: u32x16) -> (u32x16, u32x16) {
    let (r6_lo, r6_hi) = rotr_u32x16::<6>(lo, hi);
    let (r11_lo, r11_hi) = rotr_u32x16::<11>(lo, hi);
    let (r25_lo, r25_hi) = rotr_u32x16::<25>(lo, hi);
    (r6_lo ^ r11_lo ^ r25_lo, r6_hi ^ r11_hi ^ r25_hi)
}

pub const fn small_sigma0(x: u32) -> u32 {
    x.rotate_right(7) ^ x.rotate_right(18) ^ (x >> 3)
}

pub const fn small_sigma1(x: u32) -> u32 {
    x.rotate_right(17) ^ x.rotate_right(19) ^ (x >> 10)
}

pub const fn big_sigma0(x: u32) -> u32 {
    x.rotate_right(2) ^ x.rotate_right(13) ^ x.rotate_right(22)
}

pub const fn big_sigma1(x: u32) -> u32 {
    x.rotate_right(6) ^ x.rotate_right(11) ^ x.rotate_right(25)
}

pub const fn ch_left(e: u32, f: u32) -> u32 {
    e & f
}

pub const fn ch_right(e: u32, g: u32) -> u32 {
    (0xFFFFFFFF - e) & g
}

pub const fn maj(a: u32, b: u32, c: u32) -> u32 {
    (a & b) ^ (a & c) ^ (b & c)
}

pub fn process_chunk(chunk: [u32; 16], mut hash: [u32; 8]) -> [u32; 8] {
    let mut w: [u32; 64] = [0; 64];
    w[..16].copy_from_slice(&chunk);

    // Schedule
    for t in 16..64 {
        w[t] = w[t - 16] + small_sigma0(w[t - 15]) + w[t - 7] + small_sigma1(w[t - 2])
    }

    // Compression
    let mut a = hash[0];
    let mut b = hash[1];
    let mut c = hash[2];
    let mut d = hash[3];
    let mut e = hash[4];
    let mut f = hash[5];
    let mut g = hash[6];
    let mut h = hash[7];
    for &wt in w.iter() {
        let temp1 = h + big_sigma1(e) + ch_left(e, f) + ch_right(e, g) + wt;
        let temp2 = big_sigma0(a) + maj(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }
    hash[0] += a;
    hash[1] += b;
    hash[2] += c;
    hash[3] += d;
    hash[4] += e;
    hash[5] += f;
    hash[6] += g;
    hash[7] += h;
    hash
}
