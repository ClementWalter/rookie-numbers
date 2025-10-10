//! SHA-256 functions

use core::simd::Simd;
use std::simd::u32x16;

#[inline(always)]
fn rotr_u32x16(x: u32x16, n: u32) -> u32x16 {
    let n = Simd::splat(n);
    let n2 = Simd::splat(32 - n[0]); // n is constant here, so n[0] is fine
    (x >> n) | (x << n2)
}

#[inline(always)]
pub fn small_sigma0_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 7) ^ rotr_u32x16(x, 18) ^ (x >> Simd::splat(3))
}

#[inline(always)]
pub fn small_sigma1_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 17) ^ rotr_u32x16(x, 19) ^ (x >> Simd::splat(10))
}

#[inline(always)]
pub fn big_sigma0_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 2) ^ rotr_u32x16(x, 13) ^ rotr_u32x16(x, 22)
}

#[inline(always)]
pub fn big_sigma1_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 6) ^ rotr_u32x16(x, 11) ^ rotr_u32x16(x, 25)
}

#[inline(always)]
pub fn ch_left_u32x16(e: u32x16, f: u32x16) -> u32x16 {
    e & f
}

#[inline(always)]
pub fn ch_right_u32x16(e: u32x16, g: u32x16) -> u32x16 {
    (!e) & g
}

#[inline(always)]
pub fn maj_u32x16(a: u32x16, b: u32x16, c: u32x16) -> u32x16 {
    (a & b) ^ (a & c) ^ (b & c)
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

#[cfg(test)]
mod tests {
    use itertools::izip;

    use super::*;

    #[test]
    fn test_rotr_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        assert_eq!(
            base.map(|x| x.rotate_right(7)),
            rotr_u32x16(u32x16::from_array(base), 7).to_array()
        );
    }

    #[test]
    fn test_small_sigma0_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(small_sigma0),
            small_sigma0_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_small_sigma1_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(small_sigma1),
            small_sigma1_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_big_sigma0_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(big_sigma0),
            big_sigma0_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_big_sigma1_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(big_sigma1),
            big_sigma1_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_ch_left_u32x16() {
        let x: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let y: [u32; 16] = x
            .into_iter()
            .rev()
            .collect::<Vec<u32>>()
            .try_into()
            .unwrap();
        assert_eq!(
            izip!(x, y)
                .map(|(_x, _y)| ch_left(_x, _y))
                .collect::<Vec<u32>>(),
            ch_left_u32x16(u32x16::from_array(x), u32x16::from_array(y)).to_array()
        );
    }

    #[test]
    fn test_ch_right_u32x16() {
        let x: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let y: [u32; 16] = x
            .into_iter()
            .rev()
            .collect::<Vec<u32>>()
            .try_into()
            .unwrap();
        assert_eq!(
            izip!(x, y)
                .map(|(_x, _y)| ch_right(_x, _y))
                .collect::<Vec<u32>>(),
            ch_right_u32x16(u32x16::from_array(x), u32x16::from_array(y)).to_array()
        );
    }
}
