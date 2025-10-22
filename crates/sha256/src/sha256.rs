//! SHA-256 functions

use core::simd::Simd;
use std::simd::u32x16;

pub const K: [u32; 64] = [
    1116352408, 1899447441, 3049323471, 3921009573, 961987163, 1508970993, 2453635748, 2870763221,
    3624381080, 310598401, 607225278, 1426881987, 1925078388, 2162078206, 2614888103, 3248222580,
    3835390401, 4022224774, 264347078, 604807628, 770255983, 1249150122, 1555081692, 1996064986,
    2554220882, 2821834349, 2952996808, 3210313671, 3336571891, 3584528711, 113926993, 338241895,
    666307205, 773529912, 1294757372, 1396182291, 1695183700, 1986661051, 2177026350, 2456956037,
    2730485921, 2820302411, 3259730800, 3345764771, 3516065817, 3600352804, 4094571909, 275423344,
    430227734, 506948616, 659060556, 883997877, 958139571, 1322822218, 1537002063, 1747873779,
    1955562222, 2024104815, 2227730452, 2361852424, 2428436474, 2756734187, 3204031479, 3329325298,
];

pub const H: [u32; 8] = [
    1779033703, 3144134277, 1013904242, 2773480762, 1359893119, 2600822924, 528734635, 1541459225,
];

pub const CHUNK_SIZE: usize = 32; // 16 u32 = 32 u16
pub const N_SCHEDULING_ROUNDS: usize = 48; // 16..48
pub const N_COMPRESSION_ROUNDS: usize = 64;

#[inline(always)]
fn rotr_u32x16(x: u32x16, n: u32) -> u32x16 {
    let n = Simd::splat(n);
    let n2 = Simd::splat(32 - n[0]); // n is constant here, so n[0] is fine
    (x >> n) | (x << n2)
}

#[inline(always)]
pub fn small_sigma_0_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 7) ^ rotr_u32x16(x, 18) ^ (x >> Simd::splat(3))
}

#[inline(always)]
pub fn small_sigma_1_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 17) ^ rotr_u32x16(x, 19) ^ (x >> Simd::splat(10))
}

#[inline(always)]
pub fn big_sigma_0_u32x16(x: u32x16) -> u32x16 {
    rotr_u32x16(x, 2) ^ rotr_u32x16(x, 13) ^ rotr_u32x16(x, 22)
}

#[inline(always)]
pub fn big_sigma_1_u32x16(x: u32x16) -> u32x16 {
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

#[inline(always)]
pub fn process_chunk_u32x16(chunk: [u32x16; 16], mut hash: [u32x16; 8]) -> [u32x16; 8] {
    let mut w: [u32x16; 64] = [u32x16::splat(0); 64];
    w[..16].copy_from_slice(&chunk);

    // Schedule
    for t in 16..64 {
        w[t] =
            w[t - 16] + small_sigma_0_u32x16(w[t - 15]) + w[t - 7] + small_sigma_1_u32x16(w[t - 2])
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
    for round in 0..64 {
        let temp1 = h
            + big_sigma_1_u32x16(e)
            + ch_left_u32x16(e, f)
            + ch_right_u32x16(e, g)
            + u32x16::splat(K[round])
            + w[round];
        let temp2 = big_sigma_0_u32x16(a) + maj_u32x16(a, b, c);
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

pub const fn small_sigma_0(x: u32) -> u32 {
    x.rotate_right(7) ^ x.rotate_right(18) ^ (x >> 3)
}

pub const fn small_sigma_1(x: u32) -> u32 {
    x.rotate_right(17) ^ x.rotate_right(19) ^ (x >> 10)
}

pub const fn big_sigma_0(x: u32) -> u32 {
    x.rotate_right(2) ^ x.rotate_right(13) ^ x.rotate_right(22)
}

pub const fn big_sigma_1(x: u32) -> u32 {
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
        w[t] = w[t - 16] + small_sigma_0(w[t - 15]) + w[t - 7] + small_sigma_1(w[t - 2])
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
    for round in 0..64 {
        let temp1 = h + big_sigma_1(e) + ch_left(e, f) + ch_right(e, g) + w[round] + K[round];
        let temp2 = big_sigma_0(a) + maj(a, b, c);
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
    fn test_reference_implementation() {
        use sha2::{Digest, Sha256};

        use crate::sha256::{process_chunk, H};

        let input = b"hello world";

        // Reference result using sha2 crate
        let reference = Sha256::digest(input);

        // Build padded "hello world" message (11 bytes)
        let mut msg = Vec::from(input);
        msg.push(0x80);
        while !(msg.len() + 8).is_multiple_of(64) {
            msg.push(0x00);
        }
        let bit_len: u64 = input.len() as u64 * 8;
        msg.extend_from_slice(&bit_len.to_be_bytes());
        assert_eq!(msg.len(), 64);

        // Convert bytes to [u32; 16]
        let mut chunk = [0u32; 16];
        for (i, word) in chunk.iter_mut().enumerate() {
            let bytes: [u8; 4] = msg[4 * i..4 * i + 4].try_into().unwrap();
            *word = u32::from_be_bytes(bytes);
        }

        // Process single chunk
        let hash = process_chunk(chunk, H);

        // Convert [u32; 8] hash state to [u8; 32]
        let mut result = [0u8; 32];
        for (i, word) in hash.iter().enumerate() {
            result[4 * i..4 * i + 4].copy_from_slice(&word.to_be_bytes());
        }

        assert_eq!(reference[..], result[..]);
    }

    #[test]
    fn test_rotr_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        assert_eq!(
            base.map(|x| x.rotate_right(7)),
            rotr_u32x16(u32x16::from_array(base), 7).to_array()
        );
    }

    #[test]
    fn test_small_sigma_0_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(small_sigma_0),
            small_sigma_0_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_small_sigma_1_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(small_sigma_1),
            small_sigma_1_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_big_sigma_0_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(big_sigma_0),
            big_sigma_0_u32x16(u32x16::from_array(base)).to_array()
        );
    }

    #[test]
    fn test_big_sigma_1_u32x16() {
        let base: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        assert_eq!(
            base.map(big_sigma_1),
            big_sigma_1_u32x16(u32x16::from_array(base)).to_array()
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

    #[test]
    fn test_process_chunk_u32x16() {
        let data: [u32; 16] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
        let chunk: [u32x16; 16] = std::array::from_fn(|_| u32x16::from_array(data));
        // First chunk is 1, 1, ..., 1
        // Second chunk is 2, 2, ..., 2
        // ...
        // Last chunk is 16, 16, ..., 16
        let chunk_array: [[u32; 16]; 16] =
            std::array::from_fn(|i| std::array::from_fn(|_| data[i]));
        let hash: [u32x16; 8] = std::array::from_fn(|i| u32x16::splat(H[i]));
        let hash_array: [[u32; 8]; 16] = std::array::from_fn(|_| H);
        let result = process_chunk_u32x16(chunk, hash);
        let result = std::array::from_fn(|i| std::array::from_fn(|j| result[j].to_array()[i]));
        let expected: [[u32; 8]; 16] =
            std::array::from_fn(|i| process_chunk(chunk_array[i], hash_array[i]));
        assert_eq!(result, expected);
    }
}
