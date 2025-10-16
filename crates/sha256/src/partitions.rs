use std::simd::{cmp::SimdPartialEq, num::SimdUint, simd_swizzle, u16x16, u16x32, u32x16, Simd};

/// This will be used to identify the function part in the relation.
#[repr(u32)]
pub enum Partition {
    O1,
    O2,
    O3,
    I0_L,
    I0_H,
    I1_L,
    I1_H,
    I0_H0,
    I0_H1,
    I1_L0,
    I1_L1,
    ADD4,
    ADD6,
    ADD7,
}

#[allow(non_snake_case)]
pub mod Sigma0 {
    pub const I0: u32 = 0b01001010101010101011010101010101;
    pub const I1: u32 = 0b10110101010101010100101010101010;
    pub const I0_L: u32 = 0b1011010101010101;
    pub const I1_L: u32 = 0b0100101010101010;
    pub const I0_H: u32 = 0b0100101010101010;
    pub const I1_H: u32 = 0b1011010101010101;
    pub const O0: u32 = 0b10101000000101010101000000101010;
    pub const O1: u32 = 0b01010000001010101010100000010101;
    pub const O2: u32 = 0b00000111110000000000011111000000;
    pub const O0_L: u32 = 0b0101000000101010;
    pub const O1_L: u32 = 0b1010100000010101;
    pub const O2_L: u32 = 0b0000011111000000;
    pub const O0_H: u32 = 0b1010100000010101;
    pub const O1_H: u32 = 0b0101000000101010;
    pub const O2_H: u32 = 0b0000011111000000;
}

#[allow(non_snake_case)]
pub mod Sigma1 {
    pub const I0: u32 = 0b10101011010101101010100001010101;
    pub const I1: u32 = 0b01010100101010010101011110101010;
    pub const I0_L: u32 = 0b1010100001010101;
    pub const I1_L: u32 = 0b0101011110101010;
    pub const I0_H: u32 = 0b1010101101010110;
    pub const I1_H: u32 = 0b0101010010101001;
    pub const O0: u32 = 0b01010100000010101001010100101010;
    pub const O1: u32 = 0b00101010110101010000101000010100;
    pub const O2: u32 = 0b10000001001000000110000011000001;
    pub const O0_L: u32 = 0b1001010100101010;
    pub const O1_L: u32 = 0b0000101000010100;
    pub const O2_L: u32 = 0b0110000011000001;
    pub const O0_H: u32 = 0b0101010000001010;
    pub const O1_H: u32 = 0b0010101011010101;
    pub const O2_H: u32 = 0b1000000100100000;
}

#[allow(non_snake_case)]
pub mod BigSigma0 {
    pub const I0: u32 = 0b11110000011111000000111110000011;
    pub const I1: u32 = 0b00001111100000111111000001111100;
    pub const I0_L: u32 = 0b0000111110000011;
    pub const I0_L0: u32 = 0b10000011;
    pub const I0_L1: u32 = 0b00001111;
    pub const I0_H: u32 = 0b1111000001111100;
    pub const I0_H0: u32 = 0b01111100;
    pub const I0_H1: u32 = 0b11110000;
    pub const I1_L: u32 = 0b1111000001111100;
    pub const I1_L0: u32 = 0b01111100;
    pub const I1_L1: u32 = 0b11110000;
    pub const I1_H: u32 = 0b0000111110000011;
    pub const I1_H0: u32 = 0b10000011;
    pub const I1_H1: u32 = 0b00001111;
    pub const O0: u32 = 0b01110000000111100000001111000000;
    pub const O1: u32 = 0b00000011110000000111000000011110;
    pub const O2: u32 = 0b10001100001000011000110000100001;
    pub const O0_L: u32 = 0b0000001111000000;
    pub const O1_L: u32 = 0b0111000000011110;
    pub const O2_L: u32 = 0b1000110000100001;
    pub const O0_H: u32 = 0b0111000000011110;
    pub const O1_H: u32 = 0b0000001111000000;
    pub const O2_H: u32 = 0b1000110000100001;
}

#[allow(non_snake_case)]
pub mod BigSigma1 {
    pub const I0: u32 = 0b10011000110001100110011000110001;
    pub const I1: u32 = 0b01100111001110011001100111001110;
    pub const I0_L: u32 = 0b0110011000110001;
    pub const I0_H: u32 = 0b1001100011000110;
    pub const I1_L: u32 = 0b1001100111001110;
    pub const I1_H: u32 = 0b0110011100111001;
    pub const O0: u32 = 0b01000010001000110001100010001000;
    pub const O1: u32 = 0b00011000100011001110011000100011;
    pub const O2: u32 = 0b10100101010100000000000101010100;
    pub const O0_L: u32 = 0b0001100010001000;
    pub const O1_L: u32 = 0b1110011000100011;
    pub const O2_L: u32 = 0b0000000101010100;
    pub const O0_H: u32 = 0b0100001000100011;
    pub const O1_H: u32 = 0b0001100010001100;
    pub const O2_H: u32 = 0b1010010101010000;
}

/// Generates all subsets of the given bitmask `mask`
pub struct SubsetIterator {
    bit_pos: [u8; 32],
    k: u8,
    idx_fwd: u32,   // forward index
    idx_back: u32,  // backward index (exclusive upper bound)
}

impl SubsetIterator {
    pub const fn new(mask: u32) -> Self {
        let mut bit_pos = [0u8; 32];
        let mut k = 0u8;
        let mut m = mask;
        let mut pos = 0u8;
        while m != 0 {
            if (m & 1) != 0 {
                bit_pos[k as usize] = pos;
                k += 1;
            }
            m >>= 1;
            pos += 1;
        }
        let total = 1u32 << k;
        Self {
            bit_pos,
            k,
            idx_fwd: 0,
            idx_back: total,
        }
    }

    #[inline]
    fn scatter(&self, c: u32) -> u32 {
        let mut out = 0;
        for i in 0..self.k as usize {
            if (c >> i) & 1 != 0 {
                out |= 1u32 << self.bit_pos[i];
            }
        }
        out
    }
}

impl Iterator for SubsetIterator {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx_fwd >= self.idx_back {
            return None;
        }
        let val = self.scatter(self.idx_fwd);
        self.idx_fwd += 1;
        Some(val)
    }
}

impl DoubleEndedIterator for SubsetIterator {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.idx_fwd >= self.idx_back {
            return None;
        }
        self.idx_back -= 1;
        Some(self.scatter(self.idx_back))
    }
}

#[inline]
pub const fn pext_u32(x: u32, mut mask: u32) -> u32 {
    // Extract bits from x where mask has 1s, packed to the low bits (LSB-first).
    let mut out = 0u32;
    let mut bb = 1u32;
    while mask != 0 {
        let ls = mask & mask.wrapping_neg(); // lowest set bit
        if x & ls != 0 {
            out |= bb;
        }
        mask ^= ls;
        bb <<= 1;
    }
    out
}

#[inline]
pub fn pext_u32x16(x: u32x16, mut mask: u32) -> u32x16 {
    // Extract bits from each lane of x where the scalar mask has 1s, packed to LSBs.
    let mut out = Simd::splat(0u32);
    let mut bb = Simd::splat(1u32);
    while mask != 0 {
        let ls = mask & mask.wrapping_neg(); // lowest set bit in the scalar mask
        let ls_v: u32x16 = Simd::splat(ls);
        let contrib = (x & ls_v)
            .simd_ne(Simd::splat(0u32))
            .select(bb, Simd::splat(0u32));
        out |= contrib;
        mask ^= ls;
        bb <<= Simd::splat(1u32);
    }
    out
}

#[inline(always)]
pub fn widen_add_u16x32(lows: u16x32, highs: u16x32) -> (u32x16, u32x16) {
    let a_lo: u16x16 = simd_swizzle!(lows, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let b_lo: u16x16 = simd_swizzle!(
        lows,
        [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
    );
    let a_hi: u16x16 = simd_swizzle!(
        highs,
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    );
    let b_hi: u16x16 = simd_swizzle!(
        highs,
        [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
    );

    let a_lo32: u32x16 = a_lo.cast();
    let a_hi32: u32x16 = a_hi.cast();
    let a = a_lo32 + (a_hi32 << 16);
    let b_lo32: u32x16 = b_lo.cast();
    let b_hi32: u32x16 = b_hi.cast();
    let b = b_lo32 + (b_hi32 << 16);

    (a, b)
}

#[inline]
pub fn concat_u16x16(a: u16x16, b: u16x16) -> u16x32 {
    let mut buf = [0u16; 32];
    // copy the 16 lanes from each half
    buf[..16].copy_from_slice(&a.to_array());
    buf[16..].copy_from_slice(&b.to_array());
    u16x32::from_array(buf)
}

#[inline]
pub fn pext_u32_u16x32(output_low: u16x32, output_high: u16x32, mask: u32) -> u16x32 {
    let (output_0, output_1) = widen_add_u16x32(output_low, output_high);
    let output_o2_0: u16x16 = pext_u32x16(output_0, mask).cast();
    let output_o2_1: u16x16 = pext_u32x16(output_1, mask).cast();
    concat_u16x16(output_o2_0, output_o2_1)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_subset_iterator() {
        let mask = 0b101;
        let result = SubsetIterator::new(mask);
        assert_eq!(result.collect::<Vec<_>>(), vec![0, 1, 4, 5]);
    }
    #[test]
    fn test_subset_iterator_rev() {
        let mask = 0b101;
        let result = SubsetIterator::new(mask).rev();
        assert_eq!(result.collect::<Vec<_>>(), vec![5, 4, 1, 0]);
    }

    #[test]
    fn test_subset_iterator_pext() {
        let mask = 0b101;
        let result = SubsetIterator::new(mask).map(|x| pext_u32(x, mask));
        assert_eq!(result.collect::<Vec<_>>(), vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_pext_u32() {
        let x = 0b1110;
        let mask = 0b101;
        assert_eq!(pext_u32(x, mask), 0b10);
    }

    #[test]
    fn test_pext_u32x16() {
        let x = u32x16::splat(0b1110);
        let mask = 0b101;
        assert_eq!(pext_u32x16(x, mask), u32x16::splat(0b10));
    }
}
