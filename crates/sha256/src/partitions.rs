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
    current: u32,
    mask: u32,
    done: bool, // Track if we've yielded 0
}

impl SubsetIterator {
    pub const fn new(mask: u32) -> Self {
        Self {
            current: mask,
            mask,
            done: false,
        }
    }
}

impl Iterator for SubsetIterator {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let current = self.current;
        if current == 0 {
            self.done = true; // Mark done after yielding 0
        } else {
            self.current = (current - 1) & self.mask;
        }
        Some(current)
    }
}

#[cfg(test)]
mod test {
    use super::SubsetIterator;

    #[test]
    fn test_subset_iterator() {
        let mask = 0b101;
        let result = SubsetIterator::new(mask);
        assert_eq!(result.collect::<Vec<_>>(), vec![5, 4, 1, 0]);
    }
}
