use std::simd::u32x16;

#[inline(always)]
/// Generate a sequence of u32x16 values.
///
/// # Arguments
///
/// * `start` - The starting value.
/// * `len` - The length of the sequence.
///
/// # Returns
///
/// A vector of u32x16 values.
pub fn generate_simd_sequence_bulk(start: usize, len: usize) -> Vec<u32x16> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_simd_sequence_bulk() {
        let sequence = generate_simd_sequence_bulk(0, 16);
        assert_eq!(sequence.len(), 1);
        assert_eq!(
            sequence[0],
            u32x16::from_array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        );

        let sequence = generate_simd_sequence_bulk(1, 16);
        assert_eq!(sequence.len(), 1);
        assert_eq!(
            sequence[0],
            u32x16::from_array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
        );
    }
}
