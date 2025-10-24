use std::{
    alloc::{alloc_zeroed, handle_alloc_error, Layout},
    ptr::write,
};

pub mod macros;

/// Creates a Vec<T> with 64-byte alignment, filled with clones of `value`.
pub fn aligned_vec<T: Clone>(value: T, len: usize) -> Vec<T> {
    let elem_size = std::mem::size_of::<T>();
    let align = 64.max(std::mem::align_of::<T>());
    let layout = Layout::from_size_align(
        len.checked_mul(elem_size)
            .expect("Overflow in allocation size"),
        align,
    )
    .unwrap();
    unsafe {
        let ptr = alloc_zeroed(layout) as *mut T;
        if ptr.is_null() {
            handle_alloc_error(layout);
        }
        for i in 0..len {
            write(ptr.add(i), value.clone());
        }
        Vec::from_raw_parts(ptr, len, len)
    }
}

/// Creates a Vec<T> with 64-byte alignment, filled with clones of elements from the slice.
pub fn aligned_vec_from_slice<T: Clone>(elems: &[T]) -> Vec<T> {
    let len = elems.len();
    let elem_size = std::mem::size_of::<T>();
    let align = 64.max(std::mem::align_of::<T>());
    let layout = Layout::from_size_align(
        len.checked_mul(elem_size)
            .expect("Overflow in allocation size"),
        align,
    )
    .unwrap();
    unsafe {
        let ptr = alloc_zeroed(layout) as *mut T;
        if ptr.is_null() {
            handle_alloc_error(layout);
        }
        for (i, v) in elems.iter().enumerate() {
            write(ptr.add(i), v.clone());
        }
        Vec::from_raw_parts(ptr, len, len)
    }
}

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

/// Pivot a vector of `[u32x16; N]` rows into `N` column vectors.
///
///
/// # Arguments
///
/// * `rows` - A vector of `[u32x16; N]` rows.
///
/// # Returns
///
/// A vector of `N` column vectors.
#[inline(always)]
pub fn pivot<const N: usize>(rows: Vec<[u32x16; N]>) -> [Vec<u32x16>; N] {
    // Preallocate N columns, each with the same capacity as number of rows.
    let mut columns: [Vec<u32x16>; N] = core::array::from_fn(|_| Vec::with_capacity(rows.len()));

    for row in rows {
        for (col, value) in columns.iter_mut().zip(row) {
            col.push(value);
        }
    }

    columns
}

/// Flatten a vector of `u32x16` values into a vector of `u32` values.
///
/// # Arguments
///
/// * `column` - A vector of `u32x16` values.
///
/// # Returns
///
/// A vector of `u32` values.
#[inline(always)]
pub const fn flatten_simd(column: &[u32x16]) -> &[u32] {
    unsafe { std::slice::from_raw_parts(column.as_ptr() as *const u32, column.len() * 16) }
}

/// Chunk a slice of `u32` values into a slice of `u32x16` values. This is the inverse of `flatten_simd`.
///
/// # Panics
///
/// Panics if the length of the column is not a multiple of 16 or the slice is not 64-byte aligned.
///
/// # Arguments
///
/// * `column` - A slice of `u32` values.
///
/// # Returns
///
/// A slice of `u32x16` values.
#[inline(always)]
pub fn into_simd(column: &[u32]) -> &[u32x16] {
    let align = std::mem::align_of::<u32x16>();
    let ptr = column.as_ptr() as usize;

    if !column.len().is_multiple_of(16) {
        panic!(
            "into_simd: length {} is not a multiple of 16 (required for u32x16)",
            column.len()
        );
    }

    if !ptr.is_multiple_of(align) {
        panic!("into_simd: pointer {ptr:#x} is not {align}-byte aligned (required for u32x16)");
    }

    unsafe { std::slice::from_raw_parts(column.as_ptr() as *const u32x16, column.len() / 16) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aligned_vec_u32_repeat() {
        let v = aligned_vec(1u32, 32);
        let ptr = v.as_ptr() as usize;
        assert_eq!(ptr % 64, 0, "Vector pointer is not 64-byte aligned");
        assert_eq!(v.len(), 32);
        assert!(v.iter().all(|&x| x == 1));
    }

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

    #[test]
    fn test_pivot() {
        let rows = vec![
            [u32x16::splat(0), u32x16::splat(1), u32x16::splat(2)],
            [u32x16::splat(3), u32x16::splat(4), u32x16::splat(5)],
        ];
        let columns = pivot::<3>(rows);
        assert_eq!(
            columns,
            [
                vec![u32x16::splat(0), u32x16::splat(3)],
                vec![u32x16::splat(1), u32x16::splat(4)],
                vec![u32x16::splat(2), u32x16::splat(5)]
            ]
        );
    }

    #[test]
    fn test_flatten_simd() {
        let column = vec![u32x16::splat(0), u32x16::splat(1)];
        let flattened = flatten_simd(&column);
        assert_eq!(
            flattened,
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1
            ]
        );
    }

    #[test]
    fn test_chunk_simd_column() {
        let column =
            aligned_vec_from_slice(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let chunked = into_simd(&column);
        assert_eq!(chunked, [u32x16::from_slice(&column)].to_vec());
    }

    #[test]
    fn test_flatten_then_chunk() {
        let column = vec![u32x16::splat(0), u32x16::splat(1)];
        let flattened = flatten_simd(&column);
        let chunked = into_simd(flattened);
        assert_eq!(chunked, column);
    }

    #[test]
    fn test_chunk_then_flatten() {
        let column =
            aligned_vec_from_slice(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let chunked = into_simd(&column);
        let flattened = flatten_simd(chunked);
        assert_eq!(flattened, &column);
    }
}
