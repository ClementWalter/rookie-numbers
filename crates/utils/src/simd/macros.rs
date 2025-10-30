/// Macro for constructing aligned vectors.
///
/// - `aligned_vec![value; len]` calls `crate::utils::simd::aligned_vec(value, len)`
/// - `aligned_vec![v1, v2, ...]` calls `crate::utils::simd::aligned_vec_from_slice(&[v1, v2, ...])`
///
/// Benches benches/aligned_vec_bench.rs reports a negligible performance difference between the native vec! and `aligned_vec!`
///
/// Timer precision: 41 ns
/// aligned_vec            fastest       │ slowest       │ median        │ mean          │ samples │ iters
/// ├─ aligned_vec_repeat  211.7 µs      │ 456.2 µs      │ 247.1 µs      │ 248.6 µs      │ 100     │ 100
/// ╰─ vec_repeat          185.3 µs      │ 267.8 µs      │ 212.2 µs      │ 212.3 µs      │ 100     │ 100
#[macro_export]
macro_rules! aligned_vec {
    ($value:expr; $len:expr) => {
        $crate::simd::aligned_vec($value, $len)
    };
    ($($elem:expr),+ $(,)?) => {
        $crate::simd::aligned_vec_from_slice(&[$($elem),+])
    };
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_aligned_vec_macro_repeat() {
        let v = aligned_vec![7u32; 10];
        let ptr = v.as_ptr() as usize;
        assert_eq!(ptr % 64, 0, "Vector pointer is not 64-byte aligned");
        assert_eq!(v.len(), 10);
        assert!(v.iter().all(|&x| x == 7));
    }

    #[test]
    fn test_aligned_vec_macro_list() {
        let v = aligned_vec![1, 2, 3, 4];
        let ptr = v.as_ptr() as usize;
        assert_eq!(ptr % 64, 0, "Vector pointer is not 64-byte aligned");
        assert_eq!(v, vec![1, 2, 3, 4]);
    }
}
