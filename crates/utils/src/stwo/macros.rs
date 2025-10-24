#[macro_export]
macro_rules! trace_columns {
    ($name:ident, $($column:ident),* $(,)?) => {
        // ---------- Borrow version ----------
        #[derive(Debug, Clone, Copy)]
        pub struct $name<'a, T> {
            $(pub $column: &'a T),*
        }

        #[allow(dead_code)]
        impl<'a, T> $name<'a, T> {
            // ---------- Immutable "view" ----------
            #[inline(always)]
            pub fn from_slice(slice: &'a [T]) -> Self {
                assert!(
                    slice.len() == <[()]>::len(&[$(trace_columns!(@unit $column)),*]),
                    "slice length mismatch for {}",
                    stringify!($name)
                );
                let mut it = slice.iter();
                Self {
                    $(
                        $column: it.next().expect("slice too short"),
                    )*
                }
            }

            #[inline(always)]
            pub fn iter(&self) -> impl Iterator<Item = &'a T> {
                // Builds a fixed array of &T; no copies of T occur.
                [$(
                    self.$column
                ),*].into_iter()
            }

        }

        #[allow(dead_code)]
        impl $name<'static, ()> {
            pub const SIZE: usize = <[()]>::len(&[$(trace_columns!(@unit $column)),*]);

            pub fn to_ids(suffix: Option<usize>) -> Vec<
                stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId
            > {
                vec![
                    $(
                        stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId {
                            id: match suffix {
                                Some(suffix) => format!("{}_{}_{}", stringify!($name), stringify!($column), suffix.to_string()),
                                None => format!("{}_{}", stringify!($name), stringify!($column)),
                            },
                        }
                    ),*
                ]
            }
        }

        // ---------- Owned version ----------
        paste::paste! {
            #[derive(Debug, Clone)]
            #[allow(dead_code)]
            pub struct [<$name Owned>]<T> {
                $(pub $column: T),*
            }

            #[allow(dead_code)]
            impl<T> [<$name Owned>]<T> {
                #[inline(always)]
                pub fn from_eval<E>(eval: &mut E) -> Self
                where
                    E: stwo_constraint_framework::EvalAtRow<F = T>,
                {
                    Self {
                        $(
                            $column: eval.next_trace_mask(),
                        )*
                    }
                }

                pub fn from_ids<E>(eval: &mut E) -> Self
                where
                    E: stwo_constraint_framework::EvalAtRow<F = T>,
                {
                    Self {
                        $($column: eval.get_preprocessed_column(stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId { id: format!("{}_{}", stringify!($name), stringify!($column)) }),)*
                    }
                }
            }

            #[allow(dead_code)]
            impl [<$name Owned>]<()> {
                pub const SIZE: usize = <[()]>::len(&[$(trace_columns!(@unit $column)),*]);
            }
        }

    };

    // helper
    (@unit $_field:ident) => { () };
}

#[macro_export]
macro_rules! combine {
    ($relations:expr, $cols:expr $(,)?) => {{
        let cols = $cols;
        let simd_size = cols[0].len();
        let n_cols = cols.len();

        let mut combined: Vec<stwo::prover::backend::simd::qm31::PackedQM31> =
            Vec::with_capacity(simd_size);

        // Create an iterator over all columns simultaneously
        let mut col_iters: Vec<_> = cols.iter().map(|c| c.iter()).collect();

        for _ in 0..simd_size {
            // Collect one row worth of values by pulling one from each iterator
            let mut packed_m31_values = Vec::with_capacity(n_cols);
            for it in &mut col_iters {
                let v = *it.next().unwrap();
                packed_m31_values.push(unsafe {
                    stwo::prover::backend::simd::m31::PackedM31::from_simd_unchecked(v)
                });
            }
            combined.push($relations.combine(&packed_m31_values));
        }
        combined
    }};
}

#[macro_export]
macro_rules! emit_col {
    ($denom:expr, $interaction_trace:expr) => {
        use num_traits::One;
        let mut col = $interaction_trace.new_col();
        let one = stwo::prover::backend::simd::qm31::PackedQM31::one();
        for (vec_row, &d) in $denom.iter().enumerate() {
            col.write_frac(vec_row, one, d);
        }
        col.finalize_col();
    };
}

#[macro_export]
macro_rules! consume_col {
    ($denom:expr, $interaction_trace:expr) => {
        use num_traits::One;
        let mut col = $interaction_trace.new_col();
        let minus_one = -stwo::prover::backend::simd::qm31::PackedQM31::one();
        for (vec_row, &d) in $denom.iter().enumerate() {
            col.write_frac(vec_row, minus_one, d);
        }
        col.finalize_col();
    };
}

#[macro_export]
macro_rules! write_col {
    ($numerator:expr, $denom:expr, $interaction_trace:expr) => {
        let mut col = $interaction_trace.new_col();
        for (vec_row, (n, d)) in itertools::izip!($numerator, $denom).enumerate() {
            col.write_frac(vec_row, n, d);
        }
        col.finalize_col();
    };
}

#[macro_export]
macro_rules! write_pair {
    (
        $numerator_0:expr,
        $denom_0:expr,
        $numerator_1:expr,
        $denom_1:expr,
        $interaction_trace:expr
    ) => {{
        let mut col = $interaction_trace.new_col();
        for (vec_row, (n_0, d_0, n_1, d_1)) in
            itertools::izip!($numerator_0, $denom_0, $numerator_1, $denom_1).enumerate()
        {
            let numerator = n_0 * d_1 + n_1 * d_0;
            let denom = d_0 * d_1;
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! consume_pair {
    // Variant that takes a list of columns to consume in pairs
    ($interaction_trace:expr; $($col:expr),+ $(,)?) => {{
        let secure_columns = vec![$($col),+];
        for [pair0, pair1] in secure_columns.into_iter().array_chunks::<2>() {
            let mut col = $interaction_trace.new_col();
            for (vec_row, (d_0, d_1)) in itertools::izip!(pair0.iter(), pair1.iter()).enumerate() {
                let numerator = *d_0 + *d_1;
                let denom = *d_0 * *d_1;
                col.write_frac(vec_row, -numerator, denom);
            }
            col.finalize_col();
        }
    }};

    // Variant that takes two columns to write in pairs
    ($denom_0:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let mut col = $interaction_trace.new_col();
        for (vec_row, (d_0, d_1)) in itertools::izip!($denom_0, $denom_1).enumerate() {
            let numerator = d_0 + d_1;
            let denom = d_0 * d_1;
            col.write_frac(vec_row, -numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! emit_pair {
    ($denom_0:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let mut col = $interaction_trace.new_col();
        for (vec_row, (d_0, d_1)) in itertools::izip!($denom_0, $denom_1).enumerate() {
            let numerator = d_0 + d_1;
            let denom = d_0 * d_1;
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! add_to_relation {
    ($eval:expr, $relation:expr, $numerator:expr, $($col:expr),+ $(,)?) => {
        {
        $eval.add_to_relation(stwo_constraint_framework::RelationEntry::new(
            &$relation,
            $numerator.clone(),
            &[$($col.clone()),*],
        ))
        }
    };
}

#[macro_export]
macro_rules! circle_evaluation_u32 {
    ($column:expr) => {
        stwo::prover::poly::circle::CircleEvaluation::<
            stwo::prover::backend::simd::SimdBackend,
            stwo::core::fields::m31::BaseField,
            stwo::prover::poly::BitReversedOrder,
        >::new(
            stwo::core::poly::circle::CanonicCoset::new($column.len().ilog2()).circle_domain(),
            stwo::prover::backend::simd::column::BaseColumn::from_iter(
                $column
                    .iter()
                    .map(|v| stwo::core::fields::m31::BaseField::from_u32_unchecked(*v)),
            ),
        )
    };
}

#[macro_export]
macro_rules! circle_evaluation_u32x16 {
    ($column:expr) => {
        stwo::prover::poly::circle::CircleEvaluation::<
            stwo::prover::backend::simd::SimdBackend,
            stwo::core::fields::m31::BaseField,
            stwo::prover::poly::BitReversedOrder,
        >::new(
            stwo::core::poly::circle::CanonicCoset::new(
                $column.len().ilog2() + stwo::prover::backend::simd::m31::LOG_N_LANES,
            )
            .circle_domain(),
            stwo::prover::backend::simd::column::BaseColumn::from_simd(
                $column
                    .iter()
                    .map(|v| unsafe {
                        stwo::prover::backend::simd::m31::PackedM31::from_simd_unchecked(*v)
                    })
                    .collect::<Vec<stwo::prover::backend::simd::m31::PackedM31>>(),
            ),
        )
    };
}

#[macro_export]
macro_rules! column_vec_u32 {
    ($($column:expr),*) => {
        ColumnVec::from(vec![
            $(circle_evaluation_u32!($column)),*
        ])
    };
}

#[macro_export]
macro_rules! column_vec_u32x16 {
    ($($column:expr),*) => {
        ColumnVec::from(vec![$(circle_evaluation_u32x16!($column)),*])
    };
}

#[macro_export]
macro_rules! simd_vec {
    ($($column:expr),*) => {
        vec![
            $(
                $column
                .chunks(16)
                .map(u32x16::from_slice)
                .collect::<Vec<u32x16>>()
        ),*
        ]
    };
}
