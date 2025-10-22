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

                pub fn to_ids() -> Vec<
                    stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId
                > {
                    vec![
                        $(
                            stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId {
                                id: format!("{}_{}", stringify!($name), stringify!($column)),
                            }
                        ),*
                    ]
                }
            }
        }

        // ---------- Static version ----------
        #[allow(dead_code)]
        impl $name<'static, ()> {
            pub const SIZE: usize = <[()]>::len(&[$(trace_columns!(@unit $column)),*]);

            pub fn to_ids() -> Vec<
                stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId
            > {
                vec![
                    $(
                        stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId {
                            id: format!("{}_{}", stringify!($name), stringify!($column)),
                        }
                    ),*
                ]
            }
        }
    };

    // helper
    (@unit $_field:ident) => { () };
}

#[macro_export]
macro_rules! combine {
    ($relations:expr, $($col:expr),+ $(,)?) => {
        {
    #[allow(clippy::tuple_array_conversions)]
    let combined: Vec<stwo::prover::backend::simd::qm31::PackedQM31> =
        itertools::izip!($($col),+)
            .map(|($(paste::paste!{ [<_elt_ $col>] }),+)| {
                let simd_values = [$(paste::paste!{ [<_elt_ $col>] }),+];
                let packed_m31_values = unsafe {
                    simd_values.map(|&v| stwo::prover::backend::simd::m31::PackedM31::from_simd_unchecked(v))
                };
                $relations.combine(&packed_m31_values)
            })
            .collect();
    combined
        }
    };
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

#[macro_export]
macro_rules! components {
    // handle foo::bar forms
    ( $( $module:ident :: $name:ident ),+ $(,)? ) => {
        $crate::components!(@impl $( ($module, $name) ),+ );
    };
    // internal implementation
    (@impl $( ($module:ident, $name:ident) ),+ $(,)? ) => {
        use $crate::relations::Relations;
        use $crate::circle_evaluation_u32x16;

        use stwo::core::fields::qm31::SecureField;
        use std::simd::u32x16;
        use num_traits::Zero;
        use stwo::{
            core::{fields::{m31::BaseField}, ColumnVec},
            prover::{backend::simd::{SimdBackend}, poly::{circle::CircleEvaluation, BitReversedOrder}, ComponentProver},
        };
        use stwo::core::pcs::TreeVec;
        use stwo::core::air::Component;
        use stwo_constraint_framework::{
            relation_tracker::{add_to_relation_entries, RelationTrackerEntry},
            TraceLocationAllocator,
        };

        paste::paste! {
            #[derive(Clone)]
            pub struct Traces {
                $( pub [<$module _ $name>]: Vec<Vec<u32x16>>, )+
            }

            pub struct ClaimedSum {
                $( pub [<$module _ $name>]: SecureField, )+
            }

            impl ClaimedSum {
                pub fn sum(&self) -> SecureField {
                    SecureField::zero() $( + self.[<$module _ $name>] )+
                }
            }

            pub fn gen_trace(
                scheduling_lookup_data: &[Vec<u32x16>],
                compression_lookup_data: &[Vec<u32x16>],
            ) -> Traces {
                Traces {
                    $(
                        [<$module _ $name>]: $module::$name::witness::gen_trace(
                            scheduling_lookup_data,
                            compression_lookup_data
                        ),
                    )+
                }
            }

            pub fn gen_interaction_trace(
                traces: &Traces,
                relations: &Relations,
            ) -> (
                ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
                ClaimedSum,
            ) {
                let mut interaction_trace = vec![];

                $(
                    let (trace, [<$module _ $name _ claimed_sum>]) =
                        $module::$name::witness::gen_interaction_trace(&traces.[<$module _ $name>], relations);
                    interaction_trace.extend(trace);
                )+
                (
                    interaction_trace,
                    ClaimedSum {
                        $( [<$module _ $name>]: [<$module _ $name _ claimed_sum>], )+
                    },
                )
            }

            impl Traces {
                pub fn len(&self) -> usize {
                    0 $( + self.[<$module _ $name>].iter().map(|v| v.len()).sum::<usize>() )+
                }

                pub fn is_empty(&self) -> bool {
                    self.len() == 0
                }
            }

            impl IntoIterator for Traces {
                type Item = CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>;
                type IntoIter = std::vec::IntoIter<Self::Item>;

                fn into_iter(self) -> Self::IntoIter {
                    let cols = vec![
                        $( self.[<$module _ $name>] ),+
                    ];

                    cols
                        .into_iter()
                        .flatten()
                        .map(|col| circle_evaluation_u32x16!(col))
                        .collect::<Vec<_>>()
                        .into_iter()
                }
            }

            pub struct Components {
                $( pub [<$module _ $name>]: $module::$name::air::Component, )+
            }

            impl Components {
                pub fn new(
                    location_allocator: &mut TraceLocationAllocator,
                    relations: Relations,
                    claimed_sum: &ClaimedSum,
                ) -> Self {
                    Self {
                        $(
                            [<$module _ $name>]: $module::$name::air::Component::new(
                                location_allocator,
                                $module::$name::air::Eval { relations: relations.clone() },
                                claimed_sum.[<$module _ $name>],
                            ),
                        )+
                    }
                }

                pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
                    vec![ $(&self.[<$module _ $name>],)+ ]
                }

                pub fn relation_entries(
                    &self,
                    trace: &TreeVec<Vec<&Vec<BaseField>>>,
                ) -> Vec<RelationTrackerEntry> {
                    itertools::chain!(
                        $( add_to_relation_entries(&self.[<$module _ $name>], trace) ),+
                    )
                    .collect()
                }

                pub fn trace_log_degree_bounds(&self) -> Vec<TreeVec<ColumnVec<u32>>> {
                    vec![
                        $( self.[<$module _ $name>].trace_log_degree_bounds(), )+
                    ]
                }
            }
        }
    };
}
