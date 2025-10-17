use std::simd::u32x16;

use stwo::{
    core::{
        fields::{m31::BaseField, qm31::SecureField},
        ColumnVec,
    },
    prover::{
        backend::simd::{m31::LOG_N_LANES, SimdBackend},
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use tracing::{span, Level};

use crate::relations::{LookupData, Relations};
pub const W_SIZE: usize = 128; // 128 u16 = 64 u32

pub mod compression;
pub mod preprocessed;
pub mod scheduling;

pub struct ClaimedSum {
    pub scheduling: SecureField,
    pub compression: SecureField,
}

pub fn gen_trace(
    log_size: u32,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    LookupData,
) {
    assert!(log_size >= LOG_N_LANES);

    let span = span!(Level::INFO, "Scheduling").entered();
    let (scheduling_trace, scheduling_lookup_data) = scheduling::witness::gen_trace(log_size);
    span.exit();

    let span = span!(Level::INFO, "Compression").entered();
    let (compression_trace, compression_lookup_data) =
        compression::witness::gen_trace(&scheduling_trace);
    span.exit();

    let lookup_data = LookupData {
        scheduling: scheduling_lookup_data,
        compression: compression_lookup_data,
    };

    let mut trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(scheduling_trace.len() + compression_trace.len());
    trace.extend(scheduling_trace);
    trace.extend(compression_trace);

    (trace, lookup_data)
}

pub fn gen_interaction_trace(
    lookup_data: LookupData,
    relations: &Relations,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    ClaimedSum,
) {
    let span = span!(Level::INFO, "Scheduling").entered();
    let (scheduling_interaction_trace, scheduling_claimed_sum) =
        scheduling::witness::gen_interaction_trace(&lookup_data.scheduling, relations);
    span.exit();

    let span = span!(Level::INFO, "Compression").entered();
    let (compression_interaction_trace, compression_claimed_sum) =
        compression::witness::gen_interaction_trace(&lookup_data.compression, relations);
    span.exit();

    let mut interaction_trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(
            scheduling_interaction_trace.len() + compression_interaction_trace.len(),
        );
    interaction_trace.extend(scheduling_interaction_trace);
    interaction_trace.extend(compression_interaction_trace);

    (
        interaction_trace,
        ClaimedSum {
            scheduling: scheduling_claimed_sum,
            compression: compression_claimed_sum,
        },
    )
}
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
                        $($column: eval.get_preprocessed_column(stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId { id: format!("{}", stringify!($column)) }),)*
                    }
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
                            id: format!("{}", stringify!($column)),
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
    ($relations:expr, $($col:expr),+ $(,)?) => {{
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
    }};
}

#[inline(always)]
pub fn combine_w(
    relations: &Relations,
    data: &[Vec<u32x16>],
) -> Vec<stwo::prover::backend::simd::qm31::PackedQM31> {
    use stwo_constraint_framework::Relation;

    use crate::components::W_SIZE;

    let simd_size = data[0].len();
    let mut combined = Vec::with_capacity(simd_size);
    for vec_row in 0..simd_size {
        unsafe {
            let values: [stwo::prover::backend::simd::m31::PackedM31; W_SIZE] = (0..W_SIZE)
                .map(|i| {
                    stwo::prover::backend::simd::m31::PackedM31::from_simd_unchecked(
                        data[i][vec_row],
                    )
                })
                .collect::<Vec<_>>()
                .try_into()
                .unwrap();
            let denom: stwo::prover::backend::simd::qm31::PackedQM31 = relations.w.combine(&values);
            combined.push(denom);
        }
    }
    combined
}

#[macro_export]
macro_rules! emit_col {
    ($denom: expr, $interaction_trace:expr) => {
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
    ($denom: expr, $interaction_trace:expr) => {
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
macro_rules! write_pair {
    // --- Case 1: only denoms, numerators are just 1 ---
    ($denom_0:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = $denom_0[vec_row] + $denom_1[vec_row];
            let denom = $denom_0[vec_row] * $denom_1[vec_row];
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    }};

    // --- Case 2: with explicit numerators ---
    ($numerator_0:expr, $denom_0:expr, $numerator_1:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = $numerator_0[vec_row] * $denom_1[vec_row]
                + $numerator_1[vec_row] * $denom_0[vec_row];
            let denom = $denom_0[vec_row] * $denom_1[vec_row];
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! consume_pair {
    ($denom_0:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = $denom_0[vec_row] + $denom_1[vec_row];
            let denom = $denom_0[vec_row] * $denom_1[vec_row];
            col.write_frac(vec_row, -numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! emit_pair {
    ($denom_0:expr, $denom_1:expr, $interaction_trace:expr) => {{
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = $denom_0[vec_row] + $denom_1[vec_row];
            let denom = $denom_0[vec_row] * $denom_1[vec_row];
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    }};
}

#[macro_export]
macro_rules! add_to_relation {
    ($eval:expr, $relation:expr, $numerator:expr, $($col:expr),+ $(,)?) => {{
        $eval.add_to_relation(stwo_constraint_framework::RelationEntry::new(
            &$relation,
            $numerator.clone(),
            &[$($col.clone()),*],
        ))
    }};
}

#[macro_export]
macro_rules! column_vec {
    ($($column:ident),*) => {
        ColumnVec::from(vec![
            $(CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
                CanonicCoset::new($column.len().ilog2()).circle_domain(),
                BaseColumn::from_iter($column.iter().map(|v| BaseField::from_u32_unchecked(*v))),
            )),*
        ])
    };
}

#[macro_export]
macro_rules! simd_vec {
    ($($column:ident),*) => {
        vec![
            $(
                $column
                .chunks(16)
                .into_iter()
                .map(u32x16::from_slice)
                .collect::<Vec<u32x16>>()
        ),*
        ]
    };
}
