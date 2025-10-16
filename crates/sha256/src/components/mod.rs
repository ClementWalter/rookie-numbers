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

    let mut trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(scheduling::witness::N_COLUMNS + compression::witness::N_COLUMNS);

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
    let mut interaction_trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(
            scheduling::witness::N_INTERACTION_COLUMNS
                + compression::witness::N_INTERACTION_COLUMNS,
        );

    let span = span!(Level::INFO, "Scheduling").entered();
    let (scheduling_interaction_trace, scheduling_claimed_sum) =
        scheduling::witness::gen_interaction_trace(&lookup_data.scheduling, relations);
    interaction_trace.extend(scheduling_interaction_trace);
    span.exit();

    let span = span!(Level::INFO, "Compression").entered();
    let (compression_interaction_trace, compression_claimed_sum) =
        compression::witness::gen_interaction_trace(&lookup_data.compression, relations);
    interaction_trace.extend(compression_interaction_trace);
    span.exit();

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
    ($name:ident,$($column:ident),* $(,)?) => {
        #[derive(Debug)]
        pub struct $name<T> {
            $(pub $column: T),*
        }

        #[allow(dead_code)]
        impl<T: Clone + Copy> $name<T> {
            pub fn to_vec(&self) -> Vec<T> {
                vec![$(self.$column),*]
            }
        }

        #[allow(dead_code)]
        impl<T> $name<T> {
            pub fn from_eval<E>(eval: &mut E) -> Self
            where
                E: stwo_constraint_framework::EvalAtRow<F = T>,
            {
                Self {
                    $($column: eval.next_trace_mask(),)*
                }
            }
        }

        #[allow(dead_code)]
        impl<T> $name<T> {
            pub fn from_ids<E>(eval: &mut E) -> Self
            where
                E: stwo_constraint_framework::EvalAtRow<F = T>,
            {
                use stwo_constraint_framework::preprocessed_columns::PreProcessedColumnId;
                Self {
                    $($column: eval.get_preprocessed_column(PreProcessedColumnId { id: format!("{}", stringify!($column)) }),)*
                }
            }
        }

        impl<T> std::iter::FromIterator<T> for $name<T> {
            fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
                let mut it = iter.into_iter();
                Self {
                    $($column: it.next().unwrap(),)*
                }
            }
        }

        paste::paste! {
            #[allow(dead_code)]
            #[repr(usize)]
            pub enum [<$name Index>] {
                $($column),*
            }
        }
    };
}

#[macro_export]
macro_rules! combine {
    ($relations:expr, $data:expr, $base_index:expr, $($col:ident),+ $(,)?) => {{
        let simd_sizes = [$(
            $data[$base_index + InteractionColumnsIndex::$col as usize].len(),
        )+];

        let simd_size = *simd_sizes.iter().max().unwrap();
        let mut combined = Vec::with_capacity(simd_size);
        for vec_row in 0..simd_size {
            let simd_values = [
                $(
                    $data[$base_index + InteractionColumnsIndex::$col as usize][vec_row],
                )+
            ];
            let packed_m31_values = unsafe {
                simd_values.map(|value| stwo::prover::backend::simd::m31::PackedM31::from_simd_unchecked(value))
            };

            let denom: stwo::prover::backend::simd::qm31::PackedQM31 =
                $relations.combine(&packed_m31_values);
            combined.push(denom);
        }
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
macro_rules! write_col {
    // --- Case 1: only denoms, numerators are just 1 ---
    ($denom: expr, $interaction_trace:expr) => {
        use num_traits::One;
        let simd_size = $denom.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = stwo::prover::backend::simd::qm31::PackedQM31::one();
            let denom = $denom[vec_row];
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    };

    // --- Case 2: with explicit numerators ---
    ($num: expr, $denom: expr, $interaction_trace:expr) => {
        use num_traits::One;
        let simd_size = $denom.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = stwo::prover::backend::simd::qm31::PackedQM31::one() * $num[vec_row];
            let denom = $denom[vec_row];
            col.write_frac(vec_row, numerator, denom);
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
        use num_traits::One;
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator_0 =
                stwo::prover::backend::simd::qm31::PackedQM31::one() * $numerator_0[vec_row];
            let numerator_1 =
                stwo::prover::backend::simd::qm31::PackedQM31::one() * $numerator_1[vec_row];
            let numerator = numerator_0 * $denom_1[vec_row] + numerator_1 * $denom_0[vec_row];
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
                BaseColumn::from_iter($column.into_iter().map(BaseField::from_u32_unchecked)),
            )),*
        ])
    };
}
