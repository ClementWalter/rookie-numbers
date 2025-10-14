use std::simd::u32x16;

use stwo_prover::core::{
    backend::simd::{m31::LOG_N_LANES, SimdBackend},
    fields::{m31::BaseField, qm31::SecureField},
    poly::{circle::CircleEvaluation, BitReversedOrder},
    ColumnVec,
};
use tracing::{span, Level};

use crate::relations::{LookupData, Relations};
pub const W_SIZE: usize = 128; // 128 u16 = 64 u32

pub mod compression;
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
    let _span = span!(Level::INFO, "Generation").entered();
    assert!(log_size >= LOG_N_LANES);

    let mut trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(scheduling::witness::N_COLUMNS + compression::witness::N_COLUMNS);
    let (scheduling_trace, scheduling_lookup_data) = scheduling::witness::gen_trace(log_size);
    let (compression_trace, compression_lookup_data) =
        compression::witness::gen_trace(&scheduling_trace);

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
    let _span = span!(Level::INFO, "Generate interaction trace").entered();

    let mut interaction_trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(
            scheduling::witness::N_INTERACTION_COLUMNS
                + compression::witness::N_INTERACTION_COLUMNS,
        );

    let (scheduling_interaction_trace, scheduling_claimed_sum) =
        scheduling::witness::gen_interaction_trace(&lookup_data.scheduling, relations);
    interaction_trace.extend(scheduling_interaction_trace);

    let (compression_interaction_trace, compression_claimed_sum) =
        compression::witness::gen_interaction_trace(&lookup_data.compression, relations);
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
macro_rules! round_columns {
    ($name:ident,$($column:ident),*) => {
        #[derive(Debug)]
        pub struct $name<T> {
            $(pub $column: T),*
        }

        impl<T: Clone + Copy> $name<T> {
            pub fn to_vec(&self) -> Vec<T> {
                vec![$(self.$column),*]
            }
        }

        #[allow(dead_code)]
        impl<T> $name<T> {
            pub fn from_eval<E>(eval: &mut E) -> Self
            where
                E: stwo_prover::constraint_framework::EvalAtRow<F = T>,
            {
                Self {
                    $($column: eval.next_trace_mask(),)*
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
        let simd_size = $data[0].len();
        let mut combined = Vec::with_capacity(simd_size);
        for vec_row in 0..simd_size {
            unsafe {
                let denom: stwo_prover::core::backend::simd::qm31::PackedQM31 =
                    $relations.combine(
                        [
                            $(
                                stwo_prover::core::backend::simd::m31::PackedM31::from_simd_unchecked(
                                    $data[$base_index + InteractionColumnsIndex::$col as usize][vec_row],
                                ),
                            )+
                        ]
                        .as_slice(),
                    );
                combined.push(denom);
            }
        }
        combined
    }};
}

#[inline(always)]
pub fn combine_w(
    relations: &Relations,
    data: &[Vec<u32x16>],
) -> Vec<stwo_prover::core::backend::simd::qm31::PackedQM31> {
    use crate::components::W_SIZE;
    use stwo_prover::constraint_framework::Relation;

    let simd_size = data[0].len();
    let mut combined = Vec::with_capacity(simd_size);
    for vec_row in 0..simd_size {
        unsafe {
            let values: [stwo_prover::core::backend::simd::m31::PackedM31; W_SIZE] = (0..W_SIZE)
                .map(|i| {
                    stwo_prover::core::backend::simd::m31::PackedM31::from_simd_unchecked(
                        data[i][vec_row],
                    )
                })
                .collect::<Vec<_>>()
                .try_into()
                .unwrap();
            let denom: stwo_prover::core::backend::simd::qm31::PackedQM31 =
                relations.w.combine(&values);
            combined.push(denom);
        }
    }
    combined
}

#[macro_export]
macro_rules! write_col {
    ($denom_0: expr, $denom_1: expr, $interaction_trace:expr) => {
        let simd_size = $denom_0.len();
        let mut col = $interaction_trace.new_col();
        for vec_row in 0..simd_size {
            let numerator = $denom_0[vec_row] + $denom_1[vec_row];
            let denom = $denom_0[vec_row] * $denom_1[vec_row];
            col.write_frac(vec_row, numerator, denom);
        }
        col.finalize_col();
    };
}

#[macro_export]
macro_rules! add_to_relation {
    ($eval:expr, $relation:expr, $numerator:expr, $($col:expr),+ $(,)?) => {{
        $eval.add_to_relation(stwo_prover::constraint_framework::RelationEntry::new(
            &$relation,
            $numerator.clone(),
            &[$($col.clone()),*],
        ))
    }};
}
