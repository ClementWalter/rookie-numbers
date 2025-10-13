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

pub fn gen_trace(
    log_size: u32,
) -> (
    ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
    LookupData,
) {
    let _span = span!(Level::INFO, "Generation").entered();
    assert!(log_size >= LOG_N_LANES);

    let mut trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(scheduling::N_COLUMNS + compression::N_COLUMNS);
    let (scheduling_trace, scheduling_lookup_data) = scheduling::gen_trace(log_size);
    let (compression_trace, compression_lookup_data) = compression::gen_trace(&scheduling_trace);

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
    SecureField,
) {
    let _span = span!(Level::INFO, "Generate interaction trace").entered();

    let mut interaction_trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(scheduling::N_INTERACTION_COLUMNS + compression::N_INTERACTION_COLUMNS);
    let mut claimed_sum = SecureField::from_u32_unchecked(0, 0, 0, 0);

    let (scheduling_interaction_trace, scheduling_claimed_sum) =
        scheduling::gen_interaction_trace(&lookup_data.scheduling, relations);
    interaction_trace.extend(scheduling_interaction_trace);
    claimed_sum += scheduling_claimed_sum;

    // // interaction_trace.extend(compression_interaction_trace);

    (interaction_trace, claimed_sum)
}

#[macro_export]
macro_rules! round_columns {
    ($name:ident,$($column:ident),*) => {
        pub struct $name<T> {
            $(pub $column: T),*
        }

        impl<T: Clone + Copy> $name<T> {
            pub fn to_vec(&self) -> Vec<T> {
                vec![$(self.$column),*]
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
    ($relations:expr, $data:expr, $base_index:expr, $simd_row:expr, $($col:ident),+ $(,)?) => {{
        unsafe {
            let denom: stwo_prover::core::backend::simd::qm31::PackedQM31 =
                $relations.combine(
                    [
                        $(
                            stwo_prover::core::backend::simd::m31::PackedM31::from_simd_unchecked(
                                $data[$base_index + InteractionColumnsIndex::$col as usize][$simd_row],
                            ),
                        )+
                    ]
                    .as_slice(),
                );
            denom
        }
    }};
}
