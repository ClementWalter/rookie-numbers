use stwo_prover::core::{
    backend::simd::{m31::LOG_N_LANES, SimdBackend},
    fields::m31::BaseField,
    poly::{circle::CircleEvaluation, BitReversedOrder},
    ColumnVec,
};
use tracing::{span, Level};

use crate::relations::LookupData;
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

    let lookup_data = LookupData::new(log_size);

    trace.extend(scheduling_trace);
    trace.extend(compression_trace);

    (trace, lookup_data)
}
