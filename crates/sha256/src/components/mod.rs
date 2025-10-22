use std::simd::u32x16;

use stwo::{
    core::{
        air::Component,
        channel::MerkleChannel,
        fields::{m31::BaseField, qm31::SecureField},
        pcs::TreeVec,
        poly::circle::CanonicCoset,
        ColumnVec,
    },
    prover::{
        backend::{
            simd::{m31::LOG_N_LANES, SimdBackend},
            BackendForChannel, Column,
        },
        poly::{circle::CircleEvaluation, BitReversedOrder},
        CommitmentSchemeProver, ComponentProver,
    },
};
use stwo_constraint_framework::{
    relation_tracker::{add_to_relation_entries, RelationSummary, RelationTrackerEntry},
    TraceLocationAllocator,
};
use tracing::{span, Level};

use crate::relations::Relations;
pub const W_SIZE: usize = 128; // 128 u16 = 64 u32

pub mod compression;
pub mod preprocessed;
pub mod scheduling;

pub struct LookupData {
    pub scheduling: Vec<Vec<u32x16>>,
    pub compression: Vec<Vec<u32x16>>,
    pub preprocessed: preprocessed::Traces,
}

pub struct ClaimedSum {
    pub scheduling: SecureField,
    pub compression: SecureField,
    pub preprocessed: preprocessed::ClaimedSum,
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

    let span = span!(Level::INFO, "Preprocessed").entered();
    let preprocessed_trace =
        preprocessed::gen_trace(&scheduling_lookup_data, &compression_lookup_data);
    span.exit();

    let lookup_data = LookupData {
        scheduling: scheduling_lookup_data,
        compression: compression_lookup_data,
        preprocessed: preprocessed_trace.clone(),
    };

    let mut trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(
            scheduling_trace.len() + compression_trace.len() + preprocessed_trace.len(),
        );
    trace.extend(scheduling_trace);
    trace.extend(compression_trace);
    trace.extend(preprocessed_trace);

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

    let span = span!(Level::INFO, "Preprocessed").entered();
    let (preprocessed_interaction_trace, preprocessed_claimed_sum) =
        preprocessed::gen_interaction_trace(&lookup_data.preprocessed, relations);
    span.exit();

    let mut interaction_trace: Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> =
        Vec::with_capacity(
            scheduling_interaction_trace.len()
                + compression_interaction_trace.len()
                + preprocessed_interaction_trace.len(),
        );
    interaction_trace.extend(scheduling_interaction_trace);
    interaction_trace.extend(compression_interaction_trace);
    interaction_trace.extend(preprocessed_interaction_trace);
    (
        interaction_trace,
        ClaimedSum {
            scheduling: scheduling_claimed_sum,
            compression: compression_claimed_sum,
            preprocessed: preprocessed_claimed_sum,
        },
    )
}

pub struct Components {
    scheduling: scheduling::air::Component,
    compression: compression::air::Component,
    preprocessed: preprocessed::Components,
}

impl Components {
    pub fn new(
        log_size: u32,
        location_allocator: &mut TraceLocationAllocator,
        relations: &Relations,
        claimed_sum: &ClaimedSum,
    ) -> Self {
        Self {
            scheduling: scheduling::air::Component::new(
                location_allocator,
                scheduling::air::Eval {
                    log_size,
                    relations: relations.clone(),
                },
                claimed_sum.scheduling,
            ),
            compression: compression::air::Component::new(
                location_allocator,
                compression::air::Eval {
                    log_size,
                    relations: relations.clone(),
                },
                claimed_sum.compression,
            ),
            preprocessed: preprocessed::Components::new(
                location_allocator,
                relations.clone(),
                &claimed_sum.preprocessed,
            ),
        }
    }
}

impl Components {
    pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
        let mut provers: Vec<&dyn ComponentProver<SimdBackend>> =
            vec![&self.scheduling, &self.compression];
        provers.extend(self.preprocessed.provers());
        provers
    }

    pub fn track_relations<MC: MerkleChannel>(
        &self,
        commitment_scheme: &CommitmentSchemeProver<'_, SimdBackend, MC>,
    ) -> RelationSummary
    where
        SimdBackend: BackendForChannel<MC>,
    {
        let evals = commitment_scheme.trace().polys.map(|tree| {
            tree.iter()
                .map(|poly| {
                    poly.evaluate(CanonicCoset::new(poly.log_size()).circle_domain())
                        .values
                        .to_cpu()
                })
                .collect()
        });
        let evals = &evals.as_ref();
        let trace = &evals.into();

        let entries: Vec<RelationTrackerEntry> = itertools::chain!(
            add_to_relation_entries(&self.scheduling, trace),
            add_to_relation_entries(&self.compression, trace),
        )
        .chain(self.preprocessed.relation_entries(trace))
        .collect();

        RelationSummary::summarize_relations(&entries).cleaned()
    }

    pub fn trace_log_degree_bounds(&self) -> Vec<TreeVec<ColumnVec<u32>>> {
        let mut log_degree_bounds: Vec<TreeVec<ColumnVec<u32>>> = Vec::new();
        log_degree_bounds.push(self.scheduling.trace_log_degree_bounds());
        log_degree_bounds.push(self.compression.trace_log_degree_bounds());
        log_degree_bounds.extend(self.preprocessed.trace_log_degree_bounds());
        log_degree_bounds
    }
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
