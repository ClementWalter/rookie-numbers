use num_traits::Zero;
use stwo::{
    core::{
        channel::MerkleChannel, fields::qm31::SecureField, pcs::PcsConfig,
        poly::circle::CanonicCoset, proof::StarkProof,
    },
    prover::{
        backend::{simd::SimdBackend, BackendForChannel},
        poly::circle::PolyOps,
        prove, CommitmentSchemeProver, ProvingError as StwoProvingError,
    },
};
use stwo_constraint_framework::TraceLocationAllocator;
use tracing::info;
use utils::{circle_evaluation_u32x16, simd::generate_simd_sequence_bulk};

use crate::air::components::{DummyComponent, DummyEval};

pub fn prove_dummy<MC: MerkleChannel, const N: usize>(
    log_size: u32,
) -> Result<StarkProof<MC::H>, StwoProvingError>
where
    SimdBackend: BackendForChannel<MC>,
{
    // Setup protocol.
    let channel = &mut MC::C::default();

    let pcs_config = PcsConfig::default();
    pcs_config.mix_into(channel);

    info!("twiddles");
    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(log_size + pcs_config.fri_config.log_blowup_factor + 2)
            .circle_domain()
            .half_coset,
    );
    let mut commitment_scheme =
        CommitmentSchemeProver::<SimdBackend, MC>::new(pcs_config, &twiddles);

    // Preprocessed trace
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(vec![]);
    tree_builder.commit(channel);

    // Generate trace
    info!("trace");
    let col = generate_simd_sequence_bulk(0, 1 << log_size);
    let col = circle_evaluation_u32x16!(col);
    let trace = vec![col; N];

    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(trace);
    tree_builder.commit(channel);

    // Prove stark.
    info!("prove stark");

    let component = DummyComponent::new(
        &mut TraceLocationAllocator::default(),
        DummyEval {
            log_size,
            n_cols: N,
        },
        SecureField::zero(),
    );
    info!("Dummy component info:\n{}", component);

    prove::<SimdBackend, _>(&[&component], channel, commitment_scheme)
}
