use air::components::{Claim, Components};
use air::preprocessed::PreProcessedTrace;
use air::Proof;
use std::time::Instant;
use stwo_prover::constraint_framework::TraceLocationAllocator;

use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::backend::BackendForChannel;
use stwo_prover::core::channel::MerkleChannel;
use stwo_prover::core::pcs::{CommitmentSchemeProver, CommitmentSchemeVerifier, PcsConfig};
use stwo_prover::core::poly::circle::{CanonicCoset, PolyOps};
use stwo_prover::core::prover::{prove, verify, ProvingError, VerificationError};
use tracing::{info, span, Level};

pub fn prove_rookie<MC: MerkleChannel, const N: usize>(
    log_size: u32,
) -> Result<Proof<N, MC::H>, ProvingError>
where
    SimdBackend: BackendForChannel<MC>,
{
    let _span = span!(Level::INFO, "prove_rookie").entered();

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

    // Preprocessed traces
    info!("preprocessed trace");
    let preprocessed_trace = PreProcessedTrace::default();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.gen_trace());
    tree_builder.commit(channel);

    // Execution traces
    info!("execution trace");
    let claim = Claim::new(log_size);
    claim.mix_into(channel);

    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(claim.write_trace());
    tree_builder.commit(channel);

    // Prove stark.
    info!("prove stark");
    let mut tree_span_provider =
        TraceLocationAllocator::new_with_preproccessed_columns(&preprocessed_trace.ids());
    let components = Components::new(&mut tree_span_provider, &claim);
    let proving_start = Instant::now();
    let stark_proof = prove::<SimdBackend, _>(&components.provers(), channel, commitment_scheme)?;
    let proving_duration = proving_start.elapsed();
    let proving_mhz = ((1 << log_size) as f64) / proving_duration.as_secs_f64() / 1_000_000.0;
    info!("Trace size: {:?}", 1 << log_size);
    info!("Proving time: {:?}", proving_duration);
    info!("Proving speed: {:.2} MHz", proving_mhz);

    Ok(Proof { claim, stark_proof })
}

pub fn verify_rookie<MC: MerkleChannel, const N: usize>(
    proof: Proof<N, MC::H>,
) -> Result<(), VerificationError>
where
    SimdBackend: BackendForChannel<MC>,
{
    let _span = span!(Level::INFO, "verify_rookie").entered();

    // Setup protocol.
    let channel = &mut MC::C::default();

    let pcs_config = PcsConfig::default();
    pcs_config.mix_into(channel);

    let commitment_scheme_verifier = &mut CommitmentSchemeVerifier::<MC>::new(pcs_config);

    // Preprocessed trace.
    info!("preprocessed trace");
    let preprocessed_trace = PreProcessedTrace::default();
    commitment_scheme_verifier.commit(
        proof.stark_proof.commitments[0],
        &proof.claim.log_sizes()[0],
        channel,
    );

    // Execution traces
    info!("execution trace");
    proof.claim.mix_into(channel);
    commitment_scheme_verifier.commit(
        proof.stark_proof.commitments[1],
        &proof.claim.log_sizes()[1],
        channel,
    );

    // Verify stark.
    info!("verify stark");
    let mut tree_span_provider =
        TraceLocationAllocator::new_with_preproccessed_columns(&preprocessed_trace.ids());
    let components = Components::new(&mut tree_span_provider, &proof.claim);
    verify(
        &components.verifiers(),
        channel,
        commitment_scheme_verifier,
        proof.stark_proof,
    )
}
