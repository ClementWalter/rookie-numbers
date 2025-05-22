use num_traits::identities::Zero;
use air::preprocessed::PreProcessedTrace;
use air::components::{ Claim, Component, Eval };
use air::Proof;
use stwo_prover::constraint_framework::TraceLocationAllocator;
use stwo_prover::core::air::ComponentProver;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::backend::BackendForChannel;
use stwo_prover::core::channel::{ MerkleChannel };
use stwo_prover::core::fields::qm31::SecureField;
use stwo_prover::core::pcs::{ CommitmentSchemeProver, PcsConfig };
use stwo_prover::core::poly::circle::{ CanonicCoset, PolyOps };
use stwo_prover::core::prover::{ prove, ProvingError };
use tracing::{ span, Level };

pub(crate) const LOG_MAX_ROWS: u32 = 26;

pub fn prove_rookie<MC: MerkleChannel>(log_size: u32) -> Result<Proof<MC::H>, ProvingError>
    where SimdBackend: BackendForChannel<MC>
{
    let _span = span!(Level::INFO, "prove_rookie").entered();

    // Setup protocol.
    let channel = &mut MC::C::default();

    let pcs_config = PcsConfig::default();
    pcs_config.mix_into(channel);

    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(
            LOG_MAX_ROWS + pcs_config.fri_config.log_blowup_factor + 2
        ).circle_domain().half_coset
    );
    let mut commitment_scheme = CommitmentSchemeProver::<SimdBackend, MC>::new(
        pcs_config,
        &twiddles
    );

    // Preprocessed traces
    let span = span!(Level::INFO, "Preprocessed trace").entered();
    let preprocessed_trace = PreProcessedTrace::default();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.gen_trace());
    tree_builder.commit(channel);
    span.exit();

    // Execution traces
    let span = span!(Level::INFO, "Execution trace").entered();
    let claim = Claim::new(log_size);
    claim.mix_into(channel);

    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(claim.write_trace().to_evals());
    tree_builder.commit(channel);
    span.exit();

    // Interaction traces
    let span = span!(Level::INFO, "Interaction trace").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(vec![]);
    tree_builder.commit(channel);
    span.exit();

    // Prove stark.
    let span = span!(Level::INFO, "Prove STARKs").entered();
    let mut tree_span_provider = TraceLocationAllocator::new_with_preproccessed_columns(
        &preprocessed_trace.ids()
    );
    let eval = Eval { claim };
    let component = Component::new(&mut tree_span_provider, eval, SecureField::zero());
    let stark_proof = prove::<SimdBackend, _>(
        &[&component as &dyn ComponentProver<SimdBackend>],
        channel,
        commitment_scheme
    )?;
    span.exit();

    Ok(Proof {
        claim,
        stark_proof,
    })
}
