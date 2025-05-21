use air::preprocessed::PreProcessedTrace;
use air::interactions::{ InteractionElements, INTERACTION_POW_BITS };
use air::{ Components, Proof };
use stwo_prover::constraint_framework::TraceLocationAllocator;
use stwo_prover::core::backend::simd::SimdBackend;
use stwo_prover::core::backend::BackendForChannel;
use stwo_prover::core::channel::{ Channel, MerkleChannel };
use stwo_prover::core::pcs::{ CommitmentSchemeProver, PcsConfig };
use stwo_prover::core::poly::circle::{ CanonicCoset, PolyOps };
use stwo_prover::core::proof_of_work::GrindOps;
use stwo_prover::core::prover::{ prove, ProvingError };
use tracing::{ span, Level };

use air::components::ClaimGenerator;

pub(crate) const LOG_MAX_ROWS: u32 = 26;

pub fn prove_rookie<MC: MerkleChannel>(pcs_config: PcsConfig) -> Result<Proof<MC::H>, ProvingError>
    where SimdBackend: BackendForChannel<MC>
{
    let _span = span!(Level::INFO, "prove_rookie").entered();

    let twiddles = SimdBackend::precompute_twiddles(
        CanonicCoset::new(
            LOG_MAX_ROWS + pcs_config.fri_config.log_blowup_factor + 2
        ).circle_domain().half_coset
    );

    // Setup protocol.
    let channel = &mut MC::C::default();
    pcs_config.mix_into(channel);
    let mut commitment_scheme = CommitmentSchemeProver::<SimdBackend, MC>::new(
        pcs_config,
        &twiddles
    );

    // Preprocessed trace.
    let preprocessed_trace = PreProcessedTrace::default();
    let mut tree_builder = commitment_scheme.tree_builder();
    tree_builder.extend_evals(preprocessed_trace.gen_trace());
    tree_builder.commit(channel);

    // Run Trace
    let claim_generator = ClaimGenerator::new(1_000);

    let mut tree_builder = commitment_scheme.tree_builder();
    let span = span!(Level::INFO, "Base trace").entered();
    let (claim, interaction_generator) = claim_generator.write_trace(&mut tree_builder);
    span.exit();

    claim.mix_into(channel);
    tree_builder.commit(channel);

    // Draw interaction elements.
    let interaction_pow = SimdBackend::grind(channel, INTERACTION_POW_BITS);
    channel.mix_u64(interaction_pow);
    let interaction_elements = InteractionElements::draw(channel);

    // Interaction trace.
    let span = span!(Level::INFO, "Interaction trace").entered();
    let mut tree_builder = commitment_scheme.tree_builder();
    let interaction_claim = interaction_generator.write_interaction_trace(
        &mut tree_builder,
        &interaction_elements
    );
    span.exit();

    interaction_claim.mix_into(channel);
    tree_builder.commit(channel);

    // Component provers.
    let mut tree_span_provider = TraceLocationAllocator::new_with_preproccessed_columns(
        &preprocessed_trace.ids()
    );

    let component_builder = Components::new(
        &mut tree_span_provider,
        &claim,
        &interaction_elements,
        &interaction_claim
    );

    let components = component_builder.provers();

    // Prove stark.
    let span = span!(Level::INFO, "Prove STARKs").entered();
    let proof = prove::<SimdBackend, _>(&components, channel, commitment_scheme)?;
    span.exit();

    Ok(Proof {
        claim,
        interaction_pow,
        interaction_claim,
        stark_proof: proof,
    })
}
