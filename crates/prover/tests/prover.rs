use prover::prover::prove_rookie;
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[test]
fn test_prove_rookie() {
    // Initialize tracing for this test.
    // Use try_init to avoid panic if already initialized.
    // The .ok() silences the unused result warning.
    let _ = tracing_subscriber::fmt().try_init();

    let result = prove_rookie::<Blake2sMerkleChannel>((1_000u32).ilog2());
    assert!(result.is_ok());
}
