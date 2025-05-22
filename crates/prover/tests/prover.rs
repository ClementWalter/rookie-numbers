use prover::prover::prove_rookie;
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[test]
fn test_prove_rookie() {
    tracing_subscriber
        ::fmt()
        .with_max_level(tracing::Level::INFO) // Set the log level
        .init();

    let result = prove_rookie::<Blake2sMerkleChannel>((100u32).ilog2());
    assert!(result.is_ok());
}
