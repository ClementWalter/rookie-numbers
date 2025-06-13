use prover::prover::{prove_rookie, verify_rookie};
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[test]
fn test_prove_rookie() {
    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .init();

    let result = prove_rookie::<Blake2sMerkleChannel, 3>(13);
    assert!(result.is_ok());
}

#[test]
fn test_verify_rookie() {
    let proof = prove_rookie::<Blake2sMerkleChannel, 3>(13).unwrap();
    let result = verify_rookie::<Blake2sMerkleChannel, 3>(proof);
    assert!(result.is_ok());
}
