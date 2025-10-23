use rookie::prover::prove_dummy;
use stwo::core::vcs::blake2_merkle::Blake2sMerkleChannel;

#[test_log::test]
fn test_prove_dummy() {
    let result = prove_dummy::<Blake2sMerkleChannel, 10>(4);
    assert!(result.is_ok());
}
