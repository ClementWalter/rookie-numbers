use stwo::core::vcs::blake2_merkle::Blake2sMerkleChannel;
use vm::prover::{prove_rookie, verify_rookie};

#[test_log::test]
fn test_prove_rookie() {
    let result = prove_rookie::<Blake2sMerkleChannel, 3>(13);
    assert!(result.is_ok());
}

#[test_log::test]
fn test_verify_rookie() {
    let proof = prove_rookie::<Blake2sMerkleChannel, 3>(13).unwrap();
    let result = verify_rookie::<Blake2sMerkleChannel, 3>(proof);
    assert!(result.is_ok());
}
