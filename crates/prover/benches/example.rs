use prover::prover::prove_rookie;
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

fn main() {
    divan::main();
}

#[divan::bench(args = [15, 16, 17, 18, 19, 20], sample_count = 20)]
fn single_component(n: u32) {
    let result = prove_rookie::<Blake2sMerkleChannel>(n);
    assert!(result.is_ok());
}
