use prover::prover::prove_rookie;
use stwo_prover::core::vcs::blake2_merkle::Blake2sMerkleChannel;

fn main() {
    divan::main();
}

const N: &[usize] = &[3, 3 * 10, 3 * 50, 3 * 150, 3 * 500, 3 * 1000, 3 * 2000];

#[divan::bench(consts = N, args = [15, 16, 17, 18, 19, 20], sample_count = 20)]
fn single_component<const N: usize>(log_size: u32) {
    let result = prove_rookie::<Blake2sMerkleChannel, N>(log_size);
    assert!(result.is_ok());
}
