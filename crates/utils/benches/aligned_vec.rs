use utils::aligned_vec;

const LEN: usize = 1_000_000;
const VALUE: u32 = 42;

#[divan::bench(name = "vec_repeat")]
fn vec_repeat() {
    let v = vec![VALUE; LEN];
    divan::black_box(v);
}

#[divan::bench(name = "aligned_vec_repeat")]
fn aligned_vec_repeat() {
    let v = aligned_vec![VALUE; LEN];
    divan::black_box(v);
}

fn main() {
    divan::main();
}
