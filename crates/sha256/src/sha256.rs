//! SHA-256 functions
pub const fn small_sigma0(x: u32) -> u32 {
    x.rotate_right(7) ^ x.rotate_right(18) ^ (x >> 3)
}

pub const fn small_sigma1(x: u32) -> u32 {
    x.rotate_right(17) ^ x.rotate_right(19) ^ (x >> 10)
}

pub const fn big_sigma0(x: u32) -> u32 {
    x.rotate_right(2) ^ x.rotate_right(13) ^ x.rotate_right(22)
}

pub const fn big_sigma1(x: u32) -> u32 {
    x.rotate_right(6) ^ x.rotate_right(11) ^ x.rotate_right(25)
}

pub const fn ch_left(e: u32, f: u32) -> u32 {
    e & f
}

pub const fn ch_right(e: u32, g: u32) -> u32 {
    (0xFFFFFFFF - e) & g
}

pub const fn maj(a: u32, b: u32, c: u32) -> u32 {
    (a & b) ^ (a & c) ^ (b & c)
}

pub fn process_chunk(chunk: [u32; 16], mut hash: [u32; 8]) -> [u32; 8] {
    let mut w: [u32; 64] = [0; 64];
    w[..16].copy_from_slice(&chunk);

    // Schedule
    for t in 16..64 {
        w[t] = w[t - 16] + small_sigma0(w[t - 15]) + w[t - 7] + small_sigma1(w[t - 2])
    }

    // Compression
    let mut a = hash[0];
    let mut b = hash[1];
    let mut c = hash[2];
    let mut d = hash[3];
    let mut e = hash[4];
    let mut f = hash[5];
    let mut g = hash[6];
    let mut h = hash[7];
    for &wt in w.iter() {
        let temp1 = h + big_sigma1(e) + ch_left(e, f) + ch_right(e, g) + wt;
        let temp2 = big_sigma0(a) + maj(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }
    hash[0] += a;
    hash[1] += b;
    hash[2] += c;
    hash[3] += d;
    hash[4] += e;
    hash[5] += f;
    hash[6] += g;
    hash[7] += h;
    hash
}
