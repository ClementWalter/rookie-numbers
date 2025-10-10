#!/usr/bin/env python3.12

"""
Notebook to explain and actually compute parameters for the SHA-256 AIR.

Summary:
- s0, s1 and S1 cost 10T + 3L
- S0 costs 12T + 3L
- ch costs 10T + 4L
- maj costs 14T + 6L
- addition costs 2T + 2L
- schedule: 48 * (s0 + s1 + add) = 48 * (10T + 3L + 10T + 3L + 2T + 2L) = 48 * (22T + 8L) = 1056T + 384L
- compression: 64 * (S0 + S1 + ch + maj + 2*add) = 64 * (12T + 3L + 10T + 3L + 10T + 4L + 14T + 6L + 2T + 2L + 2T + 2L) = 64 * (50T + 20L) = 3200T + 1280L

Total: 4256T + 1664L

Using L = 2T, we get 7584 cells
"""

# %% Imports
import random

# %% SHA256 functions
SIGMA_PARAMS = {
    "small_sigma0": (7, 18, 3),
    "small_sigma1": (17, 19, 10),
    "big_sigma0": (2, 13, 22),
    "big_sigma1": (6, 11, 25),
}


def rotr(x, n):
    """
    Rotate right by n bits
    """
    return (x >> n) | (x << (32 - n)) & 0xFFFFFFFF


def small_sigma(x, *params):
    """
    SHA-256 σ₀ function used in message scheduling.
    """
    return rotr(x, params[0]) ^ rotr(x, params[1]) ^ (x >> params[2])


def big_sigma(x, *params):
    """
    SHA-256 Σ₁ function used in the compression function.
    """
    return rotr(x, params[0]) ^ rotr(x, params[1]) ^ rotr(x, params[2])


def sigma(x, *params):
    """
    SHA-256 σ function used in message scheduling.
    """
    if set(params) == {7, 18, 3}:
        return small_sigma(x, 7, 18, 3)
    elif set(params) == {17, 19, 10}:
        return small_sigma(x, 17, 19, 10)
    elif set(params) == {6, 11, 25}:
        return big_sigma(x, 6, 11, 25)
    elif set(params) == {2, 13, 22}:
        return big_sigma(x, 2, 13, 22)
    else:
        raise ValueError(f"Invalid parameters: {params}")


def ch(x, y, z):
    """
    SHA-256 Ch (choice) function: if x then y else z.
    Usual implementation is (x & y) ^ ((0xFFFFFFFF - x) & z), which is actually
    the same (x & y) + ((0xFFFFFFFF - x) & z).
    For the AIR, we split the function into two parts as it reduces the size of the
    lookup tables.
    """
    return ch_left(x, y) + ch_right(x, z)


def ch_left(x, y):
    """
    The left part of the Ch function.
    """
    return x & y


def ch_right(x, z):
    """
    The right part of the Ch function.
    """
    return (0xFFFFFFFF - x) & z


def maj(x, y, z):
    """
    SHA-256 Maj (majority) function: majority vote of x, y, z.
    """
    return (x & y) ^ (x & z) ^ (y & z)


def process_chunk(chunk, H):
    """
    Process a single 512-bit (16-word) chunk of the message.
    """

    assert len(chunk) == 16
    assert len(H) == 8

    W = chunk + [0] * 48
    for t in range(16, 64):
        W[t] = (
            W[t - 16]
            + sigma(W[t - 15], *SIGMA_PARAMS["small_sigma0"])
            + W[t - 7]
            + sigma(W[t - 2], *SIGMA_PARAMS["small_sigma1"])
        )

    a, b, c, d, e, f, g, h = H

    for t in range(64):
        temp1 = (
            h
            + sigma(e, *SIGMA_PARAMS["big_sigma1"])
            + ch_left(e, f)
            + ch_right(e, g)
            + W[t]
        )
        temp2 = sigma(a, *SIGMA_PARAMS["big_sigma0"]) + maj(a, b, c)
        h = g
        g = f
        f = e
        e = d + temp1
        d = c
        c = b
        b = a
        a = temp1 + temp2
    H[0] = H[0] + a
    H[1] = H[1] + b
    H[2] = H[2] + c
    H[3] = H[3] + d
    H[4] = H[4] + e
    H[5] = H[5] + f
    H[6] = H[6] + g
    H[7] = H[7] + h
    return H


# %% Compute input space splits
#'
#' The main goal with splitting the 32-bits input space is to be able to compute bit operations
#' with lookups no bigger than 2**21. Actually, there are four expressions that need to be emulated
#' using lookups:
#' - small_sigma0(w)
#' - small_sigma1(w)
#' - big_sigma0(a) + maj(a, b, c)
#' - big_sigma1(e) + ch(e, f, g)
#'
#' Sigma computations are defined as XOR over the bit shifted or bit rotated input. Given the periodicity
#' of the problem (bit i+n0 is XORed with bit i+n1 and i+n2),
#' we can derive 2 separated partitions I0 and I1 following the 3D lattice:
def lattice(*params):
    """
    Given the 3 parameters, return the best partition of the 32-bits
    input space according to the 3D periodicity constraint, in order to
    compute small and big sigma in chunks.
    """
    n1, n2, n3 = params
    selected = sum(
        2**i
        for i in {
            ((n2 - n1) * a + (n3 - n1) * b) % 32 for a in range(4) for b in range(4)
        }
    )
    remaining = (2**32 - 1) ^ selected

    return selected, remaining


def get_subset(mask):
    """
    Generator of all the numbers in the mask
    """
    current = mask
    while True:
        yield current
        if current == 0:
            break
        current = (current - 1) & mask


def get_output_bits(partition, *params):
    """
    Iterate over the whole partition to return the set of bit of the output space
    actually modified by the partition
    """
    res = 0
    for x in get_subset(partition):
        res |= sigma(x, *params)
    return res


def valid_split(i0, i1, o2):
    """
    Split needs to satisfy:
    - at least 2 of the 4 subsets with less than 8 bits
    - the total number of bits in the subsets is less than 21
    - the number of bits in the output space is less than 10
    """
    return (
        (
            ((i0 % 2**16).bit_count() < 8)
            + ((i0 >> 16).bit_count() < 8)
            + ((i1 % 2**16).bit_count() < 8)
            + ((i1 >> 16).bit_count() < 8)
            >= 2
        )
        and max(i0.bit_count(), i1.bit_count()) <= 21
        and o2.bit_count() <= 10
    )


def get_split(*params):
    for _ in range(3):
        I0, I1 = lattice(*params)
        O0 = get_output_bits(I0, *params)
        O1 = get_output_bits(I1, *params)
        O2 = O0 & O1
        O0 = O0 & ~O2
        O1 = O1 & ~O2
        if valid_split(I0, I1, O2):
            break
        params = (*params[1:], params[0])

    assert valid_split(I0, I1, O2)
    return I0, I1, O0, O1, O2


def pext(x, mask):
    out = 0
    bb = 1
    while mask != 0:
        ls = mask & (2**32 - mask)
        if x & ls != 0:
            out |= bb
        mask ^= ls
        bb <<= 1
    return out


# %% Scheduling

#' Then, for each partition I0 and I1, we evaluate how it propagates through the sigma functions in terms of output bits.
#' In other words, we find which output bits are affected by each partition I0 and I1.
#'
#' The 32-bits output space is consequently split into 3 parts:
#' - O0: bits only affected by partition I0
#' - O1: bits only affected by partition I1
#' - O2: bits affected by both partitions I0 and I1
params = SIGMA_PARAMS["small_sigma0"]
I0, I1, O0, O1, O2 = get_split(*params)

#' All these partitions also need to be split with low and high bits in order to be able
#' to go back to the usual u32 = (u16, u16) representation.
#' I0 = I0_L | I0_H
#' I1 = I1_L | I1_H
#' O0 = O0_L | O0_H
#' O1 = O1_L | O1_H
#' O2 = O2_L | O2_H

I0_L = I0 & 0xFFFF
I0_H = I0 >> 16
I1_L = I1 & 0xFFFF
I1_H = I1 >> 16

O0_L = O0 & 0xFFFF
O0_H = O0 >> 16
O1_L = O1 & 0xFFFF
O1_H = O1 >> 16
O2_L = O2 & 0xFFFF
O2_H = O2 >> 16

print(
    f"Params: {params}\n\n"
    f"I0: {I0:032b}: {I0.bit_count()} bits\n"
    f"I1: {I1:032b}: {I1.bit_count()} bits\n\n"
    f"I0_L: {I0_L:016b}: {I0_L.bit_count()} bits\n"
    f"I1_L: {I1_L:016b}: {I1_L.bit_count()} bits\n"
    f"I0_H: {I0_H:016b}: {I0_H.bit_count()} bits\n"
    f"I1_H: {I1_H:016b}: {I1_H.bit_count()} bits\n"
    f"O0: {O0:032b}: {O0.bit_count()} bits\n"
    f"O1: {O1:032b}: {O1.bit_count()} bits\n"
    f"O2: {O2:032b}: {O2.bit_count()} bits\n\n"
    f"O0_L: {O0_L:016b}: {O0_L.bit_count()} bits\n"
    f"O1_L: {O1_L:016b}: {O1_L.bit_count()} bits\n"
    f"O2_L: {O2_L:016b}: {O2_L.bit_count()} bits\n"
    f"O0_H: {O0_H:016b}: {O0_H.bit_count()} bits\n"
    f"O1_H: {O1_H:016b}: {O1_H.bit_count()} bits\n"
    f"O2_H: {O2_H:016b}: {O2_H.bit_count()} bits"
)

#' At this point small sigmas can be computed as follows:
#' - split w in 4 chunks: I0_L | I1_L | I0_H | I1_H
#' - O0_L, O0_H, O20 = LookupSmallSigma0(I0_L, I0_H)
#' - O1_L, O1_H, O21 = LookupSmallSigma1(I1_L, I1_H)
#' - O2_L, O2_H = LookupXOR(O20, O21)
lookup_sigma_0 = {
    (low, high): (
        sigma(low + 2**16 * high, *params) & O0_L,
        (sigma(low + 2**16 * high, *params) >> 16) & O0_H,
        sigma(low + 2**16 * high, *params) & O2,
    )
    for low in get_subset(I0_L)
    for high in get_subset(I0_H)
}
lookup_sigma_1 = {
    (low, high): (
        sigma(low + 2**16 * high, *params) & O1_L,
        (sigma(low + 2**16 * high, *params) >> 16) & O1_H,
        sigma(low + 2**16 * high, *params) & O2,
    )
    for low in get_subset(I1_L)
    for high in get_subset(I1_H)
}
lookup_xor = {
    (o20, o21): ((o20 ^ o21) % 2**16, (o20 ^ o21) >> 16)
    for o20 in get_subset(O2)
    for o21 in get_subset(O2)
}

#' Wrap up and test the scheduling
x = random.randint(0, 2**32 - 1)
#' Values are all stores as (u16, u16)
x_low, x_high = (x % 2**16, x >> 16)

#' Guess split in (input_0, input_1) for each part: 2T
x_0_low, x_0_high = (x_low & I0_L, x_high & I0_H)
x_1_low = x_low - x_0_low
x_1_high = x_high - x_0_high

#' Lookup intermediate values (this also range check the 4 limbs): 8T + 3L
o0_l, o0_h, o20 = lookup_sigma_0[(x_0_low, x_0_high)]
o1_l, o1_h, o21 = lookup_sigma_1[(x_1_low, x_1_high)]
o2_l, o2_h = lookup_xor[(o20, o21)]

#' Debug test
assert sigma(x, *params) == o0_l + o1_l + o2_l + 2**16 * (o0_h + o1_h + o2_h)

#' These are range-checked by design
s_l = o0_l + o1_l + o2_l
s_h = o0_h + o1_h + o2_h

#' Fake W values for linter
w_l = w_h = 0

#' Guess carry_low and carry_high for the final sum: 2T
carry_low = (s_l + w_l + s_l + w_l) >> 16
carry_high = (s_h + w_h + s_h + w_h) >> 16

range_check_16_carry = {(v, c) for v in range(2**16) for c in range(4)}

#' Range-check final parts: 2L
assert (s_l + w_l + s_l + w_l - 2**16 * carry_low, carry_low) in range_check_16_carry
assert (s_h + w_h + s_h + w_h - 2**16 * carry_high, carry_high) in range_check_16_carry

# %% Compression: 30T + 13L
#' Given the fact that maj and ch take 3 operands, we need for these terms to split
#' the input space into chunks smaller than 7 bits (3*7 = 21 bits).
# params = (11, 25, 6)
params = SIGMA_PARAMS["big_sigma0"]
I0, I1, O0, O1, O2 = get_split(*params)

I0_L = I0 & 0xFFFF
I0_L0 = I0_L & 0xFF
I0_L1 = (I0_L >> 8) & 0xFF

I0_H = I0 >> 16
I0_H0 = I0_H & 0xFF
I0_H1 = (I0_H >> 8) & 0xFF

I1_L = I1 & 0xFFFF
I1_L0 = I1_L & 0xFF
I1_L1 = (I1_L >> 8) & 0xFF

I1_H = I1 >> 16
I1_H0 = I1_H & 0xFF
I1_H1 = (I1_H >> 8) & 0xFF

O0_L = O0 & 0xFFFF
O0_H = O0 >> 16
O1_L = O1 & 0xFFFF
O1_H = O1 >> 16
O2_L = O2 & 0xFFFF
O2_H = O2 >> 16

print(
    f"Params: {params}\n\n"
    f"I0: {I0:032b}: {I0.bit_count()} bits\n"
    f"I1: {I1:032b}: {I1.bit_count()} bits\n\n"
    f"I0_L: {I0_L:016b}: {I0_L.bit_count()} bits\n"
    f"I0_L0: {I0_L0:08b}: {I0_L0.bit_count()} bits\n"
    f"I0_L1: {I0_L1:08b}: {I0_L1.bit_count()} bits\n"
    f"I0_H: {I0_H:016b}: {I0_H.bit_count()} bits\n"
    f"I0_H0: {I0_H0:08b}: {I0_H0.bit_count()} bits\n"
    f"I0_H1: {I0_H1:08b}: {I0_H1.bit_count()} bits\n"
    f"I1_L: {I1_L:016b}: {I1_L.bit_count()} bits\n"
    f"I1_L0: {I1_L0:08b}: {I1_L0.bit_count()} bits\n"
    f"I1_L1: {I1_L1:08b}: {I1_L1.bit_count()} bits\n"
    f"I1_H: {I1_H:016b}: {I1_H.bit_count()} bits\n"
    f"I1_H0: {I1_H0:08b}: {I1_H0.bit_count()} bits\n"
    f"I1_H1: {I1_H1:08b}: {I1_H1.bit_count()} bits\n"
    f"O0: {O0:032b}: {O0.bit_count()} bits\n"
    f"O1: {O1:032b}: {O1.bit_count()} bits\n"
    f"O2: {O2:032b}: {O2.bit_count()} bits\n\n"
    f"O0_L: {O0_L:016b}: {O0_L.bit_count()} bits\n"
    f"O1_L: {O1_L:016b}: {O1_L.bit_count()} bits\n"
    f"O2_L: {O2_L:016b}: {O2_L.bit_count()} bits\n"
    f"O0_H: {O0_H:016b}: {O0_H.bit_count()} bits\n"
    f"O1_H: {O1_H:016b}: {O1_H.bit_count()} bits\n"
    f"O2_H: {O2_H:016b}: {O2_H.bit_count()} bits"
)

#' At this point big sigmas can be computed as follows:
#' - split e in 6 chunks: I0_L | I1_L0 | I1_L1 | I0_H0 | I0_H1 | I1_H
#' - O0_L, O0_H, O20 = LookupBigSigma0(I0_L, I0_H0, I0_H1)
#' - O1_L, O1_H, O21 = LookupBigSigma1(I1_L0, I1_L1, I1_H)
#' - O2_L, O2_H = LookupXOR(O20, O21)
lookup_big_sigma0 = {
    (low, high_0, high_1): (
        sigma(low + 2**16 * (high_0 + 2**8 * high_1), *params) & O0_L,
        (sigma(low + 2**16 * (high_0 + 2**8 * high_1), *params) >> 16) & O0_H,
        sigma(low + 2**16 * (high_0 + 2**8 * high_1), *params) & O2,
    )
    for low in get_subset(I0_L)
    for high_0 in get_subset(I0_H0)
    for high_1 in get_subset(I0_H1)
}
lookup_big_sigma1 = {
    (low_0, low_1, high): (
        sigma(low_0 + 2**8 * low_1 + 2**16 * high, *params) & O1_L,
        (sigma(low_0 + 2**8 * low_1 + 2**16 * high, *params) >> 16) & O1_H,
        sigma(low_0 + 2**8 * low_1 + 2**16 * high, *params) & O2,
    )
    for low_0 in get_subset(I1_L0)
    for low_1 in get_subset(I1_L1)
    for high in get_subset(I1_H)
}
lookup_xor = {
    (o20, o21): ((o20 ^ o21) % 2**16, (o20 ^ o21) >> 16)
    for o20 in get_subset(O2)
    for o21 in get_subset(O2)
}

#' ch and maj are just parallelized over the split input space
lookup_ch = {
    partition: {
        (x, y, z): ch(x, y, z)
        for x in get_subset(partition)
        for y in get_subset(partition)
        for z in get_subset(partition)
    }
    for partition in [I0_L, I0_H0, I0_H1, I1_L0, I1_L1, I1_H]
}
lookup_maj = {
    partition: {
        (x, y, z): maj(x, y, z)
        for x in get_subset(partition)
        for y in get_subset(partition)
        for z in get_subset(partition)
    }
    for partition in [I0_L, I0_H0, I0_H1, I1_L0, I1_L1, I1_H]
}

#' Wrap up and test the scheduling
e = random.randint(0, 2**32 - 1)
#' Values are all stores as (u16, u16)
e_low, e_high = (e % 2**16, e >> 16)

#' Guess split in (input_0, input_1) for each part: 4T
e_0_high_0, e_0_high_1 = e_high & I0_H0, (e_high >> 8) & I0_H1
e_1_low_0, e_1_low_1 = e_low & I1_L0, (e_low >> 8) & I1_L1
e_1_high = e_high - e_0_high_0 - 2**8 * e_0_high_1
e_0_low = e_low - e_1_low_0 - 2**8 * e_1_low_1

#' Lookup intermediate values (this also range check the 6 limbs): 8T + 3L
o0_l, o0_h, o20 = lookup_big_sigma0[(e_0_low, e_0_high_0, e_0_high_1)]
o1_l, o1_h, o21 = lookup_big_sigma1[(e_1_low_0, e_1_low_1, e_1_high)]
o2_l, o2_h = lookup_xor[(o20, o21)]

#' Debug test
assert sigma(e, *params) == o0_l + o1_l + o2_l + 2**16 * (o0_h + o1_h + o2_h)

#' These are range-checked by design
s_l = o0_l + o1_l + o2_l
s_h = o0_h + o1_h + o2_h

#' Get f and g values
f = random.randint(0, 2**32 - 1)
g = random.randint(0, 2**32 - 1)
#' Values are all stores as (u16, u16)
f_low, f_high = (f % 2**16, f >> 16)
g_low, g_high = (g % 2**16, g >> 16)

#' Guess split in (input_0, input_1) for each part: 2*4T = 8T
f_0_high_0, f_0_high_1 = f_high & I0_H0, (f_high >> 8) & I0_H1
f_1_low_0, f_1_low_1 = f_low & I1_L0, (f_low >> 8) & I1_L1
f_1_high = f_high - f_0_high_0 - 2**8 * f_0_high_1
f_0_low = f_low - f_1_low_0 - 2**8 * f_1_low_1
g_0_high_0, g_0_high_1 = g_high & I0_H0, (g_high >> 8) & I0_H1
g_1_low_0, g_1_low_1 = g_low & I1_L0, (g_low >> 8) & I1_L1
g_1_high = g_high - g_0_high_0 - 2**8 * g_0_high_1
g_0_low = g_low - g_1_low_0 - 2**8 * g_1_low_1

#' Lookup ch values: 6T + 6L
ch_0_low = lookup_ch[I0_L][(e_0_low, f_0_low, g_0_low)]
ch_0_high_0 = lookup_ch[I0_H0][(e_0_high_0, f_0_high_0, g_0_high_0)]
ch_0_high_1 = lookup_ch[I0_H1][(e_0_high_1, f_0_high_1, g_0_high_1)]
ch_1_low_0 = lookup_ch[I1_L0][(e_1_low_0, f_1_low_0, g_1_low_0)]
ch_1_low_1 = lookup_ch[I1_L1][(e_1_low_1, f_1_low_1, g_1_low_1)]
ch_1_high = lookup_ch[I1_H][(e_1_high, f_1_high, g_1_high)]

#' These are range-checked by design
ch_l = ch_0_low + ch_1_low_0 + 2**8 * ch_1_low_1
ch_h = ch_0_high_0 + ch_0_high_1 + 2**8 * ch_1_high

#' Fake d, h, K and W values for linter
d_l = d_h = 0
h_l = h_h = 0
k_l = k_h = 0
w_l = w_h = 0

#' Guess the carries for the two final sums: 2*2T = 4T
#' Don't range-check temp1 and temp2, just the final sum
carry_low = 0
carry_high = 0

range_check_16_carry = {(v, c) for v in range(2**16) for c in range(5)}

#' Range-check final low and high: 2*2L = 4L
assert (
    d_l + h_l + s_l + ch_l + k_l + w_l - 2**16 * carry_low,
    carry_low,
) in range_check_16_carry
assert (
    d_h + h_h + s_h + ch_h + k_h + w_h - 2**16 * carry_high,
    carry_high,
) in range_check_16_carry
