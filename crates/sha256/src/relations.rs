//! Relations for the SHA-256 AIR.
//! To alleviate the amount of relation, we use one single relation for each function and keep
//! one slot on the relation for the function part identifier, e.g.
//!
//! Sigma0(O1, ...values)
//! Sigma0(O2, ...values)
//! Sigma0(O3, ...values)

//! instead of managing 3 independent relations
//! Sigma0_O1(...values)
//! Sigma0_O2(...values)
//! Sigma0_O3(...values)

use stwo_prover::core::backend::simd::column::BaseColumn;
use stwo_prover::core::backend::Column;
use stwo_prover::core::channel::Channel;
use stwo_prover::relation;

use crate::CHUNK_SIZE;
use crate::H_SIZE;

/// This will be used to identify the function part in the relation.
#[repr(u32)]
enum Partition {
    O1,
    O2,
    O3,
    I0_L,
    I0_H,
    I1_L,
    I1_H,
    I0_H0,
    I0_H1,
    I1_L0,
    I1_L1,
    ADD4,
    ADD6,
    ADD7,
}

// Scheduling lookups
relation!(Sigma0, 6);
pub struct Sigma0LookupData {
    pub o0: [BaseColumn; 5], // [i0_l, i0_h, o0_l, o0_h, o20]
    pub o1: [BaseColumn; 5], // [i1_l, i1_h, o1_l, o1_h, o21]
    pub o2: [BaseColumn; 4], // [o20, o21, o2_l, o2_h]
}

relation!(Sigma1, 6);
pub struct Sigma1LookupData {
    pub o0: [BaseColumn; 5], // [i0_l, i0_h, o0_l, o0_h, o20]
    pub o1: [BaseColumn; 5], // [i1_l, i1_h, o1_l, o1_h, o21]
    pub o2: [BaseColumn; 4], // [o20, o21, o2_l, o2_h]
}

// Compression lookups
relation!(BigSigma0, 7);
pub struct BigSigma0LookupData {
    pub o0: [BaseColumn; 6], // [i0_l, i0_h_0, i0_h_1, o0_l, o0_h, o20]
    pub o1: [BaseColumn; 6], // [i1_l_0, i1_l_1, i1_h, o1_l, o1_h, o21]
    pub o2: [BaseColumn; 4], // [o20, o21, o2_l, o2_h]
}

relation!(BigSigma1, 6);
pub struct BigSigma1LookupData {
    pub o0: [BaseColumn; 5], // [i0_l, i0_h, o0_l, o0_h, o20]
    pub o1: [BaseColumn; 5], // [i1_l, i1_h, o1_l, o1_h, o21]
    pub o2: [BaseColumn; 4], // [o20, o21, o2_l, o2_h]
}

relation!(ChLeft, 4);
pub struct ChLeftLookupData {
    pub i0_l: [BaseColumn; 3], // [e, f, result]
    pub i0_h: [BaseColumn; 3], // [e, f, result]
    pub i1_l: [BaseColumn; 3], // [e, f, result]
    pub i1_h: [BaseColumn; 3], // [e, f, result]
}

relation!(ChRight, 4);
pub struct ChRightLookupData {
    pub i0_l: [BaseColumn; 3], // [e, g, result]
    pub i0_h: [BaseColumn; 3], // [e, g, result]
    pub i1_l: [BaseColumn; 3], // [e, g, result]
    pub i1_h: [BaseColumn; 3], // [e, g, result]
}

relation!(Maj, 5);
pub struct MajLookupData {
    pub i0_l: [BaseColumn; 4],  // [a, b, c, result]
    pub i0_h0: [BaseColumn; 4], // [a, b, c, result]
    pub i0_h1: [BaseColumn; 4], // [a, b, c, result]
    pub i1_l0: [BaseColumn; 4], // [a, b, c, result]
    pub i1_l1: [BaseColumn; 4], // [a, b, c, result]
    pub i1_h: [BaseColumn; 4],  // [a, b, c, result]
}

// U32 range checks
relation!(RangeCheckAdd, 10);
pub struct RangeCheckAddData {
    pub add4: [BaseColumn; 6], // [*_, result, carry]
    pub add6: [BaseColumn; 8], // [*_, result, carry]
    pub add7: [BaseColumn; 9], // [*_, result, carry]
}

#[derive(Clone)]
pub struct Relations {
    pub sigma0: Sigma0,
    pub sigma1: Sigma1,
    pub big_sigma0: BigSigma0,
    pub big_sigma1: BigSigma1,
    pub ch_left: ChLeft,
    pub ch_right: ChRight,
    pub maj: Maj,
    pub range_check_add: RangeCheckAdd,
}

impl Relations {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            sigma0: Sigma0::draw(channel),
            sigma1: Sigma1::draw(channel),
            big_sigma0: BigSigma0::draw(channel),
            big_sigma1: BigSigma1::draw(channel),
            ch_left: ChLeft::draw(channel),
            ch_right: ChRight::draw(channel),
            maj: Maj::draw(channel),
            range_check_add: RangeCheckAdd::draw(channel),
        }
    }
}

pub struct LookupData {
    pub initial_state: [BaseColumn; CHUNK_SIZE],
    pub final_state: [BaseColumn; H_SIZE],
    pub sigma_0: Sigma0LookupData,
    pub sigma_1: Sigma1LookupData,
    pub big_sigma_0: BigSigma0LookupData,
    pub big_sigma_1: BigSigma1LookupData,
    pub ch_left: ChLeftLookupData,
    pub ch_right: ChRightLookupData,
    pub maj: MajLookupData,
    pub range_check_add: RangeCheckAddData,
}

impl LookupData {
    pub fn new(log_size: u32) -> Self {
        Self {
            initial_state: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            final_state: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            sigma_0: Sigma0LookupData {
                o0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            sigma_1: Sigma1LookupData {
                o0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            big_sigma_0: BigSigma0LookupData {
                o0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            big_sigma_1: BigSigma1LookupData {
                o0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            ch_left: ChLeftLookupData {
                i0_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_h: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            ch_right: ChRightLookupData {
                i0_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_h: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            maj: MajLookupData {
                i0_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_h: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            range_check_add: RangeCheckAddData {
                add4: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                add6: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                add7: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
        }
    }
}
