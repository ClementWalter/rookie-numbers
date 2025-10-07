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

use crate::preprocessed::{big_sigma_0, big_sigma_1, maj, sigma_0, sigma_1};
use crate::CHUNK_SIZE;
use crate::H_SIZE;

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

// U32 range checks
relation!(RangeCheckAdd, 10);
pub struct RangeCheckAddData {
    pub add4: [BaseColumn; 6], // [*_, result, carry]
    pub add6: [BaseColumn; 8], // [*_, result, carry]
    pub add7: [BaseColumn; 9], // [*_, result, carry]
}

#[derive(Clone)]
pub struct Relations {
    pub sigma0: sigma_0::Relation,
    pub sigma1: sigma_1::Relation,
    pub big_sigma_0: big_sigma_0::Relation,
    pub big_sigma_1: big_sigma_1::Relation,
    pub ch_left: ChLeft,
    pub ch_right: ChRight,
    pub maj: maj::Relation,
    pub range_check_add: RangeCheckAdd,
}

impl Relations {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            sigma0: sigma_0::Relation::draw(channel),
            sigma1: sigma_1::Relation::draw(channel),
            big_sigma_0: big_sigma_0::Relation::draw(channel),
            big_sigma_1: big_sigma_1::Relation::draw(channel),
            ch_left: ChLeft::draw(channel),
            ch_right: ChRight::draw(channel),
            maj: maj::Relation::draw(channel),
            range_check_add: RangeCheckAdd::draw(channel),
        }
    }

    pub fn dummy() -> Self {
        Self {
            sigma0: sigma_0::Relation::dummy(),
            sigma1: sigma_1::Relation::dummy(),
            big_sigma_0: big_sigma_0::Relation::dummy(),
            big_sigma_1: big_sigma_1::Relation::dummy(),
            ch_left: ChLeft::dummy(),
            ch_right: ChRight::dummy(),
            maj: maj::Relation::dummy(),
            range_check_add: RangeCheckAdd::dummy(),
        }
    }
}

pub struct LookupData {
    pub initial_state: [BaseColumn; CHUNK_SIZE],
    pub final_state: [BaseColumn; H_SIZE],
    pub sigma_0: sigma_0::LookupData,
    pub sigma_1: sigma_1::LookupData,
    pub big_sigma_0: big_sigma_0::LookupData,
    pub big_sigma_1: big_sigma_1::LookupData,
    pub ch_left: ChLeftLookupData,
    pub ch_right: ChRightLookupData,
    pub maj: maj::LookupData,
    pub range_check_add: RangeCheckAddData,
}

impl LookupData {
    pub fn new(log_size: u32) -> Self {
        Self {
            initial_state: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            final_state: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            sigma_0: sigma_0::LookupData {
                i0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            sigma_1: sigma_1::LookupData {
                i0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            big_sigma_0: big_sigma_0::LookupData {
                i0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                o2: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
            },
            big_sigma_1: big_sigma_1::LookupData {
                i0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
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
            maj: maj::LookupData {
                i0_l: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h_0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i0_h_1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l_0: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
                i1_l_1: std::array::from_fn(|_| BaseColumn::zeros(1 << log_size)),
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
