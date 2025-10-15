#![allow(non_camel_case_types)]
use stwo_constraint_framework::relation;

relation!(Memory, 2); // addr, value

/// Assuming a 100-bit security target, the witness may
/// contain up to 1 << (24 + INTERACTION_POW_BITS) relation terms.
pub const INTERACTION_POW_BITS: u32 = 2;
