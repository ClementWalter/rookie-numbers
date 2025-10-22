mod big_sigma_0;
mod big_sigma_1;
mod ch_left;
mod ch_right;
mod maj;
mod range_check_add;
mod sigma_0;
mod sigma_1;

use crate::components;

components!(
    sigma_0::i0_i1,
    sigma_1::i0_i1,
    sigma_0::o2,
    sigma_1::o2,
    big_sigma_0::i0_i1,
    big_sigma_0::o2,
    big_sigma_1::i0,
    big_sigma_1::i1,
    big_sigma_1::o2,
    ch_left::i0,
    ch_left::i1,
    ch_right::i0,
    ch_right::i1,
    maj::i0h0_i1l0,
    maj::i0h1_i1l1,
    maj::i0l_i1h,
    range_check_add::range_check_add,
);
