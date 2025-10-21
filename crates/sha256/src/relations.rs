use stwo::core::channel::Channel;

use crate::preprocessed::{
    big_sigma_0, big_sigma_1, ch_left, ch_right, maj, range_check_add, sigma_0, sigma_1,
};

mod w {
    use stwo_constraint_framework::relation;

    use crate::components::W_SIZE;
    relation!(Relation, W_SIZE);
}

#[derive(Clone)]
pub struct Relations {
    pub sigma_0: sigma_0::Relation,
    pub sigma_1: sigma_1::Relation,
    pub big_sigma_0: big_sigma_0::Relation,
    pub big_sigma_1: big_sigma_1::Relation,
    pub ch_left: ch_left::Relation,
    pub ch_right: ch_right::Relation,
    pub maj: maj::Relation,
    pub range_check_add: range_check_add::Relation,
    pub w: w::Relation,
}

impl Relations {
    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            sigma_0: sigma_0::Relation::draw(channel),
            sigma_1: sigma_1::Relation::draw(channel),
            big_sigma_0: big_sigma_0::Relation::draw(channel),
            big_sigma_1: big_sigma_1::Relation::draw(channel),
            ch_left: ch_left::Relation::draw(channel),
            ch_right: ch_right::Relation::draw(channel),
            maj: maj::Relation::draw(channel),
            range_check_add: range_check_add::Relation::draw(channel),
            w: w::Relation::draw(channel),
        }
    }

    pub fn dummy() -> Self {
        Self {
            sigma_0: sigma_0::Relation::dummy(),
            sigma_1: sigma_1::Relation::dummy(),
            big_sigma_0: big_sigma_0::Relation::dummy(),
            big_sigma_1: big_sigma_1::Relation::dummy(),
            ch_left: ch_left::Relation::dummy(),
            ch_right: ch_right::Relation::dummy(),
            maj: maj::Relation::dummy(),
            range_check_add: range_check_add::Relation::dummy(),
            w: w::Relation::dummy(),
        }
    }
}
