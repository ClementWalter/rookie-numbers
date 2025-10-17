use itertools::Itertools;
use stwo::{
    core::{channel::Channel, fields::m31::BaseField, poly::circle::CanonicCoset},
    prover::{
        backend::simd::{column::BaseColumn, SimdBackend},
        poly::{circle::CircleEvaluation, BitReversedOrder},
    },
};
use stwo_constraint_framework::relation;

use crate::trace_columns;

// [value, carry]
const N_COLUMNS: usize = 2;

relation!(Add4, N_COLUMNS);
relation!(Add7, N_COLUMNS);
relation!(Add8, N_COLUMNS);

trace_columns!(
    Columns,
    range_check_add_value,
    range_check_add_carry_4,
    range_check_add_carry_7,
    range_check_add_carry_8,
);

#[derive(Debug, Clone)]
pub struct Relation {
    pub add_4: Add4,
    pub add_7: Add7,
    pub add_8: Add8,
}

impl Relation {
    pub fn dummy() -> Self {
        Self {
            add_4: Add4::dummy(),
            add_7: Add7::dummy(),
            add_8: Add8::dummy(),
        }
    }

    pub fn draw(channel: &mut impl Channel) -> Self {
        Self {
            add_4: Add4::draw(channel),
            add_7: Add7::draw(channel),
            add_8: Add8::draw(channel),
        }
    }
}

pub fn gen_column_simd() -> Vec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>> {
    let mut all_columns = Vec::with_capacity(N_COLUMNS);

    let domain = CanonicCoset::new(19).circle_domain();

    let columns = (0..1 << 16).flat_map(move |value| {
        let carry_4 = [0, 1, 2, 3, 0, 0, 0, 0].into_iter();
        let carry_7 = [0, 1, 2, 3, 4, 5, 6, 0].into_iter();
        let carry_8 = [0, 1, 2, 3, 4, 5, 6, 7].into_iter();
        carry_4
            .zip(carry_7)
            .zip(carry_8)
            .map(move |((c_4, c_7), c_8)| {
                (
                    BaseField::from_u32_unchecked(value),
                    BaseField::from_u32_unchecked(c_4),
                    BaseField::from_u32_unchecked(c_7),
                    BaseField::from_u32_unchecked(c_8),
                )
            })
    });

    let (limb, rest) = columns.tee();
    let (carry_4, rest) = rest.tee();
    let (carry_7, carry_8) = rest.tee();

    all_columns.extend(vec![
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain,
            BaseColumn::from_iter(limb.map(|t| t.0)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain,
            BaseColumn::from_iter(carry_4.map(|t| t.1)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain,
            BaseColumn::from_iter(carry_7.map(|t| t.2)),
        ),
        CircleEvaluation::<SimdBackend, BaseField, BitReversedOrder>::new(
            domain,
            BaseColumn::from_iter(carry_8.map(|t| t.3)),
        ),
    ]);

    all_columns
}

#[cfg(test)]
mod tests {
    use stwo::prover::backend::Column;

    use super::*;

    #[test]
    fn test_ids() {
        assert_eq!(Columns::to_ids().len(), 4);
    }

    #[test]
    fn test_gen_column_simd() {
        let columns = gen_column_simd();
        assert_eq!(columns.len(), 4);
        assert_eq!(columns[0].values.len().ilog2(), 19);
        assert_eq!(columns[1].values.len().ilog2(), 19);
        assert_eq!(columns[2].values.len().ilog2(), 19);
        assert_eq!(columns[3].values.len().ilog2(), 19);
    }
}
