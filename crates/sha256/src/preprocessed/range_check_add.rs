use std::simd::u32x16;

use stwo::core::channel::Channel;
use stwo_constraint_framework::relation;
use utils::trace_columns;

// [value, carry]
const N_COLUMNS: usize = 2;

relation!(Add4, N_COLUMNS);
relation!(Add7, N_COLUMNS);
relation!(Add8, N_COLUMNS);

trace_columns!(RangeCheckAddColumns, value, carry_4, carry_7, carry_8,);

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

pub fn gen_column_simd() -> Vec<Vec<u32x16>> {
    const N: usize = 1 << 15;
    let mut all_columns = vec![
        Vec::with_capacity(N),
        Vec::with_capacity(N),
        Vec::with_capacity(N),
        Vec::with_capacity(N),
    ];

    let carry_4 = u32x16::from_array([0, 1, 2, 3, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0]);
    let carry_7 = u32x16::from_array([0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 2, 3, 4, 5, 6, 0]);
    let carry_8 = u32x16::from_array([0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7]);

    for i in 0..1 << 15 {
        let value = u32x16::from_array(std::array::from_fn(|j| i * 2 + (j > 7) as u32));
        all_columns[0].push(value);
        all_columns[1].push(carry_4);
        all_columns[2].push(carry_7);
        all_columns[3].push(carry_8);
    }

    all_columns
}

#[cfg(test)]
mod tests {
    use stwo::prover::backend::simd::m31::LOG_N_LANES;

    use super::*;

    #[test]
    fn test_ids() {
        assert_eq!(RangeCheckAddColumns::SIZE, 4);
    }

    #[test]
    fn test_gen_column_simd() {
        let columns = gen_column_simd();
        assert_eq!(columns.len(), RangeCheckAddColumns::SIZE);
        assert_eq!(columns[0].len().ilog2(), 19 - LOG_N_LANES);
        assert_eq!(columns[1].len().ilog2(), 19 - LOG_N_LANES);
        assert_eq!(columns[2].len().ilog2(), 19 - LOG_N_LANES);
        assert_eq!(columns[3].len().ilog2(), 19 - LOG_N_LANES);
    }
}
