use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::range_check_add::columns::ComponentColumnsOwned as ComponentColumns,
    preprocessed::range_check_add::RangeCheckAddColumnsOwned as RangeCheckAddColumns,
    relations::Relations,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumns {
        carry_4_mult,
        carry_7_mult,
        carry_8_mult,
    } = ComponentColumns::<<E as EvalAtRow>::F>::from_eval(eval);
    let RangeCheckAddColumns {
        value,
        carry_4,
        carry_7,
        carry_8,
    } = RangeCheckAddColumns::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.range_check_add.add_4,
        E::EF::from(carry_4_mult),
        value,
        carry_4,
    );
    add_to_relation!(
        eval,
        relations.range_check_add.add_7,
        E::EF::from(carry_7_mult),
        value,
        carry_7,
    );
    add_to_relation!(
        eval,
        relations.range_check_add.add_8,
        E::EF::from(carry_8_mult),
        value,
        carry_8,
    );
    eval.finalize_logup_in_pairs();
}

#[derive(Clone)]
pub struct Eval {
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        19
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        19 + 1
    }
    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        eval_constraints(&mut eval, &self.relations);
        eval
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use stwo::{
        core::{pcs::TreeVec, poly::circle::CanonicCoset},
        prover::backend::simd::m31::LOG_N_LANES,
    };
    use stwo_constraint_framework::assert_constraints_on_polys;
    use utils::circle_evaluation_u32x16;

    use super::*;
    use crate::{
        components::{
            compression::witness::gen_trace as gen_compression_trace,
            preprocessed::range_check_add::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::range_check_add::{self},
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let range_check_add_cols = range_check_add::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(&scheduling_lookup_data, &compression_lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![
            range_check_add_cols
                .iter()
                .map(|c| circle_evaluation_u32x16!(c))
                .collect::<Vec<_>>(),
            trace
                .into_iter()
                .map(|c| circle_evaluation_u32x16!(c))
                .collect::<Vec<_>>(),
            interaction_trace,
        ]);
        let log_size = (traces[0][0].data.len() * (1 << LOG_N_LANES)).ilog2();

        let trace_polys =
            traces.map(|trace| trace.into_iter().map(|c| c.interpolate()).collect_vec());

        assert_constraints_on_polys(
            &trace_polys,
            CanonicCoset::new(log_size),
            |mut eval| {
                eval_constraints(&mut eval, &relations);
            },
            claimed_sum,
        );
    }
}
