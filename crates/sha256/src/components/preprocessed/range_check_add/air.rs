use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::range_check_add::columns::ComponentColumnsOwned as ComponentColumns,
    preprocessed::range_check_add::RangeCheckAddColumnsOwned as RangeCheckAddColumns,
    relations::Relations,
};
#[cfg(feature = "dynamic-preprocessed-shape")]
use crate::preprocessed_log_size;

pub type Component = FrameworkComponent<Eval>;

#[cfg(feature = "dynamic-preprocessed-shape")]
fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations, log_size: u32) {
    let effective_log_size = preprocessed_log_size(log_size);
    let chunk_count = 1 << 19_u32.saturating_sub(effective_log_size);
    for chunk in 0..chunk_count {
        let ComponentColumns {
            carry_4_mult,
            carry_7_mult,
            carry_8_mult,
        } = ComponentColumns::<<E as EvalAtRow>::F>::from_eval(eval);
        let suffix = if chunk_count == 1 { None } else { Some(chunk) };
        let RangeCheckAddColumns {
            value,
            carry_4,
            carry_7,
            carry_8,
        } = RangeCheckAddColumns::<<E as EvalAtRow>::F>::from_ids(eval, suffix);
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
    }
    eval.finalize_logup_in_pairs();
}

#[cfg(not(feature = "dynamic-preprocessed-shape"))]
fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations, _log_size: u32) {
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
    } = RangeCheckAddColumns::<<E as EvalAtRow>::F>::from_ids(eval, None);

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
    pub log_size: u32,
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        #[cfg(feature = "dynamic-preprocessed-shape")]
        {
            19_u32.min(preprocessed_log_size(self.log_size))
        }
        #[cfg(not(feature = "dynamic-preprocessed-shape"))]
        {
            19
        }
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        #[cfg(feature = "dynamic-preprocessed-shape")]
        {
            19_u32.min(preprocessed_log_size(self.log_size)) + 1
        }
        #[cfg(not(feature = "dynamic-preprocessed-shape"))]
        {
            19 + 1
        }
    }
    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        eval_constraints(&mut eval, &self.relations, self.log_size);
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
        preprocessed::range_check_add::{
            self, RangeCheckAddColumns as RangeCheckAddColumnsBorrowed,
        },
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        #[cfg(feature = "dynamic-preprocessed-shape")]
        let max_log_size = 10;
        #[cfg(feature = "dynamic-preprocessed-shape")]
        let trace = gen_trace(
            max_log_size,
            &scheduling_lookup_data,
            &compression_lookup_data,
        );
        #[cfg(not(feature = "dynamic-preprocessed-shape"))]
        let trace = gen_trace(
            LOG_N_ROWS,
            &scheduling_lookup_data,
            &compression_lookup_data,
        );

        let simd_size = trace[0].len().ilog2();
        let log_size = simd_size + LOG_N_LANES;

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let range_check_add_cols = range_check_add::gen_column_simd();
        #[cfg(feature = "dynamic-preprocessed-shape")]
        let preprocessed_trace = RangeCheckAddColumnsBorrowed::from_slice(&range_check_add_cols)
            .chunks((1 << simd_size) as usize)
            .into_iter()
            .flat_map(|c| c.iter().map(|c| circle_evaluation_u32x16!(c)))
            .collect::<Vec<_>>();
        #[cfg(not(feature = "dynamic-preprocessed-shape"))]
        let preprocessed_trace = range_check_add_cols
            .iter()
            .map(|c| circle_evaluation_u32x16!(c))
            .collect::<Vec<_>>();

        let traces = TreeVec::new(vec![
            preprocessed_trace,
            trace
                .into_iter()
                .map(|c| circle_evaluation_u32x16!(c))
                .collect::<Vec<_>>(),
            interaction_trace,
        ]);

        let trace_polys =
            traces.map(|trace| trace.into_iter().map(|c| c.interpolate()).collect_vec());

        assert_constraints_on_polys(
            &trace_polys,
            CanonicCoset::new(log_size),
            |mut eval| {
                eval_constraints(&mut eval, &relations, log_size);
            },
            claimed_sum,
        );
    }
}
