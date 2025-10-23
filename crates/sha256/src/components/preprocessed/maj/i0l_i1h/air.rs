use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::maj::i0l_i1h::columns::ComponentColumnsOwned as ComponentColumns,
    partitions::BigSigma0 as BigSigma0Partitions,
    preprocessed::maj::MajI0LI1HColumnsOwned as MajI0LI1HColumns, relations::Relations,
};

const _: () = assert!(
    BigSigma0Partitions::I0_L.count_ones() == BigSigma0Partitions::I1_H.count_ones(),
    "BigSigma0Partitions::I0_L and I1_H must have the same number of ones"
);

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumns {
        i0_low_mult,
        i1_high_mult,
    } = ComponentColumns::<<E as EvalAtRow>::F>::from_eval(eval);
    let MajI0LI1HColumns {
        i0_low_a,
        i0_low_b,
        i0_low_c,
        i0_low_res,
        i1_high_a,
        i1_high_b,
        i1_high_c,
        i1_high_res,
    } = MajI0LI1HColumns::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.maj.i0_low,
        E::EF::from(i0_low_mult),
        i0_low_a,
        i0_low_b,
        i0_low_c,
        i0_low_res,
    );
    add_to_relation!(
        eval,
        relations.maj.i1_high,
        E::EF::from(i1_high_mult),
        i1_high_a,
        i1_high_b,
        i1_high_c,
        i1_high_res,
    );
    eval.finalize_logup_in_pairs();
}

#[derive(Clone)]
pub struct Eval {
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        BigSigma0Partitions::I0_L.count_ones() * 3
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        BigSigma0Partitions::I0_L.count_ones() * 3 + 1
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
            preprocessed::maj::i0l_i1h::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::maj,
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let maj_cols = maj::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(&scheduling_lookup_data, &compression_lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![
            maj_cols[..8]
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
