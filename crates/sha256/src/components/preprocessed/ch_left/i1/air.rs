use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::ch_left::i1::columns::ComponentColumnsOwned as ComponentColumns,
    partitions::BigSigma1 as BigSigma1Partitions,
    preprocessed::ch_left::ChLeftI1ColumnsOwned as ChLeftI1Columns, relations::Relations,
};

const _: () = assert!(
    BigSigma1Partitions::I1_L.count_ones() == BigSigma1Partitions::I1_H.count_ones(),
    "BigSigma1Partitions::I1_L and I1_H must have the same number of ones"
);

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumns {
        i1_low_mult,
        i1_high_mult,
    } = ComponentColumns::<<E as EvalAtRow>::F>::from_eval(eval);
    let ChLeftI1Columns {
        i1_low_e,
        i1_low_f,
        i1_low_res,
        i1_high_e,
        i1_high_f,
        i1_high_res,
    } = ChLeftI1Columns::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.ch_left.i1_low,
        E::EF::from(i1_low_mult),
        i1_low_e,
        i1_low_f,
        i1_low_res,
    );
    add_to_relation!(
        eval,
        relations.ch_left.i1_high,
        E::EF::from(i1_high_mult),
        i1_high_e,
        i1_high_f,
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
        BigSigma1Partitions::I1.count_ones()
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        BigSigma1Partitions::I1.count_ones() + 1
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
            preprocessed::ch_left::i1::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::ch_left::{self, ChLeftI0Columns},
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let ch_left_cols = ch_left::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(&scheduling_lookup_data, &compression_lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![
            ch_left_cols[ChLeftI0Columns::SIZE..(ChLeftI0Columns::SIZE + ChLeftI1Columns::SIZE)]
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
