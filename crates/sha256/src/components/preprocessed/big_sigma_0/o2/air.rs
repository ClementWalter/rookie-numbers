use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::sigma_0::o2::columns::ComponentColumnsOwned as ComponentColumns,
    partitions::BigSigma0 as BigSigma0Partitions,
    preprocessed::big_sigma_0::BigSigma0O2ColumnsOwned as BigSigma0O2Columns, relations::Relations,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumns { o2_mult } = ComponentColumns::<<E as EvalAtRow>::F>::from_eval(eval);

    let BigSigma0O2Columns {
        o2_0,
        o2_1,
        o2_low,
        o2_high,
    } = BigSigma0O2Columns::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.big_sigma_0.o2,
        E::EF::from(o2_mult),
        o2_0,
        o2_1,
        o2_low,
        o2_high,
    );
    eval.finalize_logup_in_pairs();
}

#[derive(Clone)]
pub struct Eval {
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        BigSigma0Partitions::O2.count_ones() * 2
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        BigSigma0Partitions::O2.count_ones() * 2 + 1
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
            preprocessed::big_sigma_0::o2::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::big_sigma_0::{self, BigSigma0I0I1Columns},
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let big_sigma_0_cols = big_sigma_0::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(&scheduling_lookup_data, &compression_lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![
            big_sigma_0_cols[BigSigma0I0I1Columns::SIZE
                ..(BigSigma0I0I1Columns::SIZE + BigSigma0O2Columns::SIZE)]
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
