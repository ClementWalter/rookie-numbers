use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::big_sigma_1::i0::columns::ComponentColumnsOwned,
    partitions::BigSigma1 as BigSigma1Partitions,
    preprocessed::big_sigma_1::BigSigma1I0ColumnsOwned, relations::Relations,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumnsOwned { i0_mult } =
        ComponentColumnsOwned::<<E as EvalAtRow>::F>::from_eval(eval);
    let BigSigma1I0ColumnsOwned {
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext,
    } = BigSigma1I0ColumnsOwned::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.big_sigma_1.i0,
        E::EF::from(i0_mult),
        i0_low,
        i0_high,
        o0_low,
        o0_high,
        o20_pext
    );
    eval.finalize_logup_in_pairs();
}

#[derive(Clone)]
pub struct Eval {
    pub relations: Relations,
}
impl FrameworkEval for Eval {
    fn log_size(&self) -> u32 {
        BigSigma1Partitions::I0.count_ones()
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        BigSigma1Partitions::I0.count_ones() + 1
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
            preprocessed::big_sigma_1::i0::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::big_sigma_1::{self, BigSigma1I0ColumnsOwned},
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let big_sigma_1_cols = big_sigma_1::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let trace = gen_trace(&scheduling_lookup_data, &compression_lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![
            big_sigma_1_cols[..BigSigma1I0ColumnsOwned::SIZE]
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
