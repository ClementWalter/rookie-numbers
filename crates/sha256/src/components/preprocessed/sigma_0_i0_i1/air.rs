use stwo_constraint_framework::{EvalAtRow, FrameworkEval};

use crate::{
    add_to_relation, components::preprocessed::sigma_0_i0_i1::columns::ComponentColumnsOwned,
    preprocessed::sigma_0::ColumnsOwned as Sigma0ColumnsOwned, relations::Relations,
};

// pub type Component = FrameworkComponent<Eval>;

fn eval_sigma_0_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let ComponentColumnsOwned {
        sigma_0_i0_mult,
        sigma_0_i1_mult,
    } = ComponentColumnsOwned::<<E as EvalAtRow>::F>::from_eval(eval);
    let Sigma0ColumnsOwned {
        sigma_0_i0_low,
        sigma_0_i0_high,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext,
        sigma_0_i1_low,
        sigma_0_i1_high,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext,
        sigma_0_o2_0: _,
        sigma_0_o2_1: _,
        sigma_0_o2_low: _,
        sigma_0_o2_high: _,
    } = Sigma0ColumnsOwned::<<E as EvalAtRow>::F>::from_ids(eval);

    add_to_relation!(
        eval,
        relations.sigma_0.i0,
        E::EF::from(sigma_0_i0_mult),
        sigma_0_i0_low,
        sigma_0_i0_high,
        sigma_0_o0_low,
        sigma_0_o0_high,
        sigma_0_o20_pext
    );
    add_to_relation!(
        eval,
        relations.sigma_0.i1,
        E::EF::from(sigma_0_i1_mult),
        sigma_0_i1_low,
        sigma_0_i1_high,
        sigma_0_o1_low,
        sigma_0_o1_high,
        sigma_0_o21_pext
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
        self.log_size
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        self.log_size() + 1
    }
    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        eval_sigma_0_constraints(&mut eval, &self.relations);
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

    use super::*;
    use crate::{
        components::{
            preprocessed::sigma_0_i0_i1::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::sigma_0,
    };

    #[cfg_attr(not(feature = "slow-tests"), ignore)]
    #[test_log::test]
    fn test_sigma_0_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let mut preprocessed_trace = sigma_0::gen_column_simd()[0..10].to_vec();

        // from_ids requires enough columns, we just extend to match the number of columns required by the relations
        preprocessed_trace.extend_from_slice(&sigma_0::gen_column_simd()[0..4]);

        let (_, lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let trace = gen_trace(&lookup_data);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let traces = TreeVec::new(vec![preprocessed_trace, trace, interaction_trace]);
        let log_size = (traces[0][0].data.len() * (1 << LOG_N_LANES)).ilog2();

        let trace_polys =
            traces.map(|trace| trace.into_iter().map(|c| c.interpolate()).collect_vec());

        assert_constraints_on_polys(
            &trace_polys,
            CanonicCoset::new(log_size),
            |mut eval| {
                eval_sigma_0_constraints(&mut eval, &relations);
            },
            claimed_sum,
        );
    }
}
