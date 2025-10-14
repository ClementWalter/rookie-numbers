use crate::{add_to_relation, components::scheduling::columns::RoundColumns, relations::Relations};
use num_traits::One;
use stwo_prover::{
    constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval},
    core::fields::m31::M31,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_scheduling_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let w: [E::F; 128] = std::array::from_fn(|_| eval.next_trace_mask());
    let one = E::EF::one();
    for round in 16..64 {
        let RoundColumns {
            w_15_i0_low,
            w_15_i0_high,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
            w_2_i0_low,
            w_2_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high,
            carry_low,
            carry_high,
        } = RoundColumns::<<E as EvalAtRow>::F>::from_eval(eval);

        let w_16_low = &w[(round - 16) * 2];
        let w_16_high = &w[(round - 16) * 2 + 1];
        let w_15_low = &w[(round - 15) * 2];
        let w_15_high = &w[(round - 15) * 2 + 1];
        let w_15_i1_low = w_15_low.clone() - w_15_i0_low.clone();
        let w_15_i1_high = w_15_high.clone() - w_15_i0_high.clone();
        let w_7_low = &w[(round - 7) * 2];
        let w_7_high = &w[(round - 7) * 2 + 1];
        let w_2_low = &w[(round - 2) * 2];
        let w_2_high = &w[(round - 2) * 2 + 1];
        let w_2_i1_low = w_2_low.clone() - w_2_i0_low.clone();
        let w_2_i1_high = w_2_high.clone() - w_2_i0_high.clone();
        let new_w_low = &w[round * 2];
        let new_w_high = &w[round * 2 + 1];

        // SIGMA 0
        add_to_relation!(
            eval,
            relations.sigma_0.i0,
            one,
            w_15_i0_low,
            w_15_i0_high,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_0.i1,
            one,
            w_15_i1_low,
            w_15_i1_high,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_0.o2,
            one,
            sigma_0_o20_pext,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high
        );
        let sigma_0_low = sigma_0_o0_low + sigma_0_o1_low + sigma_0_o2_low;
        let sigma_0_high = sigma_0_o0_high + sigma_0_o1_high + sigma_0_o2_high;

        // SIGMA 1
        add_to_relation!(
            eval,
            relations.sigma_1.i0,
            one,
            w_2_i0_low,
            w_2_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_1.i1,
            one,
            w_2_i1_low,
            w_2_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_1.o2,
            one,
            sigma_1_o20_pext,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high
        );
        let sigma_1_low = sigma_1_o0_low + sigma_1_o1_low + sigma_1_o2_low;
        let sigma_1_high = sigma_1_o0_high + sigma_1_o1_high + sigma_1_o2_high;

        eval.add_constraint(
            new_w_low.clone() + carry_low.clone() * E::F::from(M31::from(1 << 16))
                - w_16_low.clone()
                - sigma_0_low.clone()
                - w_7_low.clone()
                - sigma_1_low.clone(),
        );

        eval.add_constraint(
            new_w_high.clone() + carry_high.clone() * E::F::from(M31::from(1 << 16))
                - w_16_high.clone()
                - sigma_0_high.clone()
                - w_7_high.clone()
                - sigma_1_high.clone()
                - carry_low.clone(),
        );

        // ADD
        add_to_relation!(eval, relations.range_check_add, one, *new_w_low, carry_low);
        add_to_relation!(
            eval,
            relations.range_check_add,
            one,
            *new_w_high,
            carry_high
        );
    }
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
        eval_scheduling_constraints(&mut eval, &self.relations);
        eval
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use stwo_prover::{
        constraint_framework::assert_constraints_on_polys,
        core::{pcs::TreeVec, poly::circle::CanonicCoset},
    };

    use super::*;
    use crate::components::scheduling::witness::{gen_interaction_trace, gen_trace};

    #[test]
    fn test_scheduling_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let (trace, lookup_data) = gen_trace(LOG_N_ROWS);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&lookup_data, &relations);

        let traces = TreeVec::new(vec![vec![], trace, interaction_trace]);
        let trace_polys =
            traces.map(|trace| trace.into_iter().map(|c| c.interpolate()).collect_vec());

        assert_constraints_on_polys(
            &trace_polys,
            CanonicCoset::new(LOG_N_ROWS),
            |mut eval| {
                eval_scheduling_constraints(&mut eval, &relations);
            },
            claimed_sum,
        );
    }
}
