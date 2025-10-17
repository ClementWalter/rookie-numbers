use num_traits::One;
use stwo::core::fields::m31::M31;
use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};

use crate::{
    add_to_relation,
    components::{scheduling::columns::RoundColumnsOwned, W_SIZE},
    relations::Relations,
    sha256::N_SCHEDULING_ROUNDS,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_scheduling_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let w: [E::F; W_SIZE] = std::array::from_fn(|_| eval.next_trace_mask());
    let one = E::EF::one();
    let minus_one = -E::EF::one();
    for t in 16..(16 + N_SCHEDULING_ROUNDS) {
        let cols = RoundColumnsOwned::<<E as EvalAtRow>::F>::from_eval(eval);

        let w_16_low = &w[(t - 16) * 2];
        let w_16_high = &w[(t - 16) * 2 + 1];
        let w_15_low = &w[(t - 15) * 2];
        let w_15_high = &w[(t - 15) * 2 + 1];
        let w_15_i1_low = w_15_low.clone() - cols.w_15_i0_low.clone();
        let w_15_i1_high = w_15_high.clone() - cols.w_15_i0_high.clone();
        let w_7_low = &w[(t - 7) * 2];
        let w_7_high = &w[(t - 7) * 2 + 1];
        let w_2_low = &w[(t - 2) * 2];
        let w_2_high = &w[(t - 2) * 2 + 1];
        let w_2_i1_low = w_2_low.clone() - cols.w_2_i0_low.clone();
        let w_2_i1_high = w_2_high.clone() - cols.w_2_i0_high.clone();
        let new_w_low = &w[t * 2];
        let new_w_high = &w[t * 2 + 1];

        let sigma_0_low =
            cols.sigma_0_o0_low.clone() + cols.sigma_0_o1_low.clone() + cols.sigma_0_o2_low.clone();
        let sigma_0_high = cols.sigma_0_o0_high.clone()
            + cols.sigma_0_o1_high.clone()
            + cols.sigma_0_o2_high.clone();
        let sigma_1_low =
            cols.sigma_1_o0_low.clone() + cols.sigma_1_o1_low.clone() + cols.sigma_1_o2_low.clone();
        let sigma_1_high = cols.sigma_1_o0_high.clone()
            + cols.sigma_1_o1_high.clone()
            + cols.sigma_1_o2_high.clone();

        eval.add_constraint(
            new_w_low.clone() + cols.carry_low.clone() * E::F::from(M31::from(1 << 16))
                - w_16_low.clone()
                - sigma_0_low.clone()
                - w_7_low.clone()
                - sigma_1_low.clone(),
        );

        eval.add_constraint(
            new_w_high.clone() + cols.carry_high.clone() * E::F::from(M31::from(1 << 16))
                - w_16_high.clone()
                - sigma_0_high.clone()
                - w_7_high.clone()
                - sigma_1_high.clone()
                - cols.carry_low.clone(),
        );

        // SIGMA 0
        add_to_relation!(
            eval,
            relations.sigma_0.i0,
            minus_one,
            cols.w_15_i0_low,
            cols.w_15_i0_high,
            cols.sigma_0_o0_low,
            cols.sigma_0_o0_high,
            cols.sigma_0_o20_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_0.i1,
            minus_one,
            w_15_i1_low,
            w_15_i1_high,
            cols.sigma_0_o1_low,
            cols.sigma_0_o1_high,
            cols.sigma_0_o21_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_0.o2,
            minus_one,
            cols.sigma_0_o20_pext,
            cols.sigma_0_o21_pext,
            cols.sigma_0_o2_low,
            cols.sigma_0_o2_high
        );

        // SIGMA 1
        add_to_relation!(
            eval,
            relations.sigma_1.i0,
            minus_one,
            cols.w_2_i0_low,
            cols.w_2_i0_high,
            cols.sigma_1_o0_low,
            cols.sigma_1_o0_high,
            cols.sigma_1_o20_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_1.i1,
            minus_one,
            w_2_i1_low,
            w_2_i1_high,
            cols.sigma_1_o1_low,
            cols.sigma_1_o1_high,
            cols.sigma_1_o21_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_1.o2,
            minus_one,
            cols.sigma_1_o20_pext,
            cols.sigma_1_o21_pext,
            cols.sigma_1_o2_low,
            cols.sigma_1_o2_high
        );

        // ADD
        add_to_relation!(
            eval,
            relations.range_check_add.add_4,
            minus_one,
            *new_w_low,
            cols.carry_low
        );
        add_to_relation!(
            eval,
            relations.range_check_add.add_4,
            minus_one,
            *new_w_high,
            cols.carry_high
        );
    }

    // Emit W consumed by compression
    eval.add_to_relation(stwo_constraint_framework::RelationEntry::new(
        &relations.w,
        one,
        &w,
    ));

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
    use stwo::core::{pcs::TreeVec, poly::circle::CanonicCoset};
    use stwo_constraint_framework::assert_constraints_on_polys;

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
