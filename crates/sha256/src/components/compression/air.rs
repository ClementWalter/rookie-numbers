use crate::{
    add_to_relation,
    components::compression::columns::RoundColumns,
    relations::Relations,
    sha256::{H, K},
};
use num_traits::One;
use stwo_prover::{
    constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval},
    core::fields::m31::M31,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_compression_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations) {
    let w: [E::F; 128] = std::array::from_fn(|_| eval.next_trace_mask());

    let k: [E::F; K.len() * 2] = K
        .iter()
        .flat_map(|k| {
            [
                E::F::from(M31::from(k & 0xffff)),
                E::F::from(M31::from(k >> 16)),
            ]
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    let mut hash_buffer: [E::F; H.len() * 2] = H
        .iter()
        .flat_map(|h| {
            [
                E::F::from(M31::from(h & 0xffff)),
                E::F::from(M31::from(h >> 16)),
            ]
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    let one = E::EF::one();
    for round in 0..64 {
        let a_low = &hash_buffer[0].clone();
        let a_high = &hash_buffer[1].clone();
        let b_low = &hash_buffer[2].clone();
        let b_high = &hash_buffer[3].clone();
        let c_low = &hash_buffer[4].clone();
        let c_high = &hash_buffer[5].clone();
        let d_low = &hash_buffer[6].clone();
        let d_high = &hash_buffer[7].clone();
        let e_low = &hash_buffer[8].clone();
        let e_high = &hash_buffer[9].clone();
        let f_low = &hash_buffer[10].clone();
        let f_high = &hash_buffer[11].clone();
        let g_low = &hash_buffer[12].clone();
        let g_high = &hash_buffer[13].clone();
        let h_low = &hash_buffer[14].clone();
        let h_high = &hash_buffer[15].clone();

        let k_low = k[2 * round].clone();
        let k_high = k[2 * round + 1].clone();

        let columns = RoundColumns::<<E as EvalAtRow>::F>::from_eval(eval);
        let RoundColumns {
            e_i0_low,
            e_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high,
            f_i0_low,
            f_i0_high,
            ch_left_i0_low,
            ch_left_i0_high,
            ch_left_i1_low,
            ch_left_i1_high,
            g_i0_low,
            g_i0_high,
            ch_right_i0_low,
            ch_right_i0_high,
            ch_right_i1_low,
            ch_right_i1_high,
            a_i0_high_0,
            a_i0_high_1,
            a_i1_low_0,
            a_i1_low_1,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
            b_i0_high_0,
            b_i0_high_1,
            b_i1_low_0,
            b_i1_low_1,
            c_i0_high_0,
            c_i0_high_1,
            c_i1_low_0,
            c_i1_low_1,
            maj_i0_low,
            maj_i0_high_0,
            maj_i0_high_1,
            maj_i1_low_0,
            maj_i1_low_1,
            maj_i1_high,
            e_carry_low,
            e_carry_high,
            a_carry_low,
            a_carry_high,
        } = columns;

        let w_low = &w[round * 2];
        let w_high = &w[round * 2 + 1];
        let a_i0_low =
            a_low.clone() - a_i1_low_0.clone() - a_i1_low_1.clone() * E::F::from(M31::from(1 << 8));
        let a_i1_high = a_high.clone()
            - a_i0_high_0.clone()
            - a_i0_high_1.clone() * E::F::from(M31::from(1 << 8));
        let b_i0_low =
            b_low.clone() - b_i1_low_0.clone() - b_i1_low_1.clone() * E::F::from(M31::from(1 << 8));
        let b_i1_high = b_high.clone()
            - b_i0_high_0.clone()
            - b_i0_high_1.clone() * E::F::from(M31::from(1 << 8));
        let c_i0_low =
            c_low.clone() - c_i1_low_0.clone() - c_i1_low_1.clone() * E::F::from(M31::from(1 << 8));
        let c_i1_high = c_high.clone()
            - c_i0_high_0.clone()
            - c_i0_high_1.clone() * E::F::from(M31::from(1 << 8));
        let e_i1_low = e_low.clone() - e_i0_low.clone();
        let e_i1_high = e_high.clone() - e_i0_high.clone();
        let f_i1_low = f_low.clone() - f_i0_low.clone();
        let f_i1_high = f_high.clone() - f_i0_high.clone();
        let g_i1_low = g_low.clone() - g_i0_low.clone();
        let g_i1_high = g_high.clone() - g_i0_high.clone();

        // BIG_SIGMA1
        add_to_relation!(
            eval,
            relations.big_sigma_1.i0,
            one,
            e_i0_low,
            e_i0_high,
            sigma_1_o0_low,
            sigma_1_o0_high,
            sigma_1_o20_pext
        );
        add_to_relation!(
            eval,
            relations.big_sigma_1.i1,
            one,
            e_i1_low,
            e_i1_high,
            sigma_1_o1_low,
            sigma_1_o1_high,
            sigma_1_o21_pext
        );
        add_to_relation!(
            eval,
            relations.big_sigma_1.o2,
            one,
            sigma_1_o20_pext,
            sigma_1_o21_pext,
            sigma_1_o2_low,
            sigma_1_o2_high
        );
        let sigma_1_low = sigma_1_o0_low + sigma_1_o1_low + sigma_1_o2_low;
        let sigma_1_high = sigma_1_o0_high + sigma_1_o1_high + sigma_1_o2_high;

        // CH_LEFT
        add_to_relation!(
            eval,
            relations.ch_left.i0_low,
            one,
            e_i0_low,
            f_i0_low,
            ch_left_i0_low,
        );
        add_to_relation!(
            eval,
            relations.ch_left.i0_high,
            one,
            e_i0_high,
            f_i0_high,
            ch_left_i0_high,
        );
        add_to_relation!(
            eval,
            relations.ch_left.i1_low,
            one,
            e_i1_low,
            f_i1_low,
            ch_left_i1_low,
        );
        add_to_relation!(
            eval,
            relations.ch_left.i1_high,
            one,
            e_i1_high,
            f_i1_high,
            ch_left_i1_high,
        );
        let ch_left_low = ch_left_i0_low + ch_left_i1_low;
        let ch_left_high = ch_left_i0_high + ch_left_i1_high;

        // CH_RIGHT
        add_to_relation!(
            eval,
            relations.ch_right.i0_low,
            one,
            e_i0_low,
            g_i0_low,
            ch_right_i0_low,
        );
        add_to_relation!(
            eval,
            relations.ch_right.i0_high,
            one,
            e_i0_high,
            g_i0_high,
            ch_right_i0_high,
        );
        add_to_relation!(
            eval,
            relations.ch_right.i1_low,
            one,
            e_i1_low,
            g_i1_low,
            ch_right_i1_low,
        );
        add_to_relation!(
            eval,
            relations.ch_right.i1_high,
            one,
            e_i1_high,
            g_i1_high,
            ch_right_i1_high,
        );
        let ch_right_low = ch_right_i0_low + ch_right_i1_low;
        let ch_right_high = ch_right_i0_high + ch_right_i1_high;

        // BIG SIGMA0
        add_to_relation!(
            eval,
            relations.big_sigma_0.i0,
            one,
            a_i0_low,
            a_i0_high_0,
            a_i0_high_1,
            sigma_0_o0_low,
            sigma_0_o0_high,
            sigma_0_o20_pext,
        );
        add_to_relation!(
            eval,
            relations.big_sigma_0.i1,
            one,
            a_i1_low_0,
            a_i1_low_1,
            a_i1_high,
            sigma_0_o1_low,
            sigma_0_o1_high,
            sigma_0_o21_pext,
        );
        add_to_relation!(
            eval,
            relations.big_sigma_0.o2,
            one,
            sigma_0_o20_pext,
            sigma_0_o21_pext,
            sigma_0_o2_low,
            sigma_0_o2_high,
        );
        let sigma_0_low = sigma_0_o0_low + sigma_0_o1_low + sigma_0_o2_low;
        let sigma_0_high = sigma_0_o0_high + sigma_0_o1_high + sigma_0_o2_high;

        // MAJ
        add_to_relation!(
            eval,
            relations.maj.i0_low,
            one,
            a_i0_low,
            b_i0_low,
            c_i0_low,
            maj_i0_low,
        );
        add_to_relation!(
            eval,
            relations.maj.i0_high_0,
            one,
            a_i0_high_0,
            b_i0_high_0,
            c_i0_high_0,
            maj_i0_high_0,
        );
        add_to_relation!(
            eval,
            relations.maj.i0_high_1,
            one,
            a_i0_high_1,
            b_i0_high_1,
            c_i0_high_1,
            maj_i0_high_1,
        );
        add_to_relation!(
            eval,
            relations.maj.i1_low_0,
            one,
            a_i1_low_0,
            b_i1_low_0,
            c_i1_low_0,
            maj_i1_low_0,
        );
        add_to_relation!(
            eval,
            relations.maj.i1_low_1,
            one,
            a_i1_low_1,
            b_i1_low_1,
            c_i1_low_1,
            maj_i1_low_1,
        );
        add_to_relation!(
            eval,
            relations.maj.i1_high,
            one,
            a_i1_high,
            b_i1_high,
            c_i1_high,
            maj_i1_high,
        );
        let maj_low = maj_i0_low + maj_i1_low_0 + maj_i1_low_1 * E::F::from(M31::from(1 << 8));
        let maj_high = maj_i0_high_0 + maj_i0_high_1 * E::F::from(M31::from(1 << 8)) + maj_i1_high;

        let temp1_low =
            h_low.clone() + sigma_1_low + ch_left_low + ch_right_low + k_low + w_low.clone();
        let temp1_high =
            h_high.clone() + sigma_1_high + ch_left_high + ch_right_high + k_high + w_high.clone();
        let temp2_low = sigma_0_low + maj_low;
        let temp2_high = sigma_0_high + maj_high;
        let new_e_low = d_low.clone() + temp1_low.clone()
            - e_carry_low.clone() * E::F::from(M31::from(1 << 16));
        let new_e_high = d_high.clone() + temp1_high.clone() + e_carry_low.clone()
            - e_carry_high.clone() * E::F::from(M31::from(1 << 16));
        let new_a_low =
            temp1_low + temp2_low - a_carry_low.clone() * E::F::from(M31::from(1 << 16));
        let new_a_high = temp1_high + temp2_high + a_carry_low.clone()
            - a_carry_high.clone() * E::F::from(M31::from(1 << 16));

        // ADD
        add_to_relation!(eval, relations.range_check_add, one, new_e_low, e_carry_low);
        add_to_relation!(
            eval,
            relations.range_check_add,
            one,
            new_e_high,
            e_carry_high,
        );
        add_to_relation!(eval, relations.range_check_add, one, new_a_low, a_carry_low);
        add_to_relation!(
            eval,
            relations.range_check_add,
            one,
            new_a_high,
            a_carry_high,
        );

        hash_buffer[0] = new_a_low; // a_low
        hash_buffer[1] = new_a_high; // a_high
        hash_buffer[2] = a_low.clone(); // b_low
        hash_buffer[3] = a_high.clone(); // b_high
        hash_buffer[4] = b_low.clone(); // c_low
        hash_buffer[5] = b_high.clone(); // c_high
        hash_buffer[6] = c_low.clone(); // d_low
        hash_buffer[7] = c_high.clone(); // d_high
        hash_buffer[8] = new_e_low; // e_low
        hash_buffer[9] = new_e_high; // e_high
        hash_buffer[10] = e_low.clone(); // f_low
        hash_buffer[11] = e_high.clone(); // f_high
        hash_buffer[12] = f_low.clone(); // g_low
        hash_buffer[13] = f_high.clone(); // g_high
        hash_buffer[14] = g_low.clone(); // h_low
        hash_buffer[15] = g_high.clone(); // h_high
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
        eval_compression_constraints(&mut eval, &self.relations);
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
    use crate::components::compression::witness::{gen_interaction_trace, gen_trace};
    use crate::components::scheduling::witness::gen_trace as gen_scheduling_trace;

    #[test]
    fn test_compression_constraints() {
        const LOG_N_ROWS: u32 = 4;

        // Trace.
        let (scheduling_trace, _) = gen_scheduling_trace(LOG_N_ROWS);
        let (trace, lookup_data) = gen_trace(&scheduling_trace);

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&lookup_data, &relations);

        let traces = TreeVec::new(vec![vec![], trace, interaction_trace]);
        let trace_polys =
            traces.map(|trace| trace.into_iter().map(|c| c.interpolate()).collect_vec());

        assert_constraints_on_polys(
            &trace_polys,
            CanonicCoset::new(LOG_N_ROWS),
            |mut eval| {
                eval_compression_constraints(&mut eval, &relations);
            },
            claimed_sum,
        );
    }
}
