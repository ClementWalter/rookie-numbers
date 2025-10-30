use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};
use utils::add_to_relation;

use crate::{
    components::preprocessed::sigma_0::i0_i1::columns::ComponentColumnsOwned,
    partitions::Sigma0 as Sigma0Partitions, preprocessed::sigma_0::Sigma0I0I1ColumnsOwned,
    relations::Relations,
};

pub type Component = FrameworkComponent<Eval>;

fn eval_constraints<E: EvalAtRow>(eval: &mut E, relations: &Relations, log_size: u32) {
    let chunk_count = 1 << Sigma0Partitions::I0.count_ones().saturating_sub(log_size);
    for chunk in 0..chunk_count {
        let ComponentColumnsOwned { i0_mult, i1_mult } =
            ComponentColumnsOwned::<<E as EvalAtRow>::F>::from_eval(eval);
        let Sigma0I0I1ColumnsOwned {
            i0_low,
            i0_high,
            o0_low,
            o0_high,
            o20_pext,
            i1_low,
            i1_high,
            o1_low,
            o1_high,
            o21_pext,
        } = Sigma0I0I1ColumnsOwned::<<E as EvalAtRow>::F>::from_ids(eval, Some(chunk));
        add_to_relation!(
            eval,
            relations.sigma_0.i0,
            E::EF::from(i0_mult),
            i0_low,
            i0_high,
            o0_low,
            o0_high,
            o20_pext
        );
        add_to_relation!(
            eval,
            relations.sigma_0.i1,
            E::EF::from(i1_mult),
            i1_low,
            i1_high,
            o1_low,
            o1_high,
            o21_pext
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
        Sigma0Partitions::I0.count_ones().min(self.log_size)
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        Sigma0Partitions::I0.count_ones().min(self.log_size) + 1
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
            preprocessed::sigma_0::i0_i1::witness::{gen_interaction_trace, gen_trace},
            scheduling::witness::gen_trace as gen_scheduling_trace,
        },
        preprocessed::sigma_0::{self, Sigma0I0I1Columns as Sigma0I0I1ColumnsBorrowed},
    };

    #[test_log::test]
    fn test_constraints() {
        const LOG_N_ROWS: u32 = 8;

        // Trace.
        let sigma_0_cols = sigma_0::gen_column_simd();

        let (scheduling_trace, scheduling_lookup_data) = gen_scheduling_trace(LOG_N_ROWS);
        let (_, compression_lookup_data) = gen_compression_trace(&scheduling_trace);
        let max_log_size = 21;
        let trace = gen_trace(
            max_log_size,
            &scheduling_lookup_data,
            &compression_lookup_data,
        );

        let simd_size = trace[0].len().ilog2();
        let log_size = simd_size + LOG_N_LANES;

        let relations = Relations::dummy();
        let (interaction_trace, claimed_sum) = gen_interaction_trace(&trace, &relations);

        let sigma_0_i0_i1_cols = &sigma_0_cols[..Sigma0I0I1ColumnsBorrowed::SIZE];
        let preprocessed_trace = Sigma0I0I1ColumnsBorrowed::from_slice(sigma_0_i0_i1_cols)
            .chunks((1 << simd_size) as usize)
            .into_iter()
            .flat_map(|c| c.iter().map(|c| circle_evaluation_u32x16!(c)))
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
