use stwo_constraint_framework::{EvalAtRow, FrameworkComponent, FrameworkEval};

pub type DummyComponent = FrameworkComponent<DummyEval>;

#[derive(Clone)]
pub struct DummyEval {
    pub log_size: u32,
    pub n_cols: usize,
}
impl FrameworkEval for DummyEval {
    fn log_size(&self) -> u32 {
        self.log_size
    }
    fn max_constraint_log_degree_bound(&self) -> u32 {
        self.log_size + 1
    }
    fn evaluate<E: EvalAtRow>(&self, mut eval: E) -> E {
        for _ in 0..self.n_cols {
            eval.next_trace_mask();
        }
        eval
    }
}
