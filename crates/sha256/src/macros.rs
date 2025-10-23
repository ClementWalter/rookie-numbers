#[macro_export]
macro_rules! components {
    // handle foo::bar forms
    ( $( $module:ident :: $name:ident ),+ $(,)? ) => {
        $crate::components!(@impl $( ($module, $name) ),+ );
    };
    // internal implementation
    (@impl $( ($module:ident, $name:ident) ),+ $(,)? ) => {
        use $crate::relations::Relations;
        use utils::circle_evaluation_u32x16;

        use stwo::core::fields::qm31::SecureField;
        use std::simd::u32x16;
        use num_traits::Zero;
        use stwo::{
            core::{fields::{m31::BaseField}, ColumnVec},
            prover::{backend::simd::{SimdBackend}, poly::{circle::CircleEvaluation, BitReversedOrder}, ComponentProver},
        };
        use stwo::core::pcs::TreeVec;
        use stwo::core::air::Component;
        use stwo_constraint_framework::{
            relation_tracker::{add_to_relation_entries, RelationTrackerEntry},
            TraceLocationAllocator,
        };

        paste::paste! {
            #[derive(Clone)]
            pub struct Traces {
                $( pub [<$module _ $name>]: Vec<Vec<u32x16>>, )+
            }

            pub struct ClaimedSum {
                $( pub [<$module _ $name>]: SecureField, )+
            }

            impl ClaimedSum {
                pub fn sum(&self) -> SecureField {
                    SecureField::zero() $( + self.[<$module _ $name>] )+
                }
            }

            pub fn gen_trace(
                scheduling_lookup_data: &[Vec<u32x16>],
                compression_lookup_data: &[Vec<u32x16>],
            ) -> Traces {
                Traces {
                    $(
                        [<$module _ $name>]: $module::$name::witness::gen_trace(
                            scheduling_lookup_data,
                            compression_lookup_data
                        ),
                    )+
                }
            }

            pub fn gen_interaction_trace(
                traces: &Traces,
                relations: &Relations,
            ) -> (
                ColumnVec<CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>>,
                ClaimedSum,
            ) {
                let mut interaction_trace = vec![];

                $(
                    let (trace, [<$module _ $name _ claimed_sum>]) =
                        $module::$name::witness::gen_interaction_trace(&traces.[<$module _ $name>], relations);
                    interaction_trace.extend(trace);
                )+
                (
                    interaction_trace,
                    ClaimedSum {
                        $( [<$module _ $name>]: [<$module _ $name _ claimed_sum>], )+
                    },
                )
            }

            impl Traces {
                pub fn len(&self) -> usize {
                    0 $( + self.[<$module _ $name>].iter().map(|v| v.len()).sum::<usize>() )+
                }

                pub fn is_empty(&self) -> bool {
                    self.len() == 0
                }
            }

            impl IntoIterator for Traces {
                type Item = CircleEvaluation<SimdBackend, BaseField, BitReversedOrder>;
                type IntoIter = std::vec::IntoIter<Self::Item>;

                fn into_iter(self) -> Self::IntoIter {
                    let cols = vec![
                        $( self.[<$module _ $name>] ),+
                    ];

                    cols
                        .into_iter()
                        .flatten()
                        .map(|col| circle_evaluation_u32x16!(col))
                        .collect::<Vec<_>>()
                        .into_iter()
                }
            }

            pub struct Components {
                $( pub [<$module _ $name>]: $module::$name::air::Component, )+
            }

            impl Components {
                pub fn new(
                    location_allocator: &mut TraceLocationAllocator,
                    relations: Relations,
                    claimed_sum: &ClaimedSum,
                ) -> Self {
                    Self {
                        $(
                            [<$module _ $name>]: $module::$name::air::Component::new(
                                location_allocator,
                                $module::$name::air::Eval { relations: relations.clone() },
                                claimed_sum.[<$module _ $name>],
                            ),
                        )+
                    }
                }

                pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
                    vec![ $(&self.[<$module _ $name>],)+ ]
                }

                pub fn relation_entries(
                    &self,
                    trace: &TreeVec<Vec<&Vec<BaseField>>>,
                ) -> Vec<RelationTrackerEntry> {
                    itertools::chain!(
                        $( add_to_relation_entries(&self.[<$module _ $name>], trace) ),+
                    )
                    .collect()
                }

                pub fn trace_log_degree_bounds(&self) -> Vec<TreeVec<ColumnVec<u32>>> {
                    vec![
                        $( self.[<$module _ $name>].trace_log_degree_bounds(), )+
                    ]
                }
            }
        }
    };
}
