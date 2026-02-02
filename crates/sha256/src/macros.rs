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

        components!(@gen_structs $( ($module, $name) ),+);
        components!(@gen_fns $( ($module, $name) ),+);
    };

    (@gen_structs $( ($module:ident, $name:ident) ),+ ) => {
        #[derive(Clone)]
        pub struct Traces {
            $( pub ${concat($module, _, $name)}: Vec<Vec<u32x16>>, )+
        }

        #[derive(Clone, serde::Serialize, serde::Deserialize)]
        pub struct ClaimedSum {
            $( pub ${concat($module, _, $name)}: SecureField, )+
        }

        impl ClaimedSum {
            pub fn sum(&self) -> SecureField {
                SecureField::zero() $( + self.${concat($module, _, $name)} )+
            }

            pub fn mix_into(&self, channel: &mut impl stwo::core::channel::Channel) {
                channel.mix_felts(&[$( self.${concat($module, _, $name)} ),+]);
            }
        }

        pub struct Components {
            $( pub ${concat($module, _, $name)}: $module::$name::air::Component, )+
        }
    };

    (@gen_fns $( ($module:ident, $name:ident) ),+ ) => {
        pub fn gen_trace(
            log_size: u32,
            scheduling_lookup_data: &[Vec<u32x16>],
            compression_lookup_data: &[Vec<u32x16>],
        ) -> Traces {
            Traces {
                $(
                    ${concat($module, _, $name)}: $module::$name::witness::gen_trace(
                        log_size,
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
                let (trace, ${concat($module, _, $name, _, claimed_sum)}) =
                    $module::$name::witness::gen_interaction_trace(&traces.${concat($module, _, $name)}, relations);
                interaction_trace.extend(trace);
            )+
            (
                interaction_trace,
                ClaimedSum {
                    $( ${concat($module, _, $name)}: ${concat($module, _, $name, _, claimed_sum)}, )+
                },
            )
        }

        impl Traces {
            pub fn len(&self) -> usize {
                0 $( + self.${concat($module, _, $name)}.iter().map(|v| v.len()).sum::<usize>() )+
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
                    $( self.${concat($module, _, $name)} ),+
                ];

                cols
                    .into_iter()
                    .flatten()
                    .map(|col| circle_evaluation_u32x16!(col))
                    .collect::<Vec<_>>()
                    .into_iter()
            }
        }

        impl Components {
            pub fn new(
                log_size: u32,
                location_allocator: &mut TraceLocationAllocator,
                relations: Relations,
                claimed_sum: &ClaimedSum,
            ) -> Self {
                Self {
                    $(
                        ${concat($module, _, $name)}: $module::$name::air::Component::new(
                            location_allocator,
                            $module::$name::air::Eval { log_size, relations: relations.clone() },
                            claimed_sum.${concat($module, _, $name)},
                        ),
                    )+
                }
            }

            pub fn provers(&self) -> Vec<&dyn ComponentProver<SimdBackend>> {
                vec![ $(&self.${concat($module, _, $name)},)+ ]
            }

            pub fn components(&self) -> Vec<&dyn stwo::core::air::Component> {
                vec![ $(&self.${concat($module, _, $name)},)+ ]
            }

            pub fn relation_entries(
                &self,
                trace: &TreeVec<Vec<&Vec<BaseField>>>,
            ) -> Vec<RelationTrackerEntry> {
                itertools::chain!(
                    $( add_to_relation_entries(&self.${concat($module, _, $name)}, trace) ),+
                )
                .collect()
            }

            pub fn trace_log_degree_bounds(&self) -> Vec<TreeVec<ColumnVec<u32>>> {
                vec![
                    $( self.${concat($module, _, $name)}.trace_log_degree_bounds(), )+
                ]
            }

            pub fn n_constraints(&self) -> usize {
                0 $( + self.${concat($module, _, $name)}.n_constraints() )+
            }
        }
    };
}
