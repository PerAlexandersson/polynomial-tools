//! Tools for the flagged skew Schur Lorentzian project.
//!
//! The crate currently focuses on reusable enumeration and fiber-count checks
//! for candidate statistics in the \(2\times2\) defect-switch problem.

pub mod crystal;
pub mod defect_columns;
pub mod descent;
pub mod enumeration;
pub mod families;
pub mod fiber;
pub mod gt;
pub mod shape;

pub use crystal::{
    active_component_crystal_e_images, active_component_crystal_f_images,
    active_component_crystal_images, active_components, apply_crystal_operator,
    apply_crystal_operator_on_order, bender_knuth_e_images, bender_knuth_f_images,
    bender_knuth_unit_images, crystal_e, crystal_f, is_semistandard, CrystalDirection,
};
pub use defect_columns::{
    check_two_by_two_defect_column_count_fibers, check_two_by_two_defect_column_fibers,
    check_two_by_two_pair_defect_column_count_fibers, check_two_by_two_pair_defect_column_fibers,
    columns_containing_both_defects, number_of_columns_containing_both_defects,
    pair_columns_containing_both_defects, pair_number_of_columns_containing_both_defects,
    DefectColumnCountData, DefectColumnCountFailure, DefectColumnCountScanReport, DefectColumnData,
    DefectColumnFailure, DefectColumnScanReport, PairDefectColumnCountFailure,
    PairDefectColumnCountScanReport, PairDefectColumnFailure, PairDefectColumnScanReport,
};
pub use descent::{
    active_subword_descent_data_for_values, active_subword_descent_mask_for_order,
    descent_data_for_values, descent_mask_for_order, DescentData, DescentStatistic,
};
pub use enumeration::{
    enumerate_content_statistic_counts, enumerate_tableaux, EnumerationLimitExceeded,
    EnumerationOptions, TableauRecord,
};
pub use fiber::{
    check_two_by_two_fiber_at, check_two_by_two_fiber_inequalities, FiberFailure, FiberScanReport,
};
pub use gt::{
    active_gt_row, add_patterns, add_rows, elementary_row_exchange_neighbors, inner_row,
    is_gt_array, pair_envelope, sharp_flag, subtract_pattern_sums, subtract_rows, GtRow,
    SkewGtPattern,
};
pub use shape::{Cell, RowFlaggedSkewShape, SkewShape};
