use polynomial_lab::{LabStore, DEFAULT_LAB_ROOT};
use std::path::Path;

#[test]
fn loads_real_derangement_lab_when_available() {
    if !Path::new(DEFAULT_LAB_ROOT).exists() {
        return;
    }

    let store = LabStore::load(DEFAULT_LAB_ROOT).expect("real lab root should load");
    let report = store.project_report("derangement_descents");
    assert!(report.project.is_some());
    assert!(report
        .goals
        .iter()
        .any(|record| record.id == "derangement_descent_real_rootedness"));
}
