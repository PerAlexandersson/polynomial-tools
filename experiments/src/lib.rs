//! Shared helpers for research binaries in `src/bin`.
//!
//! These helpers are intentionally lightweight: they factor out duplicated
//! plumbing used by many exploratory binaries without promoting half-baked
//! APIs into the reusable crates too early.

pub mod peak_utils;
