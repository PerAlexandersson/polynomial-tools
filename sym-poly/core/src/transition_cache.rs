//! Generic transition matrix cache for graded algebras.
//!
//! Each algebra (Sym, QSym, ...) instantiates its own `TransitionCache<B>`
//! with its basis enum type, giving a per-algebra global cache of transition
//! matrices indexed by (source_basis, target_basis, degree).

use std::collections::BTreeMap;
use std::sync::Mutex;

/// A thread-safe cache of transition matrices for a specific algebra.
///
/// `B` is the basis enum type (e.g., `sym::Basis` or `qsym::QSymBasis`).
/// Matrices are cached by `(source, target, degree)` and stored as `Vec<Vec<i64>>`.
///
/// # Usage
///
/// Each algebra crate should create a single `static` instance:
///
/// ```ignore
/// static SYM_CACHE: TransitionCache<Basis> = TransitionCache::new();
/// ```
pub struct TransitionCache<B: Copy + Eq + Ord> {
    cache: Mutex<BTreeMap<(B, B, u32), Vec<Vec<i64>>>>,
}

impl<B: Copy + Eq + Ord> TransitionCache<B> {
    /// Create a new empty cache.
    ///
    /// This is `const` so it can be used in `static` declarations.
    pub const fn new() -> Self {
        TransitionCache {
            cache: Mutex::new(BTreeMap::new()),
        }
    }

    /// Look up or compute a transition matrix.
    ///
    /// If the matrix for `(source, target, degree)` is cached, returns a clone.
    /// Otherwise, calls `compute` to build it, caches the result, and returns it.
    ///
    /// The `compute` closure receives `(source, target, degree)` and should
    /// return the transition matrix.
    pub fn get_or_compute(
        &self,
        source: B,
        target: B,
        degree: u32,
        compute: impl FnOnce(B, B, u32) -> Vec<Vec<i64>>,
    ) -> Vec<Vec<i64>> {
        let key = (source, target, degree);

        // Fast path: check cache
        {
            let guard = self.cache.lock().unwrap();
            if let Some(mat) = guard.get(&key) {
                return mat.clone();
            }
        }

        // Compute outside the lock
        let mat = compute(source, target, degree);

        // Store and return
        {
            let mut guard = self.cache.lock().unwrap();
            guard.insert(key, mat.clone());
        }

        mat
    }

    /// Clear all cached matrices.
    pub fn clear(&self) {
        let mut guard = self.cache.lock().unwrap();
        guard.clear();
    }

    /// Number of cached matrices.
    pub fn len(&self) -> usize {
        let guard = self.cache.lock().unwrap();
        guard.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum TestBasis {
        A,
        B,
    }

    #[test]
    fn test_cache_basic() {
        let cache: TransitionCache<TestBasis> = TransitionCache::new();
        assert!(cache.is_empty());

        let mat = cache.get_or_compute(TestBasis::A, TestBasis::B, 2, |_, _, _| {
            vec![vec![1, 0], vec![0, 1]]
        });
        assert_eq!(mat, vec![vec![1, 0], vec![0, 1]]);
        assert_eq!(cache.len(), 1);

        // Second call should hit cache (compute closure would panic if called)
        let mat2 = cache.get_or_compute(TestBasis::A, TestBasis::B, 2, |_, _, _| {
            panic!("should not be called");
        });
        assert_eq!(mat2, mat);
    }

    #[test]
    fn test_cache_clear() {
        let cache: TransitionCache<TestBasis> = TransitionCache::new();
        cache.get_or_compute(TestBasis::A, TestBasis::B, 1, |_, _, _| vec![vec![1]]);
        assert_eq!(cache.len(), 1);
        cache.clear();
        assert!(cache.is_empty());
    }
}
