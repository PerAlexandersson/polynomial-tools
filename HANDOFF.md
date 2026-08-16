# Polytool handoff

## Repository ownership

- `/workspace/rust` branch `master` is the canonical monorepo history.
- All polytool changes are made under `polytool/` on `master`.
- Repository branch `main` is a generated standalone projection of this
  directory. Do not commit to it directly.
- After tested polytool changes have been committed and pushed on `master`, run
  `./scripts/sync-polytool-main.sh` from the monorepo root.

## Current state

The former independently rooted standalone history was backed up before
`main` was replaced by a `git subtree split --prefix=polytool` projection.
The projected tree includes the crate rename to `polytool` and all correctness
fixes from the August 2026 Rust review.
