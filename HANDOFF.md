# Handoff

## SymCat weighted-bond example

The main website worker owns only
`sym-poly/sym/examples/weighted_bond_symmetric_site_example.rs`.  It implements
the weighted-bond recurrence locally to verify the (P_3) and (K_3) examples
used on `symmetricfunctions.com`; it does not change the public Rust API.

Other modified Rust files predate this task and remain owned by their existing
workers.  In particular, this work does not edit `sym-poly/sym/src/chromatic.rs`
or `sym-poly/sym/src/lib.rs`.

Verification:

```text
timeout 60s nice -n 10 cargo run -q -p sym-poly-sym \
  --example weighted_bond_symmetric_site_example
```
