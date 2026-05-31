# polynomial-lab

`polynomial-lab` is a structured index over project data for polynomial
real-rootedness and interlacing experiments.

The default lab root is:

```text
/workspace/projects/polynomial-interlacing-lab
```

Set `POLY_LAB_ROOT` or pass `--root` to point at another lab directory.

Examples:

```bash
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- validate
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- list-projects
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  trace-goal derangement_descents derangement_descent_real_rootedness
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  render-markdown derangement_descents
```

