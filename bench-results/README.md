# Benchmark Results

This directory is the conventional home for longer `polytool bench` outputs.
Keep short, named Markdown reports here when a run is useful for comparing
solver or search changes.

Recommended recurrence fixture command:

```sh
polytool bench recurrence-fixtures --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/all.md
```

Recommended focused OEIS command:

```sh
polytool bench recurrence-fixtures --only oeis --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/oeis.md
```

The per-run TSV printed to stdout can still be redirected separately when raw
timings are needed. The Markdown report stores category summaries and fixture
summaries, which is usually enough for code review and handoff notes.
