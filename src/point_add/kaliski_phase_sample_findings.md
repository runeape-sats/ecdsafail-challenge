# Sample hypothesis finding

A direct targeted test was run for the strict failing case `k = 4`:
- generate the accepted 9024-point sample and post-generation RNG state from
  the **experimental** circuit's Fiat-Shamir seed,
- then run both the experimental and generic full circuits on that exact same
  sample / RNG prefix.

## Result
- experimental circuit: first phase failure at batch 10 with mask
  `0x0000040000000000`
- generic circuit on that same experimental sample / RNG prefix: **no phase
  failure** in the first 16 batches tested

## Interpretation
This rules out the strongest alternative hypothesis that the weird strict
failures were mainly caused by a bad luck point sample.

So the phase bug is genuinely introduced by the experimental circuit itself,
not merely revealed by a different sample of elliptic-curve points.
